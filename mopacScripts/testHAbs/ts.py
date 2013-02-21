import os

import logging
import cclib.parser
import openbabel
import numpy
from subprocess import Popen
from collections import defaultdict

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.molecule import Geometry
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe

import rdkit

# script to prep ts structures
actions = [
            ['BREAK_BOND', '*1', 'S', '*2'],
            ['FORM_BOND', '*2', 'S', '*3'],
            ['GAIN_RADICAL', '*1', '1'],
            ['LOSE_RADICAL', '*3', '1']
            ]

family = 'H_Abstraction'

reactRecipe = ReactionRecipe(actions)

template = KineticsFamily(forwardRecipe=reactRecipe)

trusted = open('/Users/pierreb/Code/RMG-database/input/kinetics/families/H_Abstraction/training.py')

lines = trusted.readlines()
k = 0
idx = 0
reactants1 = defaultdict(list)
reactants2 = defaultdict(list)
for line in lines:
    k += 1
    if line.startswith('    reactant1 ='):
        idx += 1
        for num in range(k, k+100):
            if lines[num].find('{') != -1 or lines[num].find('*') != -1:
                reactants1[idx].append(lines[num])
            elif lines[num].find(',') != -1:
                break
    elif line.startswith('    reactant2 ='):
        for num in range(k, k+100):
            if lines[num].find('{') != -1 or lines[num].find('*') != -1:
                reactants2[idx].append(lines[num])
            elif lines[num].find(',') != -1:
                break

tsStructures = list()
for idx in range(1, len(reactants1) + 1):
    r1 = ''
    r2 = ''
    for line in reactants1[idx]:
        r1 = r1 + line
    for line in reactants2[idx]:
        r2 = r2 + line
    r1 = Molecule().fromAdjacencyList(r1)
    r2 = Molecule().fromAdjacencyList(r2)
    rStruct = [r1, r2]
    pStruct, tsStruct = template.applyRecipe(rStruct, getTS=True)
    tsStructures.append(tsStruct)

########################################################################################    
inputFileExtension = '.mop'
outputFileExtension = '.out'
executablePath = os.path.join(os.getenv('MOPAC_DIR', default="/opt/mopac") , 'MOPAC2012.exe')
attempt = 1

usePolar = False

keywords = [
            {'top':"precise nosym", 'bottom':"oldgeo thermo nosym precise "},
            {'top':"precise nosym gnorm=0.0 nonr", 'bottom':"oldgeo thermo nosym precise "},
            {'top':"precise nosym gnorm=0.0", 'bottom':"oldgeo thermo nosym precise "},
            {'top':"precise nosym gnorm=0.0 bfgs", 'bottom':"oldgeo thermo nosym precise "},
            {'top':"precise nosym recalc=10 dmax=0.10 nonr cycles=2000 t=2000", 'bottom':"oldgeo thermo nosym precise "},
            ]

multiplicityKeywords = {
                         1: '',
                         2: 'uhf doublet',
                         3: 'uhf triplet',
                         4: 'uhf quartet',
                         5: 'uhf quintet',
                         6: 'uhf sextet',
                         7: 'uhf septet',
                         8: 'uhf octet',
                         9: 'uhf nonet',
                        }

@property
def scriptAttempts():
    "The number of attempts with different script keywords"
    return len(keywords)

@property
def maxAttempts():
    "The total number of attempts to try"
    return 2 * len(keywords)

def inputFileKeywords(attempt):
    """
    Return the top keywords for attempt number `attempt`.

    NB. `attempt`s begin at 1, not 0.
    """
    assert attempt <= maxAttempts
    if attempt > scriptAttempts:
        attempt -= scriptAttempts
    return keywords[attempt-1]

def run(executablePath, inputFilePath, outputFilePath):
    # submits the input file to mopac
    process = Popen([executablePath, inputFilePath, outputFilePath])
    process.communicate()# necessary to wait for executable termination!

def fixSortLabel(molecule):
    """
    This may not be required anymore. Was needed as when molecules were created, the
    rmg sorting labels would be set after where we tried to generate the TS.
    """
    sortLbl = 0
    for vertex in molecule.vertices:
        vertex.sortingLabel = sortLbl
        sortLbl += 1
    return molecule

def getGeometry(molecule, settings):

    multiplicity = sum([i.radicalElectrons for i in molecule.atoms]) + 1
    geom = Geometry(settings, molecule.toAugmentedInChIKey(), molecule, multiplicity)

    return geom, multiplicity

def getRDKitMol(geometry):
    """
    Check there is no RDKit mol file already made. If so, use rdkit to make a rdmol from
    a mol file. If not, make rdmol from geometry.
    """ 
    geometry.generateRDKitGeometries()
    rdKitMol = rdkit.Chem.MolFromMolFile(geometry.getCrudeMolFilePath(), removeHs=False)      

    return rdKitMol

def generateBoundsMatrix(molecule, settings):
    """
    Uses rdkit to generate the bounds matrix of a rdkit molecule.
    """
    geometry, multiplicity = getGeometry(molecule, settings)
    rdKitMol = getRDKitMol(geometry)
    boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)

    return rdKitMol, boundsMatrix, multiplicity, geometry

def writeInputFile(inputFilePath, molFilePathForCalc, geometry):
    """
    Using the :class:`Geometry` object, write the input file
    for the `attmept`th attempt.
    """

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mop")
    mol = openbabel.OBMol()
    
    obConversion.ReadFile(mol, molFilePathForCalc )
    
    mol.SetTitle(geometry.uniqueIDlong)
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    bottom_keys = 'oldgeo force'
    with open(inputFilePath, 'w') as mopacFile:
        # mopacFile.write(top_keys)
        mopacFile.write(input_string)
        mopacFile.write('\n')
        mopacFile.write(bottom_keys)
        mopacFile.write('\n')
        # mopacFile.write(bottom_keys)
        # if usePolar:
        #     mopacFile.write('\n\n\n')
        #     mopacFile.write(polar_keys)

def inputFileKeywords(attempt):
    """
    Return the top, bottom, and polar keywords for attempt number `attempt`.
    
    NB. `attempt`s begin at 1, not 0.
    """
    assert attempt <= maxAttempts
    
    if attempt > scriptAttempts:
        attempt -= scriptAttempts
    
    multiplicity_keys = multiplicityKeywords[geometry.multiplicity]

    top_keys = "pm7 {0} {1}".format(
            multiplicity_keys,
            keywords[attempt-1]['top'],
            )
    bottom_keys = "{0} pm7 {1}".format(
            keywords[attempt-1]['bottom'],
            multiplicity_keys,
            )
    polar_keys = "oldgeo {0} nosym precise pm7 {1}".format(
            'polar' if geometry.multiplicity == 1 else 'static',
            multiplicity_keys,
            )

    return top_keys, bottom_keys, polar_keys

def writeReferenceFile(inputFilePath, molFilePathForCalc, geometry, attempt, outputFile=None):
    """
    Using the :class:`Geometry` object, write the input file
    for the `attmept`th attempt.
    """
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mop")
    if outputFile:
        geomLines = parseArc(outputFile)
        atomcoords = []
        atomnos = []
        for line in geomLines:
            if len(line.split('+')) > 1:
                parts = line.split('+')
                x = parts[0]
                y = parts[1]
                z = parts[2]
                atType = x.split()[0]
                if atType == 'H':
                    atNum = 1
                elif atType == 'C':
                    atNum = 6
                elif atType == 'O':
                    atNum = 8
                coords = [float(x.split()[1]),float(y.split()[1]),float(z.split()[1])]
                atomnos.append(atNum)
                atomcoords.append(coords)
        atomnos = numpy.array(atomnos, dtype=int)
        atomcoords = numpy.array(atomcoords)
        mol = cclib.bridge.makeopenbabel(atomcoords, atomnos)
    else:
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, molFilePathForCalc )
    
    mol.SetTitle(geometry.uniqueIDlong)    
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    with open(inputFilePath, 'w') as mopacFile:
        mopacFile.write(input_string)
        mopacFile.write('\n')
        
def writeSaddleInputFile(inputFilePath, reactantRefPath, productRefPath, geometryR, geometryP):
    """
    Using the :class:`Geometry` object, write the input file
    for the `attmept`th attempt.
    """
    # reactant
    geomLines = parseArc(reactantRefPath)
    atomcoords = []
    atomnos = []
    for line in geomLines:
        if len(line.split('+')) > 1:
            parts = line.split('+')
            x = parts[0]
            y = parts[1]
            z = parts[2]
            atType = x.split()[0]
            if atType == 'H':
                atNum = 1
            elif atType == 'C':
                atNum = 6
            elif atType == 'O':
                atNum = 8
            coords = [float(x.split()[1]),float(y.split()[1]),float(z.split()[1])]
            atomnos.append(atNum)
            atomcoords.append(coords)
    atomnos = numpy.array(atomnos, dtype=int)
    atomcoords = numpy.array(atomcoords)
    molR = cclib.bridge.makeopenbabel(atomcoords, atomnos)
    
    # product
    geomLines = parseArc(productRefPath)
    atomcoords = []
    atomnos = []
    for line in geomLines:
        if len(line.split('+')) > 1:
            parts = line.split('+')
            x = parts[0]
            y = parts[1]
            z = parts[2]
            atType = x.split()[0]
            if atType == 'H':
                atNum = 1
            elif atType == 'C':
                atNum = 6
            elif atType == 'O':
                atNum = 8
            coords = [float(x.split()[1]),float(y.split()[1]),float(z.split()[1])]
            atomnos.append(atNum)
            atomcoords.append(coords)
    atomnos = numpy.array(atomnos, dtype=int)
    atomcoords = numpy.array(atomcoords)
    molP = cclib.bridge.makeopenbabel(atomcoords, atomnos)
    # 
    # refR = reactantRefPath.split('.')[0] + '.new'
    # rInput = open(refR, 'r')
    # rInput.readline()
    # rString = rInput.read()
    # 
    # refP = productRefPath.split('.')[0] + '.new'
    # pInput = open(refP, 'r')
    # pInput.readline()
    # pString = pInput.read()
    
    obrConversion = openbabel.OBConversion()
    obrConversion.SetInAndOutFormats("mol", "mop")
    obpConversion = openbabel.OBConversion()
    obpConversion.SetInAndOutFormats("mol", "mop")
    
    molR.SetTitle(geometryR.uniqueIDlong)
    obrConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    obpConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    
    input_stringR = obrConversion.WriteString(molR)
    input_stringP = obpConversion.WriteString(molP)
    
    input_string = input_stringR + input_stringP[2:]
    top_keys = 'saddle'
    
    with open(inputFilePath, 'w') as mopacFile:
        mopacFile.write(top_keys)
        mopacFile.write(input_string)
        mopacFile.write('\n')

def parseArc(outputFile):
    arcFile = outputFile.split('.')[0] + '.arc'
    geomLines = list()
    readFile = file(arcFile)
    for line in reversed(readFile.readlines()):
        if line.startswith('geo_ref'):
            break
        else:
            geomLines.append(line)
    
    geomLines.reverse()
    
    return geomLines

def writeGeoRefInputFile(inputFilePath, molFilePathForCalc, refFilePath, geometry):
    refFile = 'geo_ref="' + refFilePath + '"'
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mop")
    
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, molFilePathForCalc )
    mol.SetTitle(geometry.uniqueIDlong)
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    
    with open(inputFilePath, 'w') as mopacFile:
        mopacFile.write(refFile)
        mopacFile.write(input_string)
        mopacFile.write('\n')

def writeTSInputFile(inputFilePath, saddleOutput, count):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mop")
    parseOutput = cclib.parser.Mopac(saddleOutput)
    parseOutput = parseOutput.parse()
    mol = cclib.bridge.makeopenbabel(parseOutput.atomcoords[0], parseOutput.atomnos)
    mol.SetTitle('transitionState' + str(count))
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    
    top_keys = 'ts'
    bottom_keys = 'oldgeo force'
    with open(inputFilePath, 'w') as mopacFile:
        mopacFile.write(top_keys)
        mopacFile.write(input_string)
        mopacFile.write('\n')
        mopacFile.write(bottom_keys)
        mopacFile.write('\n')

def writeIRC(inputFilePath, tsOutPath, count):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mop")
    parseOutput = cclib.parser.Mopac(tsOutPath)
    parseOutput = parseOutput.parse()
    mol = cclib.bridge.makeopenbabel(parseOutput.atomcoords[0], parseOutput.atomnos)
    mol.SetTitle('transitionState' + str(count))
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    
    top_keys = 'irc=1*'
    with open(inputFilePath, 'w') as mopacFile:
        mopacFile.write(top_keys)
        mopacFile.write(input_string)
        mopacFile.write('\n')

def editMatrix(bm, lbl1, lbl2, num, diff):
    if bm[lbl1][lbl2] > bm[lbl2][lbl1]:
        bm[lbl2][lbl1] = num
        bm[lbl1][lbl2] = bm[lbl2][lbl1] + diff
    else:
        bm[lbl1][lbl2] = num
        bm[lbl2][lbl1] = bm[lbl1][lbl2] + diff

    return bm

def checkOutput(outputFile):
    readFile = file(outputFile)
    checked = None
    count = 0
    for line in reversed(readFile.readlines()):
        if line.startswith('          DESCRIPTION OF VIBRATIONS'):
            checked = 1
        # elif line.startswith(' Error in internal coordinate system'):
        #     checked = 2
        # Don't want to read through the whole things, as success & failure keys are at the bottom
        if checked:
            return checked
        count += 1 

def parse(tsOutput, output1, output2, outputDataFile, reactant, product):
    mol1Parse = cclib.parser.Mopac(output1)
    mol2Parse = cclib.parser.Mopac(output2)
    tsParse   = cclib.parser.Mopac(tsOutput)

    parsed1 = mol1Parse.parse()
    parsed2 = mol2Parse.parse()
    tsParse = tsParse.parse()

    # In kJ/mol
    mol1E = parsed1.scfenergies[-1]
    mol2E = parsed2.scfenergies[-1]
    tsE = tsParse.scfenergies[-1]
    dE = (tsE - mol1E - mol2E) * 1.60218 * 60221.4
    tsVib = tsParse.vibfreqs[0]

    r1String = 'Reactant 1        = ' + reactant.split()[0].toSMILES()
    r2String = 'Reactant 2        = ' + reactant.split()[1].toSMILES()
    p1String = 'Product 1         = ' + product.split()[0].toSMILES()
    p2String = 'Product 2         = ' + product.split()[1].toSMILES()
    tEnergy  = 'Activation Energy = ' + str(dE)
    tVib     = 'TS vib            = ' + str(tsVib)

    with open(outputDataFile, 'w') as parseFile:
        parseFile.write('The energies of the species in kJ/mol are:')
        parseFile.write('\n')
        parseFile.write(r1String)
        parseFile.write('\n')
        parseFile.write(r2String)
        parseFile.write('\n')
        parseFile.write(p1String)
        parseFile.write('\n')
        parseFile.write(p2String)
        parseFile.write('\n')
        parseFile.write(tEnergy)
        parseFile.write('\n')
        parseFile.write(tVib)
        parseFile.write('\n')

def optimizeGeom(outPath, inputPath, qmCalc):
    converge = 0
    if os.path.exists(outPath):
        converge = checkOutput(outPath)
        return converge
    while converge != 1:
        for attempt in range(1, 3):
            if attempt == 1:
                writeInputFile(inputPath, qmCalc.getMolFilePathForCalculation(attempt), qmCalc.geometry)
                run(executablePath, inputPath, outPath)
            elif attempt == 2:
                writeInputFile(inputPath, qmCalc.getMolFilePathForCalculation(attempt + 20), qmCalc.geometry)
                run(executablePath, inputPath, outPath)
            converge = checkOutput(outPath)
            if converge == 1:
                return converge
    return converge

def calcTS(TS, count):
    quantumMechanics = QMCalculator()
    quantumMechanics.settings.software = 'mopac'
    quantumMechanics.settings.fileStore = 'QMfiles'
    quantumMechanics.settings.scratchDirectory = 'scratch'
    quantumMechanics.settings.onlyCyclics = False
    quantumMechanics.settings.maxRadicalNumber = 0

    reactant = fixSortLabel(TS[0])
    product = fixSortLabel(TS[1])

    rRDMol, rBM, rMult, rGeom = generateBoundsMatrix(reactant, quantumMechanics.settings)
    pRDMol, pBM, pMult, pGeom = generateBoundsMatrix(product, quantumMechanics.settings)

    #edit bounds distances to align reacting atoms
    if family.lower() == 'h_abstraction':
        lbl1 = reactant.getLabeledAtom('*1').sortingLabel
        lbl2 = reactant.getLabeledAtom('*2').sortingLabel
        lbl3 = reactant.getLabeledAtom('*3').sortingLabel

        rBM = editMatrix(rBM, lbl1, lbl3, 2.5, 0.1)
        rBM = editMatrix(rBM, lbl2, lbl3, 2.0, 0.1)

        pBM = editMatrix(pBM, lbl1, lbl2, 2.0, 0.1)
        pBM = editMatrix(pBM, lbl1, lbl3, 2.5, 0.1)

    rdkit.DistanceGeometry.DoTriangleSmoothing(rBM)
    rdkit.DistanceGeometry.DoTriangleSmoothing(pBM)

    rsorted_atom_list = reactant.vertices[:]
    qmcalcR = rmgpy.qm.mopac.MopacMolPM7(reactant, quantumMechanics.settings)
    reactant.vertices = rsorted_atom_list

    psorted_atom_list = product.vertices[:]
    qmcalcP = rmgpy.qm.mopac.MopacMolPM7(product, quantumMechanics.settings)
    product.vertices = psorted_atom_list

    qmcalcR.createGeometry(rBM)
    # take the reactant geometry and apply the product bounds matrix
    # this should prevent non-reacting atom overlap
    rRDMol = rdkit.Chem.MolFromMolFile(rGeom.getCrudeMolFilePath(), removeHs=False)

    for atom in reactant.atoms:
        i = atom.sortingLabel
        pRDMol.GetConformer(0).SetAtomPosition(i, rRDMol.GetConformer(0).GetAtomPosition(i))
    atoms = len(pGeom.molecule.atoms)
    distGeomAttempts=1
    if atoms > 3:#this check prevents the number of attempts from being negative
        distGeomAttempts = 5*(atoms-3) #number of conformer attempts is just a linear scaling with molecule size, due to time considerations in practice, it is probably more like 3^(n-3) or something like that

    pGeom.rd_embed(pRDMol, distGeomAttempts, pBM)

    qmcalcP.geometry = pGeom

    geometryR = qmcalcR.geometry
    geometryP = qmcalcP.geometry

    rinputFilePath = qmcalcR.inputFilePath
    routputFilePath = qmcalcR.outputFilePath
    rmolFilePathForCalc = qmcalcR.getMolFilePathForCalculation(attempt)

    pinputFilePath = qmcalcP.inputFilePath
    poutputFilePath = qmcalcR.outputFilePath
    pmolFilePathForCalc = qmcalcP.getMolFilePathForCalculation(attempt)
    inputFilePath = rinputFilePath
    outputFilePath = poutputFilePath

    # TS file paths
    rRefInPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rRef' + inputFileExtension)
    grefInPath1 = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'pGeo' + inputFileExtension)
    pRefInPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'pRef' + inputFileExtension)
    grefInPath2 = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rGeo' + inputFileExtension)
    saddleInPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'saddleCalc' + inputFileExtension)
    tsoptInPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'tsopt' + inputFileExtension)
    tsInPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'transitionState' + inputFileExtension)
    ircInput = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'irc' + inputFileExtension)
    
    grefOutPath1 = grefInPath1.split('.')[0] + outputFileExtension
    grefOutPath2 = grefInPath2.split('.')[0] + outputFileExtension
    saddleOutPath = saddleInPath.split('.')[0] + outputFileExtension
    tsoptOutPath = tsoptInPath.split('.')[0] + outputFileExtension
    tsOutPath = tsInPath.split('.')[0] + outputFileExtension
    ircOutput = ircInput.split('.')[0] + outputFileExtension
    
    if os.path.exists(tsOutPath):
        pass
    else:
        # Write the reactant and product files and have them reference each other
        writeReferenceFile(rRefInPath, rmolFilePathForCalc, geometryR, attempt)
        writeGeoRefInputFile(grefInPath1, pmolFilePathForCalc, rRefInPath, geometryP)
        run(executablePath, grefInPath1, grefOutPath1)
        
        # Write the product input file that references the reactant
        writeReferenceFile(pRefInPath, pmolFilePathForCalc, geometryP, attempt, outputFile=grefOutPath1)
        writeGeoRefInputFile(grefInPath2, rmolFilePathForCalc, pRefInPath, geometryR)
        run(executablePath, grefInPath2, grefOutPath2)
        
        # Write the saddle calculation using outputs from both geo ref calcs
        writeSaddleInputFile(saddleInPath, grefOutPath1, grefOutPath2, geometryR, geometryP)
        run(executablePath, saddleInPath, saddleOutPath)
        
        # Write TS calculation
        writeTSInputFile(tsInPath, saddleOutPath, count)
        run(executablePath, tsInPath, tsOutPath)
    
    # writeTSEval(tsInPath, tsoptOutPath, count)
    # run(executablePath, tsInPath, tsOutPath)
    tsConverge = checkOutput(tsOutPath)
    # Split the reactants and products in order to calculate their energies
    # and generate the geometries
    rct1, rct2 = reactant.split()
    # prd1, prd2 = product.split()
    
    rct1 = fixSortLabel(rct1)
    rct2 = fixSortLabel(rct2)
    # prd1 = fixSortLabel(prd1)
    # prd2 = fixSortLabel(prd2)
    
    r1Qmcalc = rmgpy.qm.mopac.MopacMolPM7(rct1, quantumMechanics.settings)
    r2Qmcalc = rmgpy.qm.mopac.MopacMolPM7(rct2, quantumMechanics.settings)
    # p1Qmcalc = rmgpy.qm.mopac.MopacMolPM7(prd1, quantumMechanics.settings)
    # p2Qmcalc = rmgpy.qm.mopac.MopacMolPM7(prd2, quantumMechanics.settings)
    
    r1Qmcalc.createGeometry()
    r2Qmcalc.createGeometry()
    # p1Qmcalc.createGeometry()
    # p2Qmcalc.createGeometry()
    
    # Reactant and product file paths
    r1InPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct1' + inputFileExtension)
    r1OutPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct1' + outputFileExtension)
    r2InPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct2' + inputFileExtension)
    r2OutPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct2' + outputFileExtension)
    # p1InPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'prd1' + inputFileExtension)
    # p1OutPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'prd1' + outputFileExtension)
    # p2InPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'prd2' + inputFileExtension)
    # p2OutPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'prd2' + outputFileExtension)
    
    # Run the optimizations    
    r1Converge = optimizeGeom(r1OutPath, r1InPath, r1Qmcalc)
    r2Converge = optimizeGeom(r2OutPath, r2InPath, r2Qmcalc)
    # p1Converge = optimizeGeom(p1OutPath, p1InPath, p1Qmcalc)
    # p2Converge = optimizeGeom(p2OutPath, p2InPath, p2Qmcalc)
    
    # Conduct IRC calculation and validate resulting geometries
    if tsConverge == 1:
        writeIRC(ircInput, tsOutPath, count)
        run(executablePath, ircInput, ircOutput)
    
    # Check outputs
    rTest = tsConverge * r1Converge * r2Converge
    # pTest = tsConverge * p1Converge * p2Converge
    
    # Data file
    rOutputDataFile = os.path.join(quantumMechanics.settings.fileStore, 'activationER' + str(count) + outputFileExtension)
    # pOutputDataFile = os.path.join(quantumMechanics.settings.fileStore, 'activationEP' + str(count) + outputFileExtension)
    
    # Parsing, so far just reading energies
    if rTest == 1:
        if os.path.exists(rOutputDataFile):
            pass
        else:
            parse(tsOutPath, r1OutPath, r2OutPath, rOutputDataFile, reactant, product)
    # if pTest == 1:
    #     if os.path.exists(pOutputDataFile):
    #         pass
    #     else:
    #         parse(tsOutPath, p1OutPath, p2OutPath, pOutputDataFile, reactant, product)       

########################################################################################

count = 0
for TS in tsStructures:
    if count > 15:
        calcTS(TS, count)
    count += 1