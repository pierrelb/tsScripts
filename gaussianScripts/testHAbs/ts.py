import os

import logging
import cclib.parser
import openbabel
from subprocess import Popen
from collections import defaultdict

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.molecule import Geometry
# from rmgpy.qm.reaction import QMReaction
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
trusted = open('/home/pierreb/Code/RMG-database/input/kinetics/families/H_Abstraction/training.py')

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
inputFileExtension = '.gjf'
qstOutputFileExtension = '.out'
outputFileExtension = '.log'
executablePath = os.path.join(os.getenv('GAUSS_EXEDIR') , 'g09')
attempt = 1

usePolar = False

keywords = [
            "# pm6 opt=(verytight,gdiis) freq IOP(2/16=3)",
            "# pm6 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",
            "# pm6 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" ,
            "# pm6 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
            "# pm6 opt=(verytight,gdiis,small) freq IOP(2/16=3)",
            "# pm6 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
            "# pm6 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",
            "# pm6 opt=tight freq IOP(2/16=3)",
            "# pm6 opt=tight freq=numerical IOP(2/16=3)",
            "# pm6 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
            "# pm6 opt freq IOP(2/16=3)",
            "# pm6 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
            "# pm6 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
            "# pm6 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
            "# pm6 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
            "# pm6 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
            "# pm6 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
            "# pm6 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
            ]

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
    # submits the input file to Gaussian
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

def writeInputFile(inputFilePath, molFilePathForCalc, geometry, attempt):
    """
    Using the :class:`Geometry` object, write the input file
    for the `attmept`th attempt.
    """

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "gjf")
    mol = openbabel.OBMol()

    obConversion.ReadFile(mol, molFilePathForCalc )

    mol.SetTitle(geometry.uniqueIDlong)
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    top_keys = "# pm6 opt=(nofreeze,calcall,tight,noeigentest) nosymm"
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)
        gaussianFile.write('\n')

def writeQST2InputFile(inputFilePath, rmolFilePathForCalc, pmolFilePathForCalc, geometryR, geometryP, trial):
    chk_file = '%chk=' + inputFilePath.split('.')[0]
    obrConversion = openbabel.OBConversion()
    obrConversion.SetInAndOutFormats("mol", "gjf")
    molR = openbabel.OBMol()
    obrConversion.ReadFile(molR, rmolFilePathForCalc )
    molR.SetTitle(geometryR.uniqueIDlong)
    obrConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)

    obpConversion = openbabel.OBConversion()
    obpConversion.SetInAndOutFormats("mol", "gjf")
    molP = openbabel.OBMol()
    obpConversion.ReadFile(molP, pmolFilePathForCalc )
    molP.SetTitle(geometryP.uniqueIDlong)
    obpConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)

    # all of the first molecule, and remove the first 2 lines (the '\n') from the second
    input_string = obrConversion.WriteString(molR) + obpConversion.WriteString(molP)[2:]
    if trial == 1:
        top_keys = "# pm6 opt=(qst2,nofreeze,calcall,tight,noeigentest) nosymm"
    elif trial == 2:
        top_keys = "# pm6 opt=(qst2,nofreeze,calcall,tight,noeigentest,cartesian) nosymm"
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)
        gaussianFile.write('\n')

def writeTSInputFile(inputFilePath, trial):
    chk_file = '%chk=' + inputFilePath.split('.')[0]
    if trial == 1:
        top_keys = "# pm6 opt=(ts,nofreeze,calcall,tight,noeigentest) geom=allcheck guess=check nosymm"
    elif trial == 2:
        top_keys = "# pm6 opt=(ts,nofreeze,calcall,tight,noeigentest,cartesian) geom=allcheck guess=check nosymm"
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write('\n\n')
        
def writeIRCInput(inputFilePath):
    chk_file = '%chk=' + inputFilePath.split('IRC')[0]
    top_keys = "# pm6 irc=(calcall,report=read) geom=allcheck guess=check nosymm"
    # Normally you need to calculate these
    chg_mult = '1 2'
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write('\n\n')
        gaussianFile.write(chg_mult)
        gaussianFile.write('\n\n')

def checkIRC(ircOutPath):
    readFile = file(ircOutPath)
    lineList = list()
    for line in readFile.readlines():
        if line.startswith(' Point Number:'):
            lineList.append(line)
    
    pth1 = 0
    pth2 = int(lineList[-1].split()[2])
    test = int(lineList[-1].split()[-1])
    for line in lineList:
        if int(line.split()[-1]) == test:
            return pth1, pth2
        else:
            pth1 = int(line.split()[2])
            
def fromXYZ(molGeom, atomnos):
    """
    Takes cartesian coordinates and the associated atoms and generates
    an RMG molecule.
    """
    obMol = cclib.bridge.makeopenbabel(molGeom, atomnos)
    rmgMol = Molecule().fromOBMol(obMol)
    
    return rmgMol

def testGeometries(react, prdct, ircOutPath, notes):
    """
    Compares IRC geometries to input geometries.
    """
    # Search IRC output for total number of steps
    stepNum1, stepNum2 = checkIRC(ircOutPath)
    
    # Compare the reactants and products
    ircParse = cclib.parser.Gaussian(ircOutPath)
    ircParse = ircParse.parse()
    
    atomnos = ircParse.atomnos
    atomcoords = ircParse.atomcoords
    # import ipdb; ipdb.set_trace()
    
    # Convert the IRC geometries into RMG molecules
    rMol = fromXYZ(atomcoords[0], atomnos)
    pMol = fromXYZ(atomcoords[-1], atomnos)
    
    if rMol.toInChI() == react.toInChI():
        if pMol.toInChI() == prdct.toInChI():
            notes = ''
            return 1, notes
        else:
            notes = ' products do not match '
            return 0, notes
    elif rMol.toInChI() == prdct.toInChI():
        if pMol.toInChI() == react.toInChI():
            notes = ''
            return 1, notes
        else:
            notes = ' reactants do not match '
            return 0, notes
    else:
        if pMol.toInChI() == prdct.toInChI():
            notes = 'reactants do not match '
            return 0, notes
        elif pMol.toInChI() == react.toInChI():
            notes = ' products do not match '
            return 0, notes
        else:
            notes = ' reactants and products do not match '
            return 0, notes

def editMatrix(bm, lbl1, lbl2, num, diff):
    if bm[lbl1][lbl2] > bm[lbl2][lbl1]:
        bm[lbl2][lbl1] = num
        bm[lbl1][lbl2] = bm[lbl2][lbl1] + diff
    else:
        bm[lbl1][lbl2] = num
        bm[lbl2][lbl1] = bm[lbl1][lbl2] + diff

    return bm

def checkOutput(outputFile):
    readFile = open(outputFile, 'r')
    readFile = readFile.readlines()
    if readFile[-1].startswith(' Normal termination of'):
        return 1
    else:
        if readFile[-4].startswith(' Error in internal coordinate system'):
            return 0
        else:
            return 2

def parse(tsOutput, output1, output2, outputDataFile, reactant, product):
    mol1Parse = cclib.parser.Gaussian(output1)
    mol2Parse = cclib.parser.Gaussian(output2)
    tsParse = cclib.parser.Gaussian(tsOutput)

    parsed1 = mol1Parse.parse()
    parsed2 = mol2Parse.parse()
    tsParse = tsParse.parse()
    
    # In kJ/mol
    mol1E = parsed1.scfenergies[-1]
    mol2E = parsed2.scfenergies[-1]
    tsE = tsParse.scfenergies[-1]
    dE = (tsE - mol1E - mol2E) * 96.4853365
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
    for attempt in range(1, 3):
        if attempt == 1:
            writeInputFile(inputPath, qmCalc.getMolFilePathForCalculation(attempt), qmCalc.geometry, attempt)
            run(executablePath, inputPath, outPath)
        elif attempt == 2:
            writeInputFile(inputPath, qmCalc.getMolFilePathForCalculation(attempt + 20), qmCalc.geometry, attempt + 20)
            run(executablePath, inputPath, outPath)
        else:
            break
        converge = checkOutput(outPath)
        if converge == 1:
            return converge
    return converge

def calcTS(TS, count):
    quantumMechanics = QMCalculator()
    quantumMechanics.settings.software = 'gaussian'
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
    qmcalcR = rmgpy.qm.gaussian.GaussianMolPM6(reactant, quantumMechanics.settings)
    reactant.vertices = rsorted_atom_list

    psorted_atom_list = product.vertices[:]
    qmcalcP = rmgpy.qm.gaussian.GaussianMolPM6(product, quantumMechanics.settings)
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
        distGeomAttempts = 5*(atoms-3)  #number of conformer attempts is just a linear scaling with molecule size, 
                                        # due to time considerations in practice, it is 
                                        # probably more like 3^(n-3) or something like that

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
    
    tsName = 'transitionState' + str(count)
    
    # TS file paths
    tsFilePath = os.path.join(quantumMechanics.settings.fileStore, tsName + inputFileExtension)
    qstOutPath = os.path.join(quantumMechanics.settings.fileStore, tsName + qstOutputFileExtension)
    tsOutPath = os.path.join(quantumMechanics.settings.fileStore, tsName + outputFileExtension)
    ircInPath = os.path.join(quantumMechanics.settings.fileStore, tsName + 'IRC' + inputFileExtension)
    ircOutPath = os.path.join(quantumMechanics.settings.fileStore, tsName + 'IRC' + outputFileExtension)

    # TS calculations
    tsConverge = 0
    if os.path.exists(tsOutPath):
        tsConverge = checkOutput(tsOutPath)
    trial = 1
    while tsConverge == 0:
        writeQST2InputFile(tsFilePath, rmolFilePathForCalc, pmolFilePathForCalc, geometryR, geometryP, trial)
        run(executablePath, tsFilePath, qstOutPath)
        writeTSInputFile(tsFilePath, trial)
        run(executablePath, tsFilePath, tsOutPath)
        tsConverge = checkOutput(tsOutPath)
        trial += 1
    
    # Do the IRC calculation
    writeIRCInput(ircInPath)
    run(executablePath, ircInPath, ircOutPath)
    
    notes = ''
    
    ircCheck, notes = testGeometries(reactant, product, ircOutPath, notes)
    
    while ircCheck == 1:
        # Split the reactants and products in order to calculate their energies
        # and generate the geometries
        rct1, rct2 = reactant.split()
        
        rct1 = fixSortLabel(rct1)
        rct2 = fixSortLabel(rct2)
        
        r1Qmcalc = rmgpy.qm.gaussian.GaussianMolPM6(rct1, quantumMechanics.settings)
        r2Qmcalc = rmgpy.qm.gaussian.GaussianMolPM6(rct2, quantumMechanics.settings)
        p1Qmcalc = rmgpy.qm.gaussian.GaussianMolPM6(prd1, quantumMechanics.settings)
        p2Qmcalc = rmgpy.qm.gaussian.GaussianMolPM6(prd2, quantumMechanics.settings)
        
        r1Qmcalc.createGeometry()
        r2Qmcalc.createGeometry()
        p1Qmcalc.createGeometry()
        p2Qmcalc.createGeometry()
        
        # Reactant and product file paths
        r1InPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct1' + inputFileExtension)
        r1OutPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct1' + outputFileExtension)
        r2InPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct2' + inputFileExtension)
        r2OutPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct2' + outputFileExtension)
        
        # Run the optimizations    
        r1Converge = optimizeGeom(r1OutPath, r1InPath, r1Qmcalc)
        r2Converge = optimizeGeom(r2OutPath, r2InPath, r2Qmcalc)
        
        # Check outputs
        rTest = tsConverge * r1Converge * r2Converge
        
        # Data file
        rOutputDataFile = os.path.join(quantumMechanics.settings.fileStore, 'activationER' + str(count) + outputFileExtension)
        
        # Parsing, so far just reading energies
        if rTest == 1:
            if os.path.exists(rOutputDataFile):
                pass
            else:
                parse(tsOutPath, r1OutPath, r2OutPath, rOutputDataFile, reactant, product)
        if pTest == 1:
            if os.path.exists(pOutputDataFile):
                pass
            else:
                parse(tsOutPath, p1OutPath, p2OutPath, pOutputDataFile, reactant, product)       

def applyTST():
    """
    For reaction A + B-H ==> A-H + B
    k = Y (k_B * T / h) * Qts/(Qa * Qb-h) * exp(-dE0/(k_B*T))
    Y   : The 'transmission coefficient' in range 0 to 1. Allows for the possibility that not all TS structures 
    go to products.
          It is frequently 1. 
    k_B : Boltzmann Constant
    T   : Temperature
    h   : Planck's Constant
    Qx  : Total partition function of species x
    dE0 : Difference between the lowest energy levels of the reactants and the transition state
    R   : Molar gas constant
    
    k = Y * (k_B * T / h) * exp(-dGts/(R * T))
    """
    pass

########################################################################################

count = 0
for TS in tsStructures:
    calcTS(TS, count)
    count += 1