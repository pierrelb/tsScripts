import os
import logging
import external.cclib.parser
import openbabel
import time
from subprocess import Popen
from collections import defaultdict, Counter
from copy import deepcopy

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.molecule import Geometry
from rmgpy.qm.reaction import QMReaction
from rmgpy.data.base import Entry
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe, saveEntry
from rmgpy.data.kinetics.transitionstates import TransitionStates, DistanceData

from rmgpy.qm.gaussian import GaussianTSB3LYP

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

trusted = open(os.path.join(os.getenv('HOME'),'Code/RMG-database/input/kinetics/families/H_Abstraction/NIST.py'))

tsPath = os.path.join(os.getenv('HOME'), 'Code/RMG-database/input/kinetics/families/H_Abstraction')
local_context = None
global_context = None
transitionStates = TransitionStates()
transitionStates.load(tsPath, local_context, global_context)

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

prevReactions = list()
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
    rxnInChI = [rStruct[0].toInChI(), rStruct[1].toInChI(), pStruct[0].toInChI(), pStruct[1].toInChI()]
    doubleChk = 0
    for pair in prevReactions:
        if Counter(pair) == Counter(rxnInChI):
            doubleChk = 1
    if doubleChk == 0:
        prevReactions.append(rxnInChI)
        tsStructures.append(tsStruct)

########################################################################################
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

def calculate(TS, count):
    quantumMechanics = QMCalculator()
    quantumMechanics.settings.software = 'mopac'
    quantumMechanics.settings.fileStore = 'QMfiles'
    quantumMechanics.settings.scratchDirectory = 'scratch'
    quantumMechanics.settings.onlyCyclics = False
    quantumMechanics.settings.maxRadicalNumber = 0

    reactant = fixSortLabel(TS[0])
    product = fixSortLabel(TS[1])

    reaction = Reaction(label='H_Abstraction', reactants=reactant.split(), products=product.split(), reversible=True)

    qmReaction = MopacPM7(reaction, quantumMechanics.settings)

    qmReaction.generateGeometry()



for count, TS in enumerate(tsStructures):
    calculate(TS, count)
    count += 1

########################################################################################        
inputFileExtension = '.mop'
inputFileExtension2 = '.gjf'
outputFileExtension = '.out'
executablePath = os.path.join(os.getenv('MOPAC_DIR', default="/opt/mopac") , 'MOPAC2012.exe')
executablePath2 = os.path.join(os.getenv('GAUSS_EXEDIR') , 'g09')
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
        mopacFile.write(input_string)
        mopacFile.write('\n')
        mopacFile.write(bottom_keys)
        mopacFile.write('\n')

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
        reload(openbabel)
        mol = external.cclib.bridge.makeopenbabel(atomcoords, atomnos)
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
    reload(openbabel)
    molR = external.cclib.bridge.makeopenbabel(atomcoords, atomnos)
    
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
    reload(openbabel)
    molP = external.cclib.bridge.makeopenbabel(atomcoords, atomnos)
    
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

def writeDFTTSInputFile(inputFilePath, outputFile, count):
    
    parseOutput = external.cclib.parser.Mopac(outputFile)
    parseOutput = parseOutput.parse()
    reload(openbabel)
    mol = external.cclib.bridge.makeopenbabel(parseOutput.atomcoords[0], parseOutput.atomnos)
    mol.SetTotalSpinMultiplicity(2)
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "gjf")
    mol.SetTitle('transitionState' + str(count))
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    numProc = "%NProcShared=4"
    top_keys = "# m062x/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest,cartesian) int=ultrafine nosymm"
    
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(numProc)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)

def writeTSInputFile(inputFilePath, saddleOutput, count):
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "mop")
    parseOutput = external.cclib.parser.Mopac(saddleOutput)
    parseOutput = parseOutput.parse()
    reload(openbabel)
    mol = external.cclib.bridge.makeopenbabel(parseOutput.atomcoords[0], parseOutput.atomnos)
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
    parseOutput = external.cclib.parser.Mopac(tsOutPath)
    parseOutput = parseOutput.parse()
    reload(openbabel)
    mol = external.cclib.bridge.makeopenbabel(parseOutput.atomcoords[0], parseOutput.atomnos)
    mol.SetTitle('transitionState' + str(count))
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    
    top_keys = 'irc=1*'
    with open(inputFilePath, 'w') as mopacFile:
        mopacFile.write(top_keys)
        mopacFile.write(input_string)
        mopacFile.write('\n')

def getIRCGeom(lines):
    geom1 = []
    for line in lines:
        if not line.startswith('  reversed'):
            geom1.append(line)
        else:
            break
    geom1.pop()
    
    geom2 = []
    for line in reversed(lines):
        if not line.startswith(' DRC'):
            geom2.insert(0, line)
        else:
            break
    
    return geom1, geom2

def convertMol(geomLines):
    atomcoords = []
    atomnos = []
    for line in geomLines:
        atType, x, y, z = line.split()
        if atType == 'H':
            atNum = 1
        elif atType == 'C':
            atNum = 6
        elif atType == 'O':
            atNum = 8
        coords = [float(x),float(y),float(z)]
        atomnos.append(atNum)
        atomcoords.append(coords)
    atomnos = numpy.array(atomnos, dtype=int)
    atomcoords = numpy.array(atomcoords)
    reload(openbabel)
    mol = external.cclib.bridge.makeopenbabel(atomcoords, atomnos)
    
    rmgMol = Molecule().fromOBMol(mol)
    
    return rmgMol

def parseIRC(ircOutput, reactant, product):
    # You want to read the `.xyz` file.
    ircXYZ = ircOutput.split('.')[0] + '.xyz'
    xyzFile = file(ircXYZ)
    readLines = xyzFile.readlines()
    
    # Remove the first 2 lines from the `.xyz`. Makes it easier to get the
    # geometries we want.
    readLines.pop(0)
    readLines.pop(0)
    
    geom1, geom2 = getIRCGeom(readLines)
    
    ircMol1 = convertMol(geom1)
    ircMol2 = convertMol(geom2)
    
    product.resetConnectivityValues()
    reactant.resetConnectivityValues()
    ircMol1.resetConnectivityValues()
    ircMol2.resetConnectivityValues()
    
    # Compare IRC geometries iwth
    if reactant.isIsomorphic(ircMol1) and product.isIsomorphic(ircMol2):
        return 1
    elif reactant.isIsomorphic(ircMol2) and product.isIsomorphic(ircMol1):
        return 1
    else:
        return 0

def editMatrix(bm, lbl1, lbl2, lbl3, num, diff):
    if bm[lbl1][lbl3] > bm[lbl3][lbl1]:
        bm[lbl3][lbl1] = num
        bm[lbl1][lbl3] = num + diff
    else:
        bm[lbl3][lbl1] = num + diff
        bm[lbl1][lbl3] = num

    # if bm[lbl2][lbl3] == 1000. or bm[lbl3][lbl2] == 1000.:
    #     checkSide = True
    # else:
    #     checkSide = False
    # 
    # if checkSide:
    #     if bm[lbl2][lbl3] > bm[lbl3][lbl2]:
    #         if bm[lbl1][lbl2] > bm[lbl2][lbl1]:
    #             bm[lbl3][lbl2] = num - bm[lbl2][lbl1] - diff / 2.0
    #             bm[lbl2][lbl3] = num - bm[lbl1][lbl2] + diff / 2.0
    #         else:
    #             bm[lbl3][lbl2] = num - bm[lbl1][lbl2] - diff / 2.0
    #             bm[lbl2][lbl3] = num - bm[lbl2][lbl1] + diff / 2.0
    #     else:
    #         if bm[lbl1][lbl2] > bm[lbl2][lbl1]:
    #             bm[lbl2][lbl3] = num - bm[lbl2][lbl1] - diff / 2.0
    #             bm[lbl3][lbl2] = num - bm[lbl1][lbl2] + diff / 2.0
    #         else:
    #             bm[lbl2][lbl3] = num - bm[lbl1][lbl2] - diff / 2.0
    #             bm[lbl3][lbl2] = num - bm[lbl2][lbl1] + diff / 2.0
    # else:
    #     if bm[lbl1][lbl2] > bm[lbl2][lbl1]:
    #         if bm[lbl2][lbl3] > bm[lbl3][lbl2]:
    #             bm[lbl2][lbl1] = num - bm[lbl3][lbl2] - diff / 2.0
    #             bm[lbl1][lbl2] = num - bm[lbl2][lbl3] + diff / 2.0
    #         else:
    #             bm[lbl2][lbl1] = num - bm[lbl2][lbl3] - diff / 2.0
    #             bm[lbl1][lbl2] = num - bm[lbl3][lbl2] + diff / 2.0
    #     else:
    #         if bm[lbl2][lbl3] > bm[lbl3][lbl2]:
    #             bm[lbl1][lbl2] = num - bm[lbl3][lbl2] - diff / 2.0
    #             bm[lbl2][lbl1] = num - bm[lbl2][lbl3] + diff / 2.0
    #         else:
    #             bm[lbl1][lbl2] = num - bm[lbl2][lbl3] - diff / 2.0
    #             bm[lbl2][lbl1] = num - bm[lbl3][lbl2] + diff / 2.0

    return bm

def checkOutput(outputFile):
    readFile = file(outputFile)
    checked = None
    count = 0
    for line in reversed(readFile.readlines()):
        if line.startswith('          DESCRIPTION OF VIBRATIONS'):
            checked = 1
        elif line.startswith(' ** GRADIENT IS TOO LARGE TO ALLOW FORCE'):
            checked = 2
        # Don't want to read through the whole things, as success & failure keys are at the bottom
        if checked:
            return checked
        count += 1 

def getAtomType(atomnum):
    if atomnum == 1:
        atType = 'H'
    elif atomnum == 6:
        atType = 'C'
    elif atomnum == 8:
        atType = 'O'
    else:
        atType = 'Could not determine atomtype'
    
    return atType

def parse(tsOutput, output1, output2, outputDataFile, labels):
    mol1Parse = external.cclib.parser.Mopac(output1)
    mol2Parse = external.cclib.parser.Mopac(output2)
    tsParse   = external.cclib.parser.Mopac(tsOutput)

    parsed1 = mol1Parse.parse()
    parsed2 = mol2Parse.parse()
    tsParse = tsParse.parse()

    # In J/mol
    mol1E = parsed1.scfenergies[-1]
    mol2E = parsed2.scfenergies[-1]
    tsE = tsParse.scfenergies[-1]
    deltaE = (tsE - mol1E - mol2E) * 1.60218 * 60221.4
    vibFreq = tsParse.vibfreqs[0]
    
    atom1 = openbabel.OBAtom()
    atom2 = openbabel.OBAtom()
    atom3 = openbabel.OBAtom()
    
    atom1.SetAtomicNum(int(tsParse.atomnos[labels[0]]))
    atom2.SetAtomicNum(int(tsParse.atomnos[labels[1]]))
    atom3.SetAtomicNum(int(tsParse.atomnos[labels[2]]))
    
    atom1coords = tsParse.atomcoords[0][labels[0]].tolist()
    atom2coords = tsParse.atomcoords[0][labels[1]].tolist()
    atom3coords = tsParse.atomcoords[0][labels[2]].tolist()
    
    atom1.SetVector(*atom1coords)
    atom2.SetVector(*atom2coords)
    atom3.SetVector(*atom3coords)
    
    at1 = getAtomType(atom1.GetAtomicNum())
    at2 = getAtomType(atom2.GetAtomicNum())
    at3 = getAtomType(atom3.GetAtomicNum())
    
    activeAts = [at1, at2, at3]
    atomDist = [str(atom1.GetDistance(atom2)), str(atom2.GetDistance(atom3)), str(atom1.GetDistance(atom3))]
    
    return deltaE, vibFreq, activeAts, atomDist
    
def writeRxnOutputFile(outputDataFile,reactant, product, deltaE, vibFreq, activeAts, atomDist, notes):
    r1String = 'Reactant 1        = ' + reactant.split()[0].toSMILES()
    r2String = 'Reactant 2        = ' + reactant.split()[1].toSMILES()
    p1String = 'Product 1         = ' + product.split()[0].toSMILES()
    p2String = 'Product 2         = ' + product.split()[1].toSMILES()
    tEnergy  = 'Activation Energy = ' + str(deltaE)
    tVib     = 'TS vib            = ' + str(vibFreq)
    define1  = '*1                = ' + activeAts[0]
    define2  = '*2                = ' + activeAts[1]
    define3  = '*3                = ' + activeAts[2]
    dist12   = '*1 to *2          = ' + atomDist[0]
    dist23   = '*2 to *3          = ' + atomDist[1]
    dist13   = '*1 to *3          = ' + atomDist[2]
    notes    = 'Notes             = ' + notes
    
    with open(outputDataFile, 'w') as parseFile:
        parseFile.write('The energies of the species in J/mol are:')
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
        parseFile.write(define1)
        parseFile.write('\n')
        parseFile.write(define2)
        parseFile.write('\n')
        parseFile.write(define3)
        parseFile.write('\n')
        parseFile.write(dist12)
        parseFile.write('\n')
        parseFile.write(dist23)
        parseFile.write('\n')
        parseFile.write(dist13)
        parseFile.write('\n')
        parseFile.write(notes)
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
    
    if count < 10:
        fileNum = '00' + str(count)
    elif count < 100:
        fileNum = '0' + str(count)
    else:
        fileNum = str(count)
    
    # Keep track of any failures
    notes = ''
    
    reactant = fixSortLabel(TS[0])
    product = fixSortLabel(TS[1])

    rRDMol, rBM, rMult, rGeom = generateBoundsMatrix(reactant, quantumMechanics.settings)
    pRDMol, pBM, pMult, pGeom = generateBoundsMatrix(product, quantumMechanics.settings)
    
    # # decrease the vdwRadii
    # for i in range(len(rBM)):
    #     for k in range(len(rBM)):
    #         if rBM[i][k] == 1000.0:
    #             rBM[k][i] -= 1.2
    #         if pBM[i][k] == 1000.0:
    #             pBM[k][i] -= 1.2
            
    #edit bounds distances to align reacting atoms
    if family.lower() == 'h_abstraction':
        at1 = reactant.getLabeledAtom('*1')
        at2 = reactant.getLabeledAtom('*2')
        at3 = reactant.getLabeledAtom('*3')
        
        lbl1 = at1.sortingLabel
        lbl2 = at2.sortingLabel
        lbl3 = at3.sortingLabel
        
        labels = [lbl1, lbl2, lbl3]

        if (at1.symbol == 'H' and at3.symbol == 'C') or (at1.symbol == 'C' and at3.symbol == 'H'):
            rBM = editMatrix(rBM, lbl1, lbl2, lbl3, 2.2, 0.1)
            pBM = editMatrix(pBM, lbl3, lbl2, lbl1, 2.2, 0.1)
        elif (at1.symbol == 'H' and at3.symbol == 'O') or (at1.symbol == 'O' and at3.symbol == 'H'):
            rBM = editMatrix(rBM, lbl1, lbl2, lbl3, 2.1, 0.1)
            pBM = editMatrix(pBM, lbl3, lbl2, lbl1, 2.1, 0.1)
        elif at1.symbol == 'O' and at3.symbol == 'O':
            rBM = editMatrix(rBM, lbl1, lbl2, lbl3, 2.2, 0.1)
            pBM = editMatrix(pBM, lbl3, lbl2, lbl1, 2.2, 0.1)
        else:
            rBM = editMatrix(rBM, lbl1, lbl2, lbl3, 2.5, 0.1)
            pBM = editMatrix(pBM, lbl3, lbl2, lbl1, 2.5, 0.1)
    
    sect = len(reactant.split()[1].atoms)
    
    for i in range(sect,len(rBM)):
        for j in range(0,sect):
            for k in range(len(rBM)):
                if k==i or k==j: continue
                Uik = rBM[i,k] if k>i else rBM[k,i]
                Ukj = rBM[j,k] if k>j else rBM[k,j]
                
                maxLij = Uik + Ukj - 0.15
                if rBM[i,j] >  maxLij:
                    # print "CHANGING {0} to {1}".format(rBM[i,j], maxLij)
                    rBM[i,j] = maxLij
    
    for i in range(sect,len(pBM)):
        for j in range(0,sect):
            for k in range(len(pBM)):
                if k==i or k==j: continue
                Uik = pBM[i,k] if k>i else pBM[k,i]
                Ukj = pBM[j,k] if k>j else pBM[k,j]
                
                maxLij = Uik + Ukj - 0.15
                if pBM[i,j] >  maxLij:
                    # print "CHANGING {0} to {1}".format(pBM[i,j], maxLij)
                    pBM[i,j] = maxLij

    setRBM = rdkit.DistanceGeometry.DoTriangleSmoothing(rBM)
    setPBM = rdkit.DistanceGeometry.DoTriangleSmoothing(pBM)

    if setRBM and setPBM:
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
        distGeomAttempts=15
        if atoms > 3:#this check prevents the number of attempts from being negative
            distGeomAttempts = 15*(atoms-3) # number of conformer attempts is just a linear scaling with molecule size, 
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
        
        # TS file paths
        rRefInPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'rRef' + inputFileExtension)
        grefInPath1 = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'pGeo' + inputFileExtension)
        pRefInPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'pRef' + inputFileExtension)
        grefInPath2 = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'rGeo' + inputFileExtension)
        saddleInPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'saddleCalc' + inputFileExtension)
        tsoptInPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'tsopt' + inputFileExtension)
        tsInPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'transitionState' + inputFileExtension)
        ircInput = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'irc' + inputFileExtension)
        
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
        
        tsConverge = checkOutput(tsOutPath)
        rightGeom = 0
        # Conduct IRC calculation and validate resulting geometries
        if tsConverge == 1:
            writeIRC(ircInput, tsOutPath, count)
            run(executablePath, ircInput, ircOutput)
            rightGeom = parseIRC(ircOutput, reactant, product)
        else:
            notes = notes + 'Transition state not converged: '
        
        r1Converge = 0
        r2Converge = 0
        if rightGeom == 1:
            # Split the reactants and products in order to calculate their energies
            # and generate the geometries
            
            rct1, rct2 = reactant.split()
            
            rct1 = fixSortLabel(rct1)
            rct2 = fixSortLabel(rct2)
            
            r1Qmcalc = rmgpy.qm.mopac.MopacMolPM7(rct1, quantumMechanics.settings)
            r2Qmcalc = rmgpy.qm.mopac.MopacMolPM7(rct2, quantumMechanics.settings)
            
            r1Qmcalc.createGeometry()
            r2Qmcalc.createGeometry()
            
            # Reactant and product file paths
            r1InPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'rct1' + inputFileExtension)
            r1OutPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'rct1' + outputFileExtension)
            r2InPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'rct2' + inputFileExtension)
            r2OutPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'rct2' + outputFileExtension)
            
            # Run the optimizations
            if len(rct1.atoms)==1:
                run(executablePath, r1InPath, r1OutPath)
                r1Converge = 1
            else:
                r1Converge = optimizeGeom(r1OutPath, r1InPath, r1Qmcalc)
                if r1Converge != 1:
                    notes = notes + 'Failure at reactant 1: '
            if len(rct2.atoms)==1:
                run(executablePath, r2InPath, r2OutPath)
                r2Converge = 1
            else:
                r2Converge = optimizeGeom(r2OutPath, r2InPath, r2Qmcalc)
                if r2Converge != 1:
                    notes = notes + 'Failure at reactant 2: '
        else:
            notes = notes + 'Failure at IRC: '
        
        # Check outputs
        rTest = tsConverge * r1Converge * r2Converge
        
        # Data file
        rOutputDataFile = os.path.join(quantumMechanics.settings.fileStore, 'data' + fileNum + outputFileExtension)
        
        # Parsing, so far just reading energies
        if rTest == 1:
            if os.path.exists(rOutputDataFile):
                pass
            else:
                deltaE, vibFreq, activeAts, atomDist = parse(tsOutPath, r1OutPath, r2OutPath, rOutputDataFile, labels)
        else:
            deltaE = 'Failed run'
            vibFreq = 'Failed run'
            activeAts = ['Failed run', 'Failed run', 'Failed run']
            atomDist = ['Failed run', 'Failed run', 'Failed run']
        
        writeRxnOutputFile(rOutputDataFile, reactant, product, deltaE, vibFreq, activeAts, atomDist, notes)
    else:
        print count

def calcDFTTS(TS, count):
    quantumMechanics = QMCalculator()
    quantumMechanics.settings.software = 'gaussian'
    quantumMechanics.settings.fileStore = 'QMfiles'
    quantumMechanics.settings.scratchDirectory = 'scratch'
    quantumMechanics.settings.onlyCyclics = False
    quantumMechanics.settings.maxRadicalNumber = 0
    
    if count < 10:
        fileNum = '00' + str(count)
    elif count < 100:
        fileNum = '0' + str(count)
    else:
        fileNum = str(count)
    
    reactant = fixSortLabel(TS[0])
    product = fixSortLabel(TS[1])
    
    # TS file paths
    tsInPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'transitionState' + inputFileExtension)
    tsDFTInPath = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'transitionStateDFT' + inputFileExtension2)
    ircInput = os.path.join(quantumMechanics.settings.fileStore, fileNum + 'irc' + inputFileExtension)
    
    tsOutPath = tsInPath.split('.')[0] + outputFileExtension
    tsDFTOutPath = tsDFTInPath.split('.')[0] + outputFileExtension
    ircOutput = ircInput.split('.')[0] + outputFileExtension
    
    # Split the reactants and products in order to calculate their energies
    # and generate the geometries
    writeDFTTSInputFile(tsDFTInPath, tsOutPath, count)
    run(executablePath2, tsDFTInPath, tsDFTOutPath)
        
########################################################################################
for count, TS in enumerate(tsStructures):
    letsRun = False
    for testMol in TS[0].split():
        if testMol.toSMILES() == '[H][H]' or testMol.toSMILES() == '[H]':
            letsRun = True
    
    if letsRun:
        calcTS(TS, count)