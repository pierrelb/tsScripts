import os
import logging
import openbabel
from subprocess import Popen
from collections import defaultdict, Counter

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.molecule import Geometry
from rmgpy.qm.reaction import QMReaction
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe
from rmgpy.data.kinetics.transitionstates import TransitionStates, DistanceData

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

inputFileExtension = '.gjf'
outputFileExtension = '.log'
executablePath = os.path.join(os.getenv('GAUSS_EXEDIR') , 'g09')
attempt = 1

usePolar = False

keywords = [
            "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3)",
            "# pm3 opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",
            "# pm3 opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" ,
            "# pm3 opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
            "# pm3 opt=(verytight,gdiis,small) freq IOP(2/16=3)",
            "# pm3 opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
            "# pm3 opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",
            "# pm3 opt=tight freq IOP(2/16=3)",
            "# pm3 opt=tight freq=numerical IOP(2/16=3)",
            "# pm3 opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
            "# pm3 opt freq IOP(2/16=3)",
            "# pm3 opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
            "# pm3 opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
            "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
            "# pm3 opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
            "# pm3 opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
            "# pm3 opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
            "# pm3 opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
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
    #if not os.path.exists(geometry.getCrudeMolFilePath()):
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

    return rdKitMol, boundsMatrix, multiplicity

def writeTSInputFile(inputFilePath, molFilePathForCalc, geometry, family):
    chk_file = '%chk=' + inputFilePath.split('.')[0]
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "gjf")
    mol = openbabel.OBMol()

    obConversion.ReadFile(mol, molFilePathForCalc )

    mol.SetTitle(geometry.uniqueIDlong)
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    numProc = "%NProcShared=4"
    top_keys = "# m062x/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest)  int=ultrafine nosymm"
    title = ' ' + geometry.uniqueIDlong + '' + family
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(numProc)
        gaussianFile.write('\n')
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)

def writeTSCartInput(inputFilePath, count):
    chk_file = '%chk=' + inputFilePath.split('.')[0]
    top_keys = "# m062x/6-31+g(d,p) opt=(ts,calcall,tight,noeigentest,cartesian)  int=ultrafine nosymm geom=allcheck guess=check nosymm"
    numProc = "%NProcShared=4"

    # Normally you need to calculate these
    chg_mult = '1 2'
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(numProc)
        gaussianFile.write('\n')
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write('\n\n')
        gaussianFile.write(chg_mult)
        gaussianFile.write('\n\n')

def writeSubFile(filename):
    fout = open(filename + '.sh', "w")
    fout.write('#!/bin/sh\n')
    fout.write('#BSUB -n 4\n')
    fout.write('#BSUB -o ' + filename + '.log\n')
    fout.write('#BSUB -e error' + filename.split('e')[-1] + '\n')
    fout.write('#BSUB -J ' + filename.split('e')[-1] + '\n\n')
    fout.write('export GAUSS_EXEDIR=/shared/g09\n')
    fout.write('export PATH=$GAUSS_EXEDIR:$PATH\n\n')
    fout.write('g09 < ' + filename + '.gjf' + '\n\n')
    fout.close()

def writeIRCInput(inputFilePath, count):
    chk_file = '%chk=' + inputFilePath.split('IRC')[0] + str(count)
    top_keys = "# m062x/6-31+g(d,p) irc=(calcall,report=read) geom=allcheck guess=check nosymm"
    numProc = "%NProcShared=4"

    # Normally you need to calculate these
    chg_mult = '1 2'
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(numProc)
        gaussianFile.write('\n')
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
    obMol = external.cclib.bridge.makeopenbabel(molGeom, atomnos)
    rmgMol = Molecule().fromOBMol(obMol)

    return rmgMol

def testGeometries(reactant,product,ircOutput,notes):
    """
    Compares IRC geometries to input geometries.
    """
    # Search IRC output for steps on each side of the path
    readFile = file(ircOutput)
    pth1 = list()
    steps = list()
    for line in readFile.readlines():
        if line.startswith(' Point Number:'):
            if int(line.split()[2]) > 0:
                if int(line.split()[-1]) == 1:
                    ptNum = int(line.split()[2])
                    pth1.append(ptNum)
                else:
                    pass
        elif line.startswith('  # OF STEPS ='):
            numStp = int(line.split()[-1])
            steps.append(numStp)

    # This indexes the coordinate to be used from the parsing
    if steps == []:
        notes = ' IRC failed '
        return 0, notes
    else:
        pth1End = sum(steps[:pth1[-1]])		
        # Compare the reactants and products
        ircParse = external.cclib.parser.Gaussian(ircOutput)
        ircParse = ircParse.parse()

        atomnos = ircParse.atomnos
        atomcoords = ircParse.atomcoords

        # Convert the IRC geometries into RMG molecules
        # We don't know which is reactant or product, so take the two at the end of the
        # paths and compare to the reactants and products

        mol1 = fromXYZ(atomcoords[pth1End], atomnos)
        mol2 = fromXYZ(atomcoords[-1], atomnos)

        # Had trouble with isIsomorphic, but resetting the connectivity seems to fix it (WHY??).
        reactant.resetConnectivityValues()
        product.resetConnectivityValues()
        mol1.resetConnectivityValues()
        mol2.resetConnectivityValues()

        if reactant.isIsomorphic(mol1) and product.isIsomorphic(mol2):
                notes = 'Verified TS'
                return 1, notes
        elif reactant.isIsomorphic(mol2) and product.isIsomorphic(mol1):
                notes = 'Verified TS'
                return 1, notes
        else:
            notes = 'Saddle found, but wrong one'
            return 0, notes

def generateKineticData():
    pass

def editMatrix(bm, lbl1, lbl2, num, diff):
    if lbl1 > lbl2:
        bm[lbl2][lbl1] = num + diff
        bm[lbl1][lbl2] = num
    else:
        bm[lbl2][lbl1] = num
        bm[lbl1][lbl2] = num + diff

    return bm

def checkOutput(outputFile):
    readFile = open(outputFile, 'r')
    readFile = readFile.readlines()
    i = 0
    for line in reversed(readFile):
        if line.startswith(' Normal termination of'):
            return 1
        elif line.startswith(' Error in internal coordinate system'):
            return 2
        elif i == 20:
            return 0
        i += 1

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

def parse(tsOutput, outputDataFile, labels):

    tsParse = external.cclib.parser.Gaussian(tsOutput)
    tsParse = tsParse.parse()

    # In J/mol
    vibFreq = tsParse.vibfreqs[0]

    atom1 = openbabel.OBAtom()
    atom2 = openbabel.OBAtom()
    atom3 = openbabel.OBAtom()

    atom1.SetAtomicNum(int(tsParse.atomnos[labels[0]]))
    atom2.SetAtomicNum(int(tsParse.atomnos[labels[1]]))
    atom3.SetAtomicNum(int(tsParse.atomnos[labels[2]]))

    atom1coords = tsParse.atomcoords[-1][labels[0]].tolist()
    atom2coords = tsParse.atomcoords[-1][labels[1]].tolist()
    atom3coords = tsParse.atomcoords[-1][labels[2]].tolist()

    atom1.SetVector(*atom1coords)
    atom2.SetVector(*atom2coords)
    atom3.SetVector(*atom3coords)

    at1 = getAtomType(atom1.GetAtomicNum())
    at2 = getAtomType(atom2.GetAtomicNum())
    at3 = getAtomType(atom3.GetAtomicNum())

    activeAts = [at1, at2, at3]
    atomDist = [str(atom1.GetDistance(atom2)), str(atom2.GetDistance(atom3)), str(atom1.GetDistance(atom3))]

    return vibFreq, activeAts, atomDist

def writeRxnOutputFile(outputDataFile,reactant, product, vibFreq, activeAts, atomDist, notes):
    # Parse the data into training set format
    
    item = Reaction(reactants=reactant.split(), products=product.split())
    distances = {'d12':atomDist[0], 'd23':atomDist[1], 'd13':atomDist[2]}
    date = time.asctime()
    user = "Pierre Bhoorasingh <bhoorasingh.p@husky.neu.edu>"
    description = "Found via direct estimation using automatic transition state generator"
    entry = Entry(
        index = 1,
        item = item,
        data = DistanceData(distances=distances, method='M06-2X/6-31+G(d,p)'),
        shortDesc = "M06-2X/6-31+G(d,p) calculation via group additive automatic TS estimator.",
        history = [(date, user, 'action', description)]
    )

    with open(outputDataFile, 'w') as parseFile:
        saveEntry(parseFile, entry)

def calculate(TS, count):
    quantumMechanics = QMCalculator()
    quantumMechanics.settings.software = 'gaussian'
    quantumMechanics.settings.fileStore = 'QMfiles'
    quantumMechanics.settings.scratchDirectory = 'scratch'
    quantumMechanics.settings.onlyCyclics = False
    quantumMechanics.settings.maxRadicalNumber = 0
    
    reactant = fixSortLabel(TS[0])
    product = fixSortLabel(TS[1])
    
    rRDMol, tsBM, tsMult = generateBoundsMatrix(reactant, quantumMechanics.settings)
    
    # edit bounds distances to align reacting atoms
    if family.lower() == 'h_abstraction':
        sect = len(reactant.split()[1].atoms)
    
        tsBM[sect:,:sect] = 1.8
    
        lbl1 = reactant.getLabeledAtom('*1').sortingLabel
        lbl2 = reactant.getLabeledAtom('*2').sortingLabel
        lbl3 = reactant.getLabeledAtom('*3').sortingLabel
        
        labels = [lbl1, lbl2, lbl3]
        atomMatch = ((lbl1,),(lbl2,),(lbl3,))
        
        reaction = Reaction(reactants=reactant.split(), products=product.split())
        distanceData = transitionStates.estimateDistances(reaction)
        
        tsBM = editMatrix(tsBM, lbl1, lbl2, distanceData.distances['d12'], 0.1)
        tsBM = editMatrix(tsBM, lbl2, lbl3, distanceData.distances['d23'], 0.001)
        tsBM = editMatrix(tsBM, lbl1, lbl3, distanceData.distances['d13'], 0.001)
    
    setBM = rdkit.DistanceGeometry.DoTriangleSmoothing(tsBM)
    
    if setBM:
        tssorted_atom_list = reactant.vertices[:]
        qmcalcTS = rmgpy.qm.gaussian.GaussianMolPM3(reactant, quantumMechanics.settings)
        reactant.vertices = tssorted_atom_list
    
        qmcalcTS.createGeometry(tsBM,atomMatch)
    
        geometryTS = qmcalcTS.geometry
        tsmolFilePathForCalc = qmcalcTS.getMolFilePathForCalculation(attempt)
    
        # Create filenames
        tsFilePath = os.path.join(quantumMechanics.settings.fileStore, 'transitionState' + str(count) + '.gjf')
        tsOutPath = os.path.join(quantumMechanics.settings.fileStore, 'transitionState' + str(count) + '.log')
        ircInPath = os.path.join(quantumMechanics.settings.fileStore, 'transitionStateIRC' + str(count) + '.gjf')
        ircOutPath = os.path.join(quantumMechanics.settings.fileStore, 'transitionStateIRC' + str(count) + '.log')
        tsOutputDataFile = os.path.join(quantumMechanics.settings.fileStore, 'data' + str(count) + outputFileExtension)
    
        # QM saddle search
        # Write and run the TS optimization
        tsConverge = 0
        if os.path.exists(tsOutPath):
            tsConverge = checkOutput(tsOutPath)
        else:
            writeTSInputFile(tsFilePath, tsmolFilePathForCalc, geometryTS, family)
            run(executablePath, tsFilePath, tsOutPath)
            tsConverge = checkOutput(tsOutPath)
    
        # Validation
        # Run in internal coodinates, if 2*pi angle is achieved, switch to cartesian coordinates
        if tsConverge == 2:
            # Error in internal coodinate system, continue calculation in cartesian
            writeTSCartInput(tsFilePath, count)
            run(executablePath, tsFilePath, tsOutPath)
            tsConverge = checkOutput(tsOutPath)
        
        # If saddle found, write and run the IRC calculation to check the TS geometry
        if tsConverge == 1:
            writeIRCInput(ircInPath, count)
            run(executablePath, ircInPath, ircOutPath)
            if os.path.exists(ircOutPath):
                ircCheck, notes = testGeometries(reactant, product, ircOutPath, notes)
                if ircCheck == 1:
                    vibFreq, activeAts, atomDist = parse(tsOutPath, tsOutputDataFile, labels)
                    writeRxnOutputFile(tsOutputDataFile, reactant, product, vibFreq, activeAts, atomDist, notes)

########################################################################################
    
count = 0
for TS in tsStructures:
    for molecule in TS[0].split():
        calculate(TS, count)
    count += 1
