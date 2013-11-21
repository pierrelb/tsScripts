import os

import logging
import external.cclib.parser
import openbabel
import time
from subprocess import Popen
from collections import defaultdict, Counter

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.molecule import Geometry, RDKitFailedError
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe, saveEntry
from rmgpy.reaction import Reaction
from rmgpy.data.base import Entry
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
qstOutputFileExtension = '.out'
outputFileExtension = '.log'
executablePath = os.path.join(os.getenv('GAUSS_EXEDIR') , 'g09')
attempt = 1

usePolar = False

keywords = [
            "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis) freq IOP(2/16=3)",
            "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis) freq IOP(2/16=3) IOP(4/21=2)",
            "# b3lyp/6-31+g(d,p) opt=(verytight,calcfc,maxcyc=200) freq IOP(2/16=3) nosymm" ,
            "# b3lyp/6-31+g(d,p) opt=(verytight,calcfc,maxcyc=200) freq=numerical IOP(2/16=3) nosymm",
            "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis,small) freq IOP(2/16=3)",
            "# b3lyp/6-31+g(d,p) opt=(verytight,nolinear,calcfc,small) freq IOP(2/16=3)",
            "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis,maxcyc=200) freq=numerical IOP(2/16=3)",
            "# b3lyp/6-31+g(d,p) opt=tight freq IOP(2/16=3)",
            "# b3lyp/6-31+g(d,p) opt=tight freq=numerical IOP(2/16=3)",
            "# b3lyp/6-31+g(d,p) opt=(tight,nolinear,calcfc,small,maxcyc=200) freq IOP(2/16=3)",
            "# b3lyp/6-31+g(d,p) opt freq IOP(2/16=3)",
            "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis) freq=numerical IOP(2/16=3) IOP(4/21=200)",
            "# b3lyp/6-31+g(d,p) opt=(calcfc,verytight,newton,notrustupdate,small,maxcyc=100,maxstep=100) freq=(numerical,step=10) IOP(2/16=3) nosymm",
            "# b3lyp/6-31+g(d,p) opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
            "# b3lyp/6-31+g(d,p) opt=(tight,gdiis,small,maxcyc=200,maxstep=100) freq=numerical IOP(2/16=3) nosymm",
            "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis,calcall,small,maxcyc=200) IOP(2/16=3) IOP(4/21=2) nosymm",
            "# b3lyp/6-31+g(d,p) opt=(verytight,gdiis,calcall,small) IOP(2/16=3) nosymm",
            "# b3lyp/6-31+g(d,p) opt=(calcall,small,maxcyc=100) IOP(2/16=3)",
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
    top_keys = "# b3lyp/6-31+g(d,p) opt=(nofreeze,calcall,tight,noeigentest) int=ultrafine nosymm"
    numProc = "%NProcShared=4"
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(numProc)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)
        gaussianFile.write('\n')

def writeQST2InputFile(inputFilePath, rmolFilePathForCalc, pmolFilePathForCalc, geometryR, geometryP):
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
    top_keys = "# b3lyp/6-31+g(d,p) opt=(qst2,nofreeze,calcall,noeigentest) nosymm"
    numProc = "%NProcShared=4"
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(numProc)
        gaussianFile.write('\n')
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)
        gaussianFile.write('\n')

def writeTSInputFile(inputFilePath, trial):
    chk_file = '%chk=' + inputFilePath.split('.')[0]
    if trial == 1:
        top_keys = "# b3lyp/6-31+g(d,p) opt=(ts,nofreeze,calcall,tight,noeigentest) int=ultrafine geom=allcheck guess=check nosymm"
    elif trial == 2:
        top_keys = "# b3lyp/6-31+g(d,p) opt=(ts,nofreeze,calcall,tight,noeigentest,cartesian) int=ultrafine geom=allcheck guess=check nosymm"
    numProc = "%NProcShared=4"    
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(numProc)
        gaussianFile.write('\n')
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write('\n\n')

def writeIRCInput(inputFilePath):
    chk_file = '%chk=' + inputFilePath.split('IRC')[0]
    top_keys = "# b3lyp/6-31+g(d,p) irc=(calcall,report=read) geom=allcheck guess=check nosymm"
    # Normally you need to calculate these
    chg_mult = '1 2'
    numProc = "%NProcShared=4"    
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

def editMatrix(bm, lbl1, lbl2, lbl3, num, diff):
    if bm[lbl1][lbl3] > bm[lbl3][lbl1]:
        bm[lbl3][lbl1] = num
        bm[lbl1][lbl3] = num + diff
    else:
        bm[lbl3][lbl1] = num + diff
        bm[lbl1][lbl3] = num

    if bm[lbl2][lbl3] == 1000. or bm[lbl3][lbl2] == 1000.:
        checkSide = True
    else:
        checkSide = False

    if checkSide:
        if bm[lbl2][lbl3] > bm[lbl3][lbl2]:
            if bm[lbl1][lbl2] > bm[lbl2][lbl1]:
                bm[lbl3][lbl2] = num - bm[lbl2][lbl1] - diff / 2.0
                bm[lbl2][lbl3] = num - bm[lbl1][lbl2] + diff / 2.0
            else:
                bm[lbl3][lbl2] = num - bm[lbl1][lbl2] - diff / 2.0
                bm[lbl2][lbl3] = num - bm[lbl2][lbl1] + diff / 2.0
        else:
            if bm[lbl1][lbl2] > bm[lbl2][lbl1]:
                bm[lbl2][lbl3] = num - bm[lbl2][lbl1] - diff / 2.0
                bm[lbl3][lbl2] = num - bm[lbl1][lbl2] + diff / 2.0
            else:
                bm[lbl2][lbl3] = num - bm[lbl1][lbl2] - diff / 2.0
                bm[lbl3][lbl2] = num - bm[lbl2][lbl1] + diff / 2.0
    else:
        if bm[lbl1][lbl2] > bm[lbl2][lbl1]:
            if bm[lbl2][lbl3] > bm[lbl3][lbl2]:
                bm[lbl2][lbl1] = num - bm[lbl3][lbl2] - diff / 2.0
                bm[lbl1][lbl2] = num - bm[lbl2][lbl3] + diff / 2.0
            else:
                bm[lbl2][lbl1] = num - bm[lbl2][lbl3] - diff / 2.0
                bm[lbl1][lbl2] = num - bm[lbl3][lbl2] + diff / 2.0
        else:
            if bm[lbl2][lbl3] > bm[lbl3][lbl2]:
                bm[lbl1][lbl2] = num - bm[lbl3][lbl2] - diff / 2.0
                bm[lbl2][lbl1] = num - bm[lbl2][lbl3] + diff / 2.0
            else:
                bm[lbl1][lbl2] = num - bm[lbl2][lbl3] - diff / 2.0
                bm[lbl2][lbl1] = num - bm[lbl3][lbl2] + diff / 2.0

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

def writeRxnOutputFile(self):outputDataFile,reactant, product, vibFreq, activeAts, atomDist, notes):

    item = Reaction(reactants=reactant.split(), products=product.split())
    distances = {'d12':atomDist[0], 'd23':atomDist[1], 'd13':atomDist[2]}
    date = time.asctime()
    user = "Pierre Bhoorasingh <bhoorasingh.p@husky.neu.edu>"
    description = "Found via group estimation strategy using automatic transition state generator"
    entry = Entry(
        index = 1,
        item = self.reaction,
        data = DistanceData(distances=distances, method='B3LYP/6-31+G(d,p)'),
        shortDesc = "B3LYP/6-31+G(d,p) calculation via double-ended TS generator.",
        history = [(date, user, 'action', description)]
    )
    star3 = product.getLabeledAtom('*1').sortingLabel
    star1 = product.getLabeledAtom('*3').sortingLabel
    product.atoms[star1].label = '*1'
    product.atoms[star3].label = '*3'

    with open(outputDataFile, 'w') as parseFile:
        saveEntry(parseFile, entry)

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

    # decrease the vdwRadii
    for i in range(len(rBM)):
        for k in range(len(rBM)):
            if rBM[i][k] == 1000.0:
                rBM[k][i] -= 1.2
            if pBM[i][k] == 1000.0:
                pBM[k][i] -= 1.2

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
    tsOutputDataFile = os.path.join(quantumMechanics.settings.fileStore, 'data' + str(count) + outputFileExtension)

    # If a data file exists, it was already successfully run
    if not os.path.exists(tsOutputDataFile):
        # TS calculations
        tsConverge = 0
        if os.path.exists(tsOutPath):
            tsConverge = checkOutput(tsOutPath)
        else:
            writeQST2InputFile(tsFilePath, rmolFilePathForCalc, pmolFilePathForCalc, geometryR, geometryP)
            run(executablePath, tsFilePath, qstOutPath)
            trial = 1
            while tsConverge == 0:
                writeTSInputFile(tsFilePath, trial)
                run(executablePath, tsFilePath, tsOutPath)
                tsConverge = checkOutput(tsOutPath)
                trial += 1

        if not os.path.exists(ircOutPath):
            # Do the IRC calculation
            writeIRCInput(ircInPath)
            run(executablePath, ircInPath, ircOutPath)

        notes = ''

        ircCheck, notes = testGeometries(reactant, product, ircOutPath, notes)

        if not os.path.exists(tsOutputDataFile) and ircCheck == 1:
            # Split the reactants and products in order to calculate their energies
            # and generate the geometries
            rct1, rct2 = reactant.split()

            # rct1 = fixSortLabel(rct1)
            # rct2 = fixSortLabel(rct2)
            # 
            # r1Qmcalc = rmgpy.qm.gaussian.GaussianMolPM6(rct1, quantumMechanics.settings)
            # r2Qmcalc = rmgpy.qm.gaussian.GaussianMolPM6(rct2, quantumMechanics.settings)
            # 
            # r1Qmcalc.createGeometry()
            # r2Qmcalc.createGeometry()
            # 
            # # Reactant and product file paths
            # r1InPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct1' + inputFileExtension)
            # r1OutPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct1' + outputFileExtension)
            # r2InPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct2' + inputFileExtension)
            # r2OutPath = os.path.join(quantumMechanics.settings.fileStore, str(count) + 'rct2' + outputFileExtension)
            # 
            # # Run the optimizations
            # r1Converge = 0 
            # r2Converge = 0
            # r1Converge = optimizeGeom(r1OutPath, r1InPath, r1Qmcalc)
            # r2Converge = optimizeGeom(r2OutPath, r2InPath, r2Qmcalc)

            # # TS calculations
            # tsConverge = 0
            # if os.path.exists(tsOutPath):
            #     tsConverge = checkOutput(tsOutPath)

            # # Check outputs
            # rTest = tsConverge * r1Converge * r2Converge
            # 
            # Data file

            # # Parsing, so far just reading energies
            # if tsConverge == 1:
            #     if os.path.exists(rOutputDataFile):
            #         pass
            #     else:
            #         parse(tsOutPath, rOutputDataFile, reactant, product, notes)
            vibFreq, activeAts, atomDist = parse(tsOutPath, tsOutputDataFile, labels)
            writeRxnOutputFile(tsOutputDataFile, reactant, product, vibFreq, activeAts, atomDist, notes)
            # parse(tsOutPath, rOutputDataFile, reactant, product, notes)

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


for count,TS in enumerate(tsStructures):
    calcTS(TS, count)
