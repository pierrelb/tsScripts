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
from rmgpy.qm.reaction import QMReaction
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

trusted = open('../../RMG-database/input/kinetics/families/H_Abstraction/depository.py')

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
    top_keys = inputFileKeywords(attempt)
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

# def parse(tsOutput, reactantOutput, productOutput, outputDataFile):
def parse(tsOutput, reactantOutput, outputDataFile):
    rParse = cclib.parser.Gaussian(reactantOutput)
    # pParse = cclib.parser.Gaussian(productOutput)
    tsParse = cclib.parser.Gaussian(tsOutput)
    
    rParse = rParse.parse()
    # pParse = pParse.parse()
    tsParse = tsParse.parse()
    
    # In Hartrees
    reactantE = rParse.scfenergies[-1]/27.2113845
    # productE = pParse.scfenergies[-1]/27.2113845
    tsE = tsParse.scfenergies[-1]/27.2113845
    tsVib = tsParse.vibfreqs[0]
    
    rString = 'Reactant energy = ' + str(reactantE)
    # pString = 'Product energy  = ' + str(productE)
    tEnergy = 'TS energy       = ' + str(tsE)
    tVib    = 'TS vib          = ' + str(tsVib)
    
    with open(outputDataFile, 'w') as parseFile:
        parseFile.write('The energies of the species in Hartree are:')
        parseFile.write('\n')
        parseFile.write(rString)
        parseFile.write('\n')
        # parseFile.write(pString)
        # parseFile.write('\n')
        parseFile.write(tEnergy)
        parseFile.write('\n')
        parseFile.write(tVib)
        parseFile.write('\n')
    
def calculate(TS, count):
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
    qmcalcR = rmgpy.qm.gaussian.GaussianMolPM3(reactant, quantumMechanics.settings)
    reactant.vertices = rsorted_atom_list

    psorted_atom_list = product.vertices[:]
    qmcalcP = rmgpy.qm.gaussian.GaussianMolPM3(product, quantumMechanics.settings)
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
    
    # Reactant and product file paths
    rInPath = os.path.join(quantumMechanics.settings.fileStore, 'reactant' + str(count) + inputFileExtension)
    rOutPath = os.path.join(quantumMechanics.settings.fileStore, 'reactant' + str(count) + outputFileExtension)
    # pInPath = os.path.join(quantumMechanics.settings.fileStore, 'product' + str(count) + inputFileExtension)
    # pOutPath = os.path.join(quantumMechanics.settings.fileStore, 'product' + str(count) + outputFileExtension)
    
    # TS file paths
    tsFilePath = os.path.join(quantumMechanics.settings.fileStore, 'transitionState' + str(count) + inputFileExtension)
    qstOutPath = os.path.join(quantumMechanics.settings.fileStore, 'transitionState' + str(count) + qstOutputFileExtension)
    tsOutPath = os.path.join(quantumMechanics.settings.fileStore, 'transitionState' + str(count) + outputFileExtension)
    
    # Data file
    outputDataFile = os.path.join(quantumMechanics.settings.fileStore, 'kinetics' + str(count) + outputFileExtension)
    
    # Reactant and product calculations
    rAttempt = 1
    rConverge = 0
    if os.path.exists(rOutPath):
        rConverge = checkOutput(rOutPath)
    while rConverge != 1:
        writeInputFile(rInPath, rmolFilePathForCalc, geometryR, rAttempt)
        run(executablePath, rInPath, rOutPath)
        rConverge = checkOutput(rOutPath)
        rAttempt += 1
    
    # pAttempt = 1    
    # pConverge = 0
    # if os.path.exists(pOutPath):
    #     pConverge = checkOutput(pOutPath)
    # while pConverge != 1:
    #     writeInputFile(pInPath, pmolFilePathForCalc, geometryP, pAttempt)
    #     run(executablePath, pInPath, pOutPath)
    #     pConverge = checkOutput(pOutPath)
    #     pAttempt += 1
    
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
            
    # Parsing, so far just reading energies
    if tsConverge == 1:
        if os.path.exists(outputDataFile):
            pass
        else:
            # parse(tsOutPath, rOutPath, pOutPath, outputDataFile)
            parse(tsOutPath, rOutPath, outputDataFile)        

########################################################################################

count = 0
for TS in tsStructures:
    calculate(TS, count)
    count += 1

