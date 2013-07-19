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

def atoms(mol):
    atoms = {}
    for atom in mol:
        args = atom.split()
        index = int(args.pop(0))
        if '*' in args[0]:
            label = args.pop(0)
        else:
            label = '  '
        type = args.pop(0)
        rad = args.pop(0)
        bonds = {}
        while args:
            bond = args.pop(0)[1:-1].split(',')
            bonds[int(bond[0])] = bond[1]
        atoms[index] = {'label': label, 'type': type,
                        'rad': rad, 'bonds': bonds}
    return atoms

def adjlist(atoms):
    str = ''
    for key in atoms:
        atom = atoms[key]
        str += '\n{0:<{1}}{2}'.format(key,
                                      len('{0}'.format(max(atoms.keys()))) + 1,
                                      atom['label'])
        str += ' {0} {1}'.format(atom['type'], atom['rad'])
        for key0 in sorted(atom['bonds'].keys()):
            str += ' {' + '{0},{1}'.format(key0, atom['bonds'][key0]) + '}'
    return str.strip() + '\n'

def bondForm(fullString, otherIdx, bond):
    fullString = fullString + ' {' + str(otherIdx) + ',' + bond + '}'
    
    return fullString
    
def bondBreak(fullString, otherIdx):
    splits = fullString.split('{')
    i = 0
    for lineSplit in splits:
        if lineSplit.split(',')[0] == str(otherIdx):
            splits.pop(i)
        i += 1
    fullString = splits[0]
    for k in range(1, len(splits)):
        fullString = fullString + '{' + splits[k]
    return fullString

def radChange(fullString, action, decrease = False):
    
    radChg = fullString[8]
    if decrease:
        radNum = int(radChg) - int(action)
    else:
        radNum = int(radChg) + int(action)
    radNum = str(radNum)
    fullString = fullString.replace(' ' + radChg + ' ', ' ' + radNum + ' ')
    
    return fullString

def matchAtoms(reactant):
    newadjlist = reactant.toAdjacencyList().strip().splitlines()
    radjlist = reactant.toAdjacencyList().strip().splitlines()
    rdict = {}
    for line in radjlist:
        if line.find('*') > -1:
            rdict[line.split()[1]] = int(line.split()[0])
    
    for action in actions:
        if action[0].lower() == 'break_bond':
            idx1 = rdict[action[1]]
            idx2 = rdict[action[3]]
            
            edit1 = newadjlist.pop(idx1 - 1)
            edit1 = bondBreak(edit1, idx2)
            newadjlist.insert(idx1 - 1, edit1)
            
            edit2 = newadjlist.pop(idx2 - 1)
            edit2 = bondBreak(edit2, idx1)
            newadjlist.insert(idx2 - 1, edit2)
        elif action[0].lower() == 'form_bond':
            idx1 = rdict[action[1]]
            idx2 = rdict[action[3]]
            
            edit1 = newadjlist.pop(idx1 - 1)
            edit1 = bondForm(edit1, idx2, action[2])
            newadjlist.insert(idx1 - 1, edit1)
            
            edit2 = newadjlist.pop(idx2 - 1)
            edit2 = bondForm(edit2, idx1, action[2])
            newadjlist.insert(idx2 - 1, edit2)
        elif action[0].lower() == 'gain_radical':
            idx = rdict[action[1]]
            
            edit = newadjlist.pop(idx - 1)
            edit = radChange(edit, action[2])
            newadjlist.insert(idx - 1, edit)
        elif action[0].lower() == 'lose_radical':
            idx = rdict[action[1]]
            
            edit = newadjlist.pop(idx - 1)
            edit = radChange(edit, action[2], decrease = True)
            newadjlist.insert(idx - 1, edit)
    return newadjlist
        
def writeInputFile():
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
    chk_file = '%chk=' + chkFilePath
    top_keys = inputFileKeywords(attempt)
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
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
    top_keys = "# pm6 opt=(qst2,calcfc,nofreeze) nosymm"
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)
        gaussianFile.write('\n')

def writeTSInputFile(inputFilePath, molFilePathForCalc, geometry, family):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "gjf")
    mol = openbabel.OBMol()
    
    obConversion.ReadFile(mol, molFilePathForCalc )
    
    mol.SetTitle(geometry.uniqueIDlong)
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    top_keys = inputFileKeywords(attempt)
    
    top_keys = "# b3lyp/6-31+g(d,p) opt=(ts,calcall,noeigentest) geom=allcheck guess=check nosymm"
    title = ' ' + geometry.uniqueIDlong + '' + family
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write('\n\n')
        gaussianFile.write(title)
        gaussianFile.write('\n\n')

def writeModRedundantFile():
    chk_file = '%chk=' + chkFilePath
    top_keys = "# pm3 opt=(modredundant) geom=(allcheck) guess=check nosymm"
    bottom_keys = 'B 1 2 += 0.1 F'
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write('\n\n')
        gaussianFile.write(bottom_keys)
        gaussianFile.write('\n')
        
def writeModRedundantFile1():
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol", "gjf")
    mol = openbabel.OBMol()
    
    obConversion.ReadFile(mol, molFilePathForCalc )
    
    mol.SetTitle(geometry.uniqueIDlong)
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    chk_file = '%chk=' + chkFilePath
    top_keys = "# pm3 opt=(modredundant) nosymm"
    bottom_keys1 = 'B 4 9 += 0.1 F'
    bottom_keys2 = 'B 5 9 += -0.1 F'
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)
        gaussianFile.write(bottom_keys1)
        gaussianFile.write('\n')
        gaussianFile.write(bottom_keys2)
        gaussianFile.write('\n')
        
def writeModRedundantFile2():
    chk_file = '%chk=' + chkFilePath
    top_keys = "# pm3 opt=(modredundant) geom=(allcheck) guess=check nosymm"
    bottom_keys1 = 'B 4 9 += 0.1 F'
    bottom_keys2 = 'B 5 9 += -0.1 F'
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write('\n\n')
        gaussianFile.write(bottom_keys1)
        gaussianFile.write('\n')
        gaussianFile.write(bottom_keys2)
        gaussianFile.write('\n')
            
def convertOutputToInput():
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("g09", "gjf")
    mol = openbabel.OBMol()
    
    obConversion.ReadFile(mol, outputFilePath)
    
    mol.SetTitle(geometry.uniqueIDlong)
    obConversion.SetOptions('k', openbabel.OBConversion.OUTOPTIONS)
    input_string = obConversion.WriteString(mol)
    chk_file = '%chk=' + chkFilePath
    top_keys = "# pm3 opt=(modredundant)"
    bottom_keys = '1 2 F'
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write(input_string)
        gaussianFile.write('\n')
        gaussianFile.write(bottom_keys)
        gaussianFile.write('\n')

def generateKineticData():
    pass

def editMatrix(bm, lbl1, lbl2, num, diff):
    if lbl1 > lbl2:
        bm[lbl2][lbl1] = bm[lbl2][lbl1] + diff
        bm[lbl1][lbl2] = num
    else:
        bm[lbl2][lbl1] = num
        bm[lbl1][lbl2] = bm[lbl2][lbl1] + diff
    
    return bm

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
