import os
import logging
import openbabel
from subprocess import Popen

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.molecule import QMMolecule
from rmgpy.qm.reaction import QMReaction

"""
degeneracy = 1
C(=C[CH2])C=C to C1CC=C[CH]1
"""
family = 'Intra_R_Add_Endocyclic'

reactant = """
1  *5 C 0 {2,S} {3,D} {7,S}
2  *2 C 0 {1,S} {4,D} {6,S}
3  *4 C 0 {1,D} {5,S} {8,S}
4  *3 C 0 {2,D} {9,S} {10,S}
5  *1 C 1 {3,S} {11,S} {12,S}
6     H 0 {2,S}
7     H 0 {1,S}
8     H 0 {3,S}
9     H 0 {4,S}
10    H 0 {4,S}
11    H 0 {5,S}
12    H 0 {5,S}
"""

product = """
1  *3 C 0 {2,S} {3,S} {6,S} {7,S}
2  *1 C 0 {1,S} {4,S} {8,S} {9,S}
3  *2 C 1 {1,S} {5,S} {10,S}
4  *4 C 0 {2,S} {5,D} {11,S}
5  *5 C 0 {3,S} {4,D} {12,S}
6     H 0 {1,S}
7     H 0 {1,S}
8     H 0 {2,S}
9     H 0 {2,S}
10    H 0 {3,S}
11    H 0 {4,S}
12    H 0 {5,S}
"""

actions = [
            ['CHANGE_BOND', '*2', '-1', '*3'],
            ['FORM_BOND', '*1', 'S', '*3'],
            ['GAIN_RADICAL', '*2', '1'],
            ['LOSE_RADICAL', '*1', '1']
            ]
#########################################################################

inputFileExtension = '.gjf'
outputFileExtension = '.out'
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

def run():
    # submits the input file to Gaussian
    process = Popen([executablePath, inputFilePath, outputFilePath])
    process.communicate()# necessary to wait for executable termination!
    
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

def writeModRedundantFile():
    chk_file = '%chk=' + chkFilePath
    top_keys = "# pm3 opt=(modredundant) geom=(allcheck) guess=check nosymm"
    bottom_keys = 'B 1 2 += 0.2 F'
    with open(inputFilePath, 'w') as gaussianFile:
        gaussianFile.write(chk_file)
        gaussianFile.write('\n')
        gaussianFile.write(top_keys)
        gaussianFile.write('\n\n')
        gaussianFile.write(bottom_keys)
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
    
#########################################################################

reactant = Molecule().fromAdjacencyList(reactant)
product = Molecule().fromAdjacencyList(product)

if reactant.isCyclic():
    tsBase = reactant
elif product.isCyclic():
    tsBase = product
else:
    logging.info('not able to generate the transition state')

quantumMechanics = QMCalculator()
quantumMechanics.settings.software = 'gaussian'
quantumMechanics.settings.fileStore = 'QMfiles'
quantumMechanics.settings.scratchDirectory = None
quantumMechanics.settings.onlyCyclics = True
quantumMechanics.settings.maxRadicalNumber = 0

qmcalc = rmgpy.qm.gaussian.GaussianMolPM3(tsBase, quantumMechanics.settings)
qmcalc.createGeometry()
geometry = qmcalc.geometry
inputFilePath = qmcalc.inputFilePath
outputFilePath = qmcalc.outputFilePath
chkFilePath = os.path.join(qmcalc.settings.fileStore, geometry.uniqueID)
molFilePathForCalc = qmcalc.getMolFilePathForCalculation(attempt)
writeInputFile()
run()
# i = 1
writeModRedundantFile()
# while i in range(1, 11):
#     
#     convertOutputToInput()
#     run()
#     i += 1