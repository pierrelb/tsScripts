import os
import sys

import rmgpy
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.data.kinetics import KineticsFamily, ReactionRecipe, KineticsDatabase
from rmgpy.data.rmg import RMGDatabase
from rmgpy.qm.main import QMCalculator
from rmgpy.qm.gaussian import GaussianTSB3LYP

if len(sys.argv)>1:
	i = int(sys.argv[-1])
elif os.getenv('LSB_JOBINDEX'):
	i = int(os.getenv('LSB_JOBINDEX'))
else:
	raise Exception("Specify a TS number!")
	
rxnFamiles = ['intra_H_migration', 'R_Addition_MultipleBond', 'H_Abstraction'] 

# Create a python dictionary of the species name and the adjacency list
# from the species dictionary
moleculeDict = {}
f = open('dictionary.txt', 'r')
adjlist = ''; label = ''
for line in f:
	if len(line.strip()) == 0:
		if len(adjlist.strip()) > 0:
			molecule = Molecule()
			molecule.fromAdjacencyList(adjlist, saturateH=True)
			moleculeDict[label] = molecule
		adjlist = ''; label = ''
	else:
		if len(adjlist.strip()) == 0:
			label = line.strip()
		adjlist += line


file_object = open('mechanism.txt', 'r')
mechLines = file_object.readlines()

rxnList = []
gotit = []
for rxnFamily in rxnFamiles:
	for k, line in enumerate(mechLines):
		if line.startswith('! {0} estimate:'.format(rxnFamily)) or line.startswith('! {0} exact:'.format(rxnFamily)):
			reaction = mechLines[k+1].split()[0]
			if reaction not in gotit:
				gotit.append(reaction)
				rxnList.append((rxnFamily, mechLines[k+1]))
		elif '{0} estimate:'.format(rxnFamily) in line or '{0} exact:'.format(rxnFamily) in line:
			reaction = mechLines[k+1].split()[0]
			if reaction not in gotit:
				gotit.append(reaction)
				rxnList.append((rxnFamily, mechLines[k+1]))

print 'Loading RMG Database ...'
rmgDatabase = RMGDatabase()
rmgDatabase.load(os.path.abspath(os.path.join(os.getenv('RMGpy'), '..', 'RMG-database', 'input')), kineticsFamilies='default')
print 'Finished loading RMG Database ...'

reactionTuple = rxnList[i-1]
rxnFamily, reactionLine = reactionTuple
rxnFormula, A, n, Ea = reactionLine.split()
reactants, products = rxnFormula.split('=')
if rxnFamily in ['H_Abstraction', 'Disproportionation']:
	reactant1, reactant2 = [moleculeDict[j] for j in reactants.split('+')]
	product1, product2 = [moleculeDict[j] for j in products.split('+')]
	rSpecies1 = Species(molecule=[reactant1])
	rSpecies2 = Species(molecule=[reactant2])
	pSpecies1 = Species(molecule=[product1])
	pSpecies2 = Species(molecule=[product2])
	rSpecies1.generateResonanceIsomers()
	rSpecies2.generateResonanceIsomers()
	pSpecies1.generateResonanceIsomers()
	pSpecies2.generateResonanceIsomers()
	testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies1, pSpecies2], reversible=True)
	reactionList = []
	for moleculeA in rSpecies1.molecule:
		for moleculeB in rSpecies2.molecule:
			tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[rxnFamily])
			for rxn0 in tempList:
				reactionList.append(rxn0)
elif rxnFamily in ['intra_H_migration']:
	reactant = moleculeDict[reactants]
	product = moleculeDict[products]
	rSpecies = Species(molecule=[reactant])
	pSpecies = Species(molecule=[product])
	rSpecies.generateResonanceIsomers()
	pSpecies.generateResonanceIsomers()
	testReaction = Reaction(reactants=[rSpecies], products=[pSpecies], reversible=True)
	reactionList = []
	for moleculeA in rSpecies.molecule:
		tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA], [], only_families=[rxnFamily])
		for rxn0 in tempList:
			reactionList.append(rxn0)
elif rxnFamily in ['R_Addition_MultipleBond']:
	if '(+M)' in reactants:
		reactants = reactants.split('(+M)')[0]
		products = products.split('(+M)')[0]
	if len(reactants.split('+'))==2:
		reactant1, reactant2 = [moleculeDict[j] for j in reactants.split('+')]
		product = moleculeDict[products]
	else:
		reactant1, reactant2 = [moleculeDict[j] for j in products.split('+')]
		product = moleculeDict[reactants]
	rSpecies1 = Species(molecule=[reactant1])
	rSpecies2 = Species(molecule=[reactant2])
	pSpecies = Species(molecule=[product])
	rSpecies1.generateResonanceIsomers()
	rSpecies2.generateResonanceIsomers()
	pSpecies.generateResonanceIsomers()
	testReaction = Reaction(reactants=[rSpecies1, rSpecies2], products=[pSpecies], reversible=False)
	reactionList = []
	for moleculeA in rSpecies1.molecule:
		for moleculeB in rSpecies2.molecule:
			tempList = rmgDatabase.kinetics.generateReactionsFromFamilies([moleculeA, moleculeB], [], only_families=[rxnFamily])
			for rxn0 in tempList:
				reactionList.append(rxn0)
gotOne=False
for reaction in reactionList:
	# Check if any of the RMG proposed reactions matches the reaction in the mechanism
	if reaction.isIsomorphic(testReaction):
		# Now add the labeled atoms to the Molecule, and check all labels were added
		atLblsR = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
		atLblsP = dict([(lbl[0], False) for lbl in reaction.labeledAtoms])
		
		for reactant in reaction.reactants:
			reactant = reactant.molecule[0]
			reactant.clearLabeledAtoms()
			for atom in reactant.atoms:
				for atomLabel in reaction.labeledAtoms:
					if atom==atomLabel[1]:
						atom.label = atomLabel[0]
						atLblsR[atomLabel[0]] = True
		for product in reaction.products:
			product = product.molecule[0]
			product.clearLabeledAtoms()
			for atom in product.atoms:
				for atomLabel in reaction.labeledAtoms:
					if atom==atomLabel[1]:
						atom.label = atomLabel[0]
						atLblsP[atomLabel[0]] = True
		if all( atLblsR.values() ) and all( atLblsP.values() ):
			gotOne=True
			break

def calculate(reaction):
	reaction = qmCalc.getKineticData(reaction)

	for files in os.listdir('./'):
		if files.startswith('core'):
			os.remove(files)

if not gotOne:
	print "No reactions found for reaction {4}: {0} + {1} = {2} + {3}".format(rSpecies1.molecule[0].toSMILES(), rSpecies2.molecule[0].toSMILES(), pSpecies1.molecule[0].toSMILES(), pSpecies2.molecule[0].toSMILES(), i)
else:
	qmCalc = QMCalculator(
									software='gaussian',
									method='b3lyp',
									fileStore='QMfiles',
									scratchDirectory='QMscratch',
									)
	calculate(reaction)
