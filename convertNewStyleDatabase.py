"""
Loads the TS_training.py file and removes duplicate reactions, then the .data files generated
by the automated TS generator are added (again removing any duplicate reactions) and printed to
a new TS_training.py in the folder where this is run.

Reverse reactions are also duplicates, so run `duplicate.py` to generate those.
"""

import os
import logging
import numpy

from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.transitionstates import DistanceData, Entry

def loadTraining(index,
			  	 reactant1=None,
			  	 reactant2=None,
			  	 reactant3=None,
			  	 product1=None,
			  	 product2=None,
			  	 product3=None,
			  	 distances=None,
			  	 degeneracy=1,
			  	 label='',
			  	 reference=None,
			  	 referenceType='',
			  	 shortDesc='',
			  	 longDesc='',
			  	 rank=None,
			  	 ):
	
	reactants = [Molecule().fromAdjacencyList(reactant1, saturateH=True)]
	if reactant2 is not None: reactants.append(Molecule().fromAdjacencyList(reactant2, saturateH=True))
	if reactant3 is not None: reactants.append(Molecule().fromAdjacencyList(reactant3, saturateH=True))
	
	products = [Molecule().fromAdjacencyList(product1, saturateH=True)]
	if product2 is not None: products.append(Molecule().fromAdjacencyList(product2, saturateH=True))
	if product3 is not None: products.append(Molecule().fromAdjacencyList(product3, saturateH=True))
	
	reaction = Reaction(reactants=reactants, products=products)
	alreadyHave = False
	# for testReaction in reactionList:
	# 	if testReaction.isIsomorphic(reaction):
	# 		alreadyHave = True
	
	if not alreadyHave:
		entry = Entry(
			index = index,
			label = label,
			item = reaction,
			data = distances,
			reference = reference,
			referenceType = referenceType,
			shortDesc = shortDesc,
			longDesc = longDesc.strip(),
			rank = rank,
		)
		reactionList.append(reaction)
		entries['{0:d}:{1}'.format(index,label)] = entry
		return entry
	
def loadEntry(index,
			  reactant1=None,
			  reactant2=None,
			  reactant3=None,
			  product1=None,
			  product2=None,
			  product3=None,
			  distances=None,
			  degeneracy=1,
			  label='',
			  reference=None,
			  referenceType='',
			  shortDesc='',
			  longDesc='',
			  rank=None,
			  ):
	
	reactants = [Molecule().fromAdjacencyList(reactant1, saturateH=True)]
	if reactant2 is not None: reactants.append(Molecule().fromAdjacencyList(reactant2, saturateH=True))
	if reactant3 is not None: reactants.append(Molecule().fromAdjacencyList(reactant3, saturateH=True))
	
	products = [Molecule().fromAdjacencyList(product1, saturateH=True)]
	if product2 is not None: products.append(Molecule().fromAdjacencyList(product2, saturateH=True))
	if product3 is not None: products.append(Molecule().fromAdjacencyList(product3, saturateH=True))
	
	reaction = Reaction(reactants=reactants, products=products)
	alreadyHave = False
	for testReaction in reactionList:
		if testReaction.isIsomorphic(reaction):
			alreadyHave = True
	
	if not alreadyHave:
		index = len(entries) + 1
		
		entry = Entry(
			index = index,
			label = label,
			item = reaction,
			data = distances,
			reference = reference,
			referenceType = referenceType,
			shortDesc = shortDesc,
			longDesc = longDesc.strip(),
			rank = rank,
		)
		entries['{0:d}:{1}'.format(index,label)] = entry
		return entry
		
def createLabels(entryList):
	"""
	Create the labels and the species dictionary.
	"""
	def compareMols(molecule, refMol):
		"""
		Check if the molecule is already in the dictionary. This is not simply and isomorphism
		check, as the atom labels are also relevant.
		"""
		# Start with an isomorphism check
		if not refMol.isIsomorphic(molecule):
			return False
		
		if refMol.getLabeledAtoms().keys()!=molecule.getLabeledAtoms().keys():
			return False
		
		for label in refMol.getLabeledAtoms():
			refAtm = refMol.getLabeledAtom(label)
			otherAtm = molecule.getLabeledAtom(label)
			if not (refAtm.connectivity1==otherAtm.connectivity1 and refAtm.connectivity2==otherAtm.connectivity2 and refAtm.connectivity3==otherAtm.connectivity3):
				return False
		
		return True
						
	specDict = {}
	for i in range(len(entryList)):
		entry = entryList[str(i+1)+':']
		label = ''
		for reactant in entry.item.reactants:
			if label != '':
				label = label + ' + '
			
			specNm = reactant.getFormula()
			
			while specNm in specDict.keys():
				dictAdjlist = specDict[specNm]
				dictMol = Molecule().fromAdjacencyList(dictAdjlist)
				sameMol = compareMols(reactant, dictMol)
				if sameMol:
					break
				try:
					specNm.index('-')
					oldString = list(specNm)
					idx = int(oldString[-1])
					oldString[-1] = str(idx + 1)
					specNm = "".join(oldString)
				except ValueError:
					specNm = specNm + "-2"
				specDict[specNm] = reactant.toAdjacencyList()
			else:
				specDict[specNm] = reactant.toAdjacencyList()
			label = label + specNm
		
		label = label + ' <=> '
		
		for product in entry.item.products:
			if not label.endswith('<=> '):
				label = label + ' + '
			
			specNm = product.getFormula()
			
			while specNm in specDict.keys():
				dictAdjlist = specDict[specNm]
				dictMol = Molecule().fromAdjacencyList(dictAdjlist)
				sameMol = compareMols(product, dictMol)
				if sameMol:
					break
				try:
					specNm.index('-')
					oldString = list(specNm)
					idx = int(oldString[-1])
					oldString[-1] = str(idx + 1)
					specNm = "".join(oldString)
				except ValueError:
					specNm = specNm + "-2"
				specDict[specNm] = product.toAdjacencyList()
			else:
				specDict[specNm] = product.toAdjacencyList()
			label = label + specNm
			
		entry.label = label
	
	return specDict
		
def saveEntries(entryList):
	"""
	Save the entries.
	"""
	with open('reactions.py', 'w') as resultFile:
		resultFile.write('#!/usr/bin/env python\n')
		resultFile.write('# encoding: utf-8\n\n')
		resultFile.write('name = "H_Abstraction/TS_training"\n')
		resultFile.write('shortDesc = u"Distances used to train group additivity values for TS geometries"\n')
		resultFile.write('longDesc = u"""\nPut interatomic distances for reactions to use as a training set for fitting\ngroup additivity values in this file.\n"""\n')
		resultFile.write('recommended = True\n\n')
		for i in range(len(entryList)):
			entry = entryList[str(i+1)+':']
			resultFile.write('entry(\n')
			resultFile.write('    index = {0},\n'.format(entry.index))
			resultFile.write('    label = {0},\n'.format(entry.label))
			resultFile.write('    distances = DistanceData(\n')
			resultFile.write('        distances = {0},\n'.format(entry.data.distances))
			resultFile.write('        method = "{0!s}",\n'.format(entry.data.method))
			resultFile.write('    ),\n')
			resultFile.write('    reference = None,\n')
			resultFile.write('    referenceType = "",\n')
			resultFile.write('    rank = 3,\n')
			resultFile.write('    shortDesc = u"""{0!s}""",\n'.format(entry.shortDesc))
			resultFile.write('    longDesc = \nu"""{0!s}\n""",\n)\n\n'.format(entry.longDesc))

def saveDictionary(specDict):
	"""
	Save the species dictionary.
	"""
	
	with open('dictionary.txt', 'w') as resultFile:
		for name, adjlist in specDict.iteritems():
			resultFile.write('{0}\n'.format(name))
			resultFile.write('{0}\n\n'.format(adjlist))
			
# global_context = None
# local_context = None

# # Set up global and local context
# if global_context is None: global_context = {}
# global_context['__builtins__'] = None
# global_context['True'] = True
# global_context['False'] = False
# if local_context is None: local_context = {}
# local_context['__builtins__'] = None
# local_context['entry'] = loadEntry
# local_context['tree'] = self.__loadTree
# local_context['name'] = self.name
# local_context['shortDesc'] = self.shortDesc
# local_context['longDesc'] = self.longDesc
# local_context['recommended'] = False
# add in anything from the Class level dictionary.
# for key, value in Database.local_context.iteritems():
#     local_context[key]=value

entries = {}
reactionList = []
filePath = os.path.abspath(os.path.join(os.getenv('RMGpy'),'../RMG-database/input/kinetics/families/H_Abstraction/TS_training.py'))

with open(filePath) as resultFile:
	global_context = { '__builtins__': None }
	local_context = {
		'__builtins__': None,
		'True': True,
		'False': False,
		'entry': loadTraining,
		'DistanceData': DistanceData,
		'array': numpy.array,
		'int32': numpy.int32,
	}
	exec resultFile in global_context, local_context

# for i in range(1507):
#     newpath = os.path.join('QMfiles', str(i+1))
#     if os.path.exists(newpath):
#         for files in os.listdir(newpath):
#             if files.endswith('.data'):
#                 with open(os.path.join(newpath, files)) as resultFile:
#                     global_context = { '__builtins__': None }
#                     local_context = {
#                         '__builtins__': None,
#                         'True': True,
#                         'False': False,
#                         'entry': loadEntry,
#                         'DistanceData': DistanceData,
#                         'array': numpy.array,
#                         'int32': numpy.int32,
#                     }
#                     exec resultFile in global_context, local_context
specDict = createLabels(entries)
saveEntries(entries)
saveDictionary(specDict)