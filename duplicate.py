"""
Run after `database.py`

Loads the TS_training.py file, then creates the reverse reactions and addes them to the list.the .data files generated
These are printed to a new TS_training_duplicates.py in the folder where this is run.
"""

import os
import logging
import numpy
from copy import deepcopy

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
	for testReaction in reactionList:
		if testReaction.isIsomorphic(reaction):
			alreadyHave = True
	
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

def duplicate(item, family):
	
	newItem = deepcopy(item)
	newDist = deepcopy(newItem.data.distances)
	if family.lower()=='h_abstraction':
		newDist['d12'] = item.data.distances['d23']
		newDist['d23'] = item.data.distances['d12']
	elif family.lower()=='intra_h_migration':
		newDist['d13'] = item.data.distances['d23']
		newDist['d23'] = item.data.distances['d13']
	
	newItem.data.distances = newDist
	
	newProd = deepcopy(item.item.reactants)
	newReact = deepcopy(item.item.products)
	reaction = Reaction(reactants=newReact, products=newProd)
		
	entry = Entry(
		index = len(entries) + 1,
		label = item.label,
		item = reaction,
		data = newItem.data,
		reference = item.reference,
		referenceType = item.referenceType,
		shortDesc = 'Reverse reaction for reaction index {0}'.format(item.index),
		longDesc = item.longDesc.strip(),
		rank = item.rank,
	)
	
	entries['{0:d}:{1}'.format(len(entries) + 1,item.label)] = entry
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

def saveEntries(entryList, family):
	"""
	Save the entries.
	"""
	with open('TS_training_duplicate.py', 'w') as resultFile:
		resultFile.write('#!/usr/bin/env python\n')
		resultFile.write('# encoding: utf-8\n\n')
		resultFile.write('name = "{0}/TS_training"\n'.format(family))
		resultFile.write('shortDesc = u"Distances used to train group additivity values for TS geometries"\n')
		resultFile.write('longDesc = u"""\nPut interatomic distances for reactions to use as a training set for fitting\ngroup additivity values in this file.\n"""\n')
		resultFile.write('recommended = True\n\n')
		for i in range(len(entryList)):
			entry = entryList[str(i+1)+':']
			resultFile.write('entry(\n')
			resultFile.write('    index = {0},\n'.format(entry.index))
			resultFile.write('    reactant1 = """\n{0!s}""",\n'.format(entry.item.reactants[0].toAdjacencyList()))
			if len(entry.item.reactants)==2:
				resultFile.write('    reactant2 = """\n{0!s}""",\n'.format(entry.item.reactants[1].toAdjacencyList()))
			resultFile.write('    product1 = """\n{0!s}""",\n'.format(entry.item.products[0].toAdjacencyList()))
			if len(entry.item.products)==2:
				resultFile.write('    product2 = """\n{0!s}""",\n'.format(entry.item.products[1].toAdjacencyList()))
			resultFile.write('    distances = DistanceData(\n')
			resultFile.write('        distances = {0},\n'.format(entry.data.distances))
			resultFile.write('        method = "{0!s}",\n'.format(entry.data.method))
			resultFile.write('    ),\n')
			resultFile.write('    reference = None,\n')
			resultFile.write('    referenceType = "",\n')
			resultFile.write('    rank = 3,\n')
			resultFile.write('    shortDesc = u"""{0!s}""",\n'.format(entry.shortDesc))
			resultFile.write('    longDesc = \nu"""{0!s}\n""",\n)\n\n'.format(entry.longDesc))

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
filePath = 'TS_training.py'
family = None
with open(filePath, 'r') as trainingFile:
	trainingLines = trainingFile.readlines()
	for line in trainingLines:
		if line.startswith('name = '):
			family = line.split('/')[0].split('"')[1]
			break


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


for i in range(len(entries)):
	item = entries[str(i+1) + ':']
	duplicate(item, family)
	
# for i in range(334):
#     newpath = os.path.join('QMfiles', str(i+1))
#     for files in os.listdir(newpath):
#         if files.endswith('.data'):
#             with open(os.path.join(newpath, files)) as resultFile:
#                 global_context = { '__builtins__': None }
#                 local_context = {
#                     '__builtins__': None,
#                     'True': True,
#                     'False': False,
#                     'entry': loadEntry,
#                     'DistanceData': DistanceData,
#                     'array': numpy.array,
#                     'int32': numpy.int32,
#                 }
#                 exec resultFile in global_context, local_context

saveEntries(entries, family)
			