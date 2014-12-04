"""
Loads the groups.py file and clears the kinetics, leaving a blank DistanceData object.

Then it prints these to a TS_group.py file.
"""

import os
import re
import logging
import numpy

from rmgpy.molecule import Molecule, Group
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.transitionstates import DistanceData, Entry
from rmgpy.data.base import makeLogicNode

# def makeLogicNode(string):
#     """
#     Creates and returns a node in the tree which is a logic node.
# 
#     String should be of the form:
# 
#     * OR{}
#     * AND{}
#     * NOT OR{}
#     * NOT AND{}
# 
#     And the returned object will be of class LogicOr or LogicAnd
#     """
# 
#     match = re.match("(?i)\s*(NOT\s)?\s*(OR|AND|UNION)\s*(\{.*\})",string)  # the (?i) makes it case-insensitive
#     if not match:
#         raise Exception("Unexpected string for Logic Node: {0}".format(string))
# 
#     if match.group(1): invert = True
#     else: invert = False
# 
#     logic = match.group(2)  # OR or AND (or Union)
# 
#     contents = match.group(3).strip()
#     while contents.startswith('{'):
#         if not contents.endswith('}'):
#             raise Exception("Unbalanced braces in Logic Node: {0}".format(string))
#         contents = contents[1:-1]
# 
#     items=[]
#     chars=[]
#     brace_depth = 0
#     for character in contents:
#         if character == '{':
#             brace_depth += 1
#         if character == '}':
#             brace_depth -= 1
#         if character == ',' and brace_depth == 0:
#             items.append(''.join(chars).lstrip().rstrip() )
#             chars = []
#         else:
#             chars.append(character)
#     if chars: # add last item
#         items.append(''.join(chars).lstrip().rstrip() )
#     if brace_depth != 0: raise Exception("Unbalanced braces in Logic Node: {0}".format(string))
# 
#     if logic.upper() in ['OR', 'UNION']:
#         return LogicOr(items, invert)
#     if logic == 'AND':
#         return LogicAnd(items, invert)
# 
#     raise Exception("Could not create Logic Node from {0}".format(string))

def loadGroup(index, label, group, kinetics, reference=None, referenceType='', shortDesc='', longDesc=''):
    	if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
    		item = makeLogicNode(group)
    	else:
    		item = Group().fromAdjacencyList(group)
    		
    	entry = Entry(
    		index = index,
    		label = label,
    		item = item,
    		data = kinetics,
    		reference = reference,
    		referenceType = referenceType,
    		shortDesc = shortDesc,
    		longDesc = longDesc.strip(),
    		)
    		
    	entries.append(entry)

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
			  	 rank=None
			  	 ):

	entry = Entry(
		index = index,
		label = label,
		item = group,
		data = kinetics,
		reference = reference,
		referenceType = referenceType,
		shortDesc = shortDesc,
		longDesc = longDesc.strip(),
	)
	
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
			  rank=None
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

def saveEntries(entryList):
	"""
	Save the entries.
	"""
	with open('TS_groups.py', 'w') as resultFile:
		resultFile.write('#!/usr/bin/env python\n')
		resultFile.write('# encoding: utf-8\n\n')
		resultFile.write('name = "H_Abstraction/TS_groups"\n')
		resultFile.write('shortDesc = u""\n')
		resultFile.write('longDesc = u"""\n\n"""\n\n')
		for entry in entryList:
			resultFile.write('entry(\n')
			resultFile.write('    index = {0},\n'.format(entry.index))
			resultFile.write('    label = "{0}",\n'.format(entry.label))
			if type(entry.item) == Group:
				resultFile.write('    group =\n"""\n{0}""",\n'.format(entry.item.toAdjacencyList()))
			else:
				resultFile.write('    group = "{0}'.format(entry.item.symbol))
				resultFile.write('{')
				for i, item in enumerate(entry.item.components):
					resultFile.write(item)
					if i == len(entry.item.components) - 1:
						resultFile.write('}",\n')
					else:
						resultFile.write(', ')
			resultFile.write('    distances = DistanceData(distances = {}),\n')
			resultFile.write('    reference = None,\n')
			resultFile.write('    referenceType = "",\n')
			resultFile.write('    shortDesc = u"""""",\n')
			resultFile.write('    longDesc = \nu"""\n""",\n')
			resultFile.write(')\n\n')

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

entries = []
filePath = os.path.join(os.getenv('HOME'),'Code/RMG-database/input/kinetics/families/H_Abstraction/groups.py')

with open(filePath) as resultFile:
	global_context = { '__builtins__': None }
	local_context = {
		'__builtins__': None,
		'True': True,
		'False': False,
		'entry': loadGroup,
		'array': numpy.array,
		'int32': numpy.int32,
	}
	exec resultFile in global_context, local_context

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

saveEntries(entries)
