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

def calculate(TS, count):
    quantumMechanics = QMCalculator()
    quantumMechanics.settings.software = 'gaussian'
    quantumMechanics.settings.fileStore = 'QMfiles'
    quantumMechanics.settings.scratchDirectory = 'scratch'
    quantumMechanics.settings.onlyCyclics = False
    quantumMechanics.settings.maxRadicalNumber = 0
    
    reaction = Reaction(label='H_Abstraction', reactants=reactant.split(), products=product.split(), reversible=True)
    
    qmReaction = GaussianTSB3LYP(reaction, quantumMechanics.settings)
    
    qmReaction.generateGeometry()
    


for count, TS in enumerate(tsStructures):
    calculate(TS, count)
    count += 1