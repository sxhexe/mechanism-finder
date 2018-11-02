import reactionroute as rr
from datetime import datetime
import os

reactants = ['CCO', 'CCCO', 'CCCCO', 'CCCCCO']
products = ['C=C.O', 'CC=C.O', 'CCC=C.O', 'CCCC=C.O']
# reactant = 'C=CCOC(=CC)C'
# product = 'CC(C(=O)C)CC=C'
# activeList = '1,2,3,4,5,6'
reactant = 'CC=O'
product = 'C=CO'
#reactant = 'CCO'
#product = 'C=C.O'
#reactant = 'C=C.Br'
#product = 'CCBr'
allowPair = 2

with open('timetest_{}_{}pairs.log'.format(reactant, allowPair), 'w') as logfile:
    for i in range(10):
        # job = rr.ReactionRoute(inputJson='{{"reactant": "{}", "product": "{}", "activeList": [{}]}}'.format(reactant, product, activeList))
        job = rr.ReactionRoute(inputJson='{{"reactant": "{}", "product": "{}"}}'.format(reactant, product))
	job._allowedPairs = allowPair
        startTime = datetime.now()
        head, target = job.isomerSearch()
        endTime = datetime.now()
        # logfile.write('time used: {}, number of nodes: {}, activeList: {}\n'.format((endTime - startTime).total_seconds(), len(job._reactionMap), activeList))
	sourceTime = job.eSourceSizeRunningSum
	targetTime = job.eTargetSizeRunningSum
	logfile.write('source size running average: {}, target size running average: {}'.format(float(sourceTime[0])/sourceTime[1], float(targetTime[0])/targetTime[1]))
        logfile.write('time used: {}, number of nodes: {}, addMolCount: {}\n'.format((endTime - startTime).total_seconds(), len(job._reactionMap), job.addMolCount))
        reactant = 'C' + reactant
        product = 'C' + product
        # activeList += ',{},{}'.format(7+2*i, 7+2*i+1)
