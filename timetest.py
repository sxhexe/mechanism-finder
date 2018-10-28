import reactionroute as rr
from datetime import datetime
import os

reactants = ['CCO', 'CCCO', 'CCCCO', 'CCCCCO']
products = ['C=C.O', 'CC=C.O', 'CCC=C.O', 'CCCC=C.O']
reactant = 'C=CCOC(=CC)C'
product = 'CC(C(=O)C)CC=C'

with open('timetest_{}.log'.format(reactant), 'w') as logfile:
    for i in range(1):
        job = rr.ReactionRoute(inputJson='{{"reactant": "{}", "product": "{}"}}'.format(reactant, product))
        startTime = datetime.now()
        head, target = job.isomerSearch()
        endTime = datetime.now()
    	logfile.write('time used: {}\nnumber of nodes: {}'.format((endTime - startTime).total_seconds(), len(job._reactionMap)))
        reactant = 'C' + reactant
        product = 'C' + product
