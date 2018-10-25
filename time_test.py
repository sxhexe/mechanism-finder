import reactionroute as rr
from datetime import datetime

reactants = ['CCO', 'CCCO', 'CCCCO', 'CCCCCO']
products = ['C=C.O', 'CC=C.O', 'CCC=C.O', 'CCCC=C.O']
reactant = 'CCO'
product = 'C=C.O'

for i in range(5):
    job = rr.ReactionRoute(inputJson='{{"reactant": "{}", "product": "{}"}}'.format(reactant, product))
    startTime = datetime.now()
    head, target = job.isomerSearch()
    endTime = datetime.now()
    print('time used: {}\nnumber of nodes: {}'.format(endTime - startTime, len(job._reactionMap)))
    reactant = 'C' + reactant
    product = 'C' + product