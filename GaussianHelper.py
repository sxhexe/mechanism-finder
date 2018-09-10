import pybel
import openbabel as ob

def errorOrNormal(fileName):
    import re
    f = open(fileName, 'r')
    f.seek(-740, 2)
    lines = f.read().split('\n')
    for line in lines:
        if "Error termination via" in line:
            m = re.search(r"l([0-9]+)\.exe", line)
            if m is None:
                return line
            return "Error via link "+m.group(1)
        elif "Normal termination" in line:
            return "Normal"

def logParser(logFile):
    # Assume that each log file only has one job, i.e., no --link1--.
    try:
        mol = pybel.readfile('log', logFile).next()
    except StopIteration:
        print 'pybel failed to read molecule from {}'.format(logFile)
        return {'result':'Error'}
    moldict = {
        'formula':mol.formula,
        'energy':mol.energy,
        'job specs':mol.data['Comment'], # or self._gaussianKeywords
        'result':errorOrNormal(logFile),
        'smiles':mol.write('can').split('\t')[0]
        }
    return moldict
    dataIter = ob.vectorpOBGenericData(mol.OBMol.GetData())