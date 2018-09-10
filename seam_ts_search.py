import pybel
import openbabel as ob
import numpy as np
import logging

logger = logging.getLogger(__name__)

class TsEstimConvergeError(Exception):
    def __init__(self, value):
        self.message = value
    def __str__(self):
        return repr(self.message)

def getGradients(ff, mol):
    gList = []
    for atom in mol:
        gList.append(ff.GetGradient(atom.OBAtom).GetX())
        gList.append(ff.GetGradient(atom.OBAtom).GetY())
        gList.append(ff.GetGradient(atom.OBAtom).GetZ())
    g = np.array(gList) # keep in mind that this gradient from openbabel is negative gradient really.
    return g

def updateCoords(mol, coords):
    molStringLines = mol.write('mol').split('\n')
    for i in range(len(mol.atoms)):
        molStringLines[i+4] = "{:10.4f}{:10.4f}{:10.4f}".format(coords[i*3], coords[i*3+1], coords[i*3+2]) + molStringLines[i+4][30:]
    mol2 = pybel.readstring('mol','\n'.join(molStringLines))
    return mol2

def getCoordinates(mol):
    cList = []
    xyz = mol.write('xyz')
    lines = xyz.split('\n')
    for line in lines[2:-1]:
        for num in line.split()[1:4]:
            cList.append(float(num))
    return np.array(cList)

def SeamTsSearch(mol1, mol2, forcefield='uff'):

    def getFAndG(molA, molB):
        logger.debug("in getFAndG")
        ff1 = ob.OBForceField.FindForceField(forcefield)
        ff2 = ob.OBForceField.FindForceField(forcefield)
        ff1.Setup(molA.OBMol) # Have to setup the molecule everytime before calculating energy/gradients
        e1 = ff1.Energy() # Energy() call is needed to generate gradients
        if e1 != e1:
            raise ValueError('e1 is nan')
        logger.debug("e1 = {}".format(e1))
        g1 = getGradients(ff1, molA)
        logger.debug("g1 = {}".format(g1))
        ff2.Setup(molB.OBMol)
        e2 = ff2.Energy()
        # logger.debug("e2 = ", e2)
        g2 = getGradients(ff2, molB)
        # logger.debug("g2 = ", g2)
        alpha = 20 # arbitrary constant. Bigger alpha means smaller weight on energy difference
        # logger.debug("(g1+g2)/2= ", (g1+g2)/2)
        # logger.debug("(2/alpha)*(e1-e2)*(g1-g2) = ", (2.0/alpha)*(e1-e2)*(g1-g2))
        return (e1+e2)/2+(e1-e2)*(e1-e2)/alpha, (g1+g2)/2+(2.0/alpha)*(e1-e2)*(g1-g2)

    def lineSearch(direction, coords, f):
        # ref: https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
        #      https://en.wikipedia.org/wiki/Backtracking_line_search
        logger.debug("\n**********in lineSearch")
        step = 0.05 # maximum step number, shrink until energy are sufficiently decreased.
        # logger.debug("direction = ", direction)
        directionNorm = np.linalg.norm(direction)
        directionUnit = direction/directionNorm
        for iteration in range(200):
            newCoords = coords + step * directionUnit
            # logger.debug("newCoords = ", newCoords)
            fNew = getFAndG(updateCoords(mol1, newCoords), updateCoords(mol2, newCoords))[0]
            # logger.debug("fNew and f - step*0.5 are ", fNew, f - step*0.5)
            if fNew < f - step*0.5:
                return newCoords
            else:
                step = step*0.5
        else:
            raise TsEstimConvergeError('Error! line search in SeamTsSearch did not converge in 200 steps')

    def getBeta(directionNew, direction):
        return directionNew.dot(directionNew.T)/direction.dot(direction.T)

    logger.info("before iterations...")
    f, g = getFAndG(mol1, mol2)
    # logger.debug("f and g are ", f, g)
    directionNew = g # g from openbabel is direction (-g)
    conjugateDirectionNew = directionNew
    c = getCoordinates(mol1)
    print("coordinates \n", c)
    cNew = lineSearch(directionNew, c, f)
    c = cNew
    conjugateDirection = conjugateDirectionNew
    direction = directionNew
    mol1 = updateCoords(mol1, c)
    mol2 = updateCoords(mol2, c)
    for iteration in xrange(50):
        # print("\n\n============calculating iteration number {}".format(iteration))
        fNew, gNew = getFAndG(mol1, mol2)
        # logger.debug("fNew and gNew = ", fNew, gNew)
        # logger.info("f and fNew = ", f, fNew)
        if f - fNew < 1:
            break
        directionNew = gNew
        beta = getBeta(directionNew, direction)
        # logger.debug("beta = ", beta)
        conjugateDirectionNew = directionNew + beta * conjugateDirection
        # logger.debug("conjugateDirectionNew = ", conjugateDirectionNew)
        cNew = lineSearch(conjugateDirectionNew, c, fNew)
        c = cNew
        direction = directionNew
        f = fNew
        g = gNew
        conjugateDirection = conjugateDirectionNew
        mol1 = updateCoords(mol1, c)
        mol2 = updateCoords(mol2, c)
    else:
        raise TsEstimConvergeError('Error! SeamTsSearch did not converge in 50 steps')
    return mol1

def main():
    import pybel
    import openbabel as ob
    import sys
    h2o_ts_string = '\n OpenBabel07171718443D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.5111    1.1929    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0943    1.2466   -0.0910 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1878    0.5054    0.6031 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\nM  END\n'
    h2o_ts = pybel.readstring('mol',h2o_ts_string)
    h2o_ts_b = pybel.readstring('sdf', h2o_ts.write('sdf'))
    h2o_ts_b.OBMol.DeleteBond(h2o_ts_b.OBMol.GetBond(1,2))
    h2o_qst3 = SeamTsSearch(h2o_ts, h2o_ts_b, sys.argv[1])
    f = open('h2o_ts.com', 'w')
    f.write(h2o_ts.write('gjf'))
    f.close()
    f_b = open("h2o_ts_b.com", 'w')
    f_b.write(h2o_ts_b.write('gjf'))
    f_b.close()
    f_qst3 = open('h2o_qst3.com', 'w')
    h2o_qst3.OBMol.SetTitle('seam_ts_search generated')
    f_qst3.write(h2o_qst3.write('gjf'))
    f_qst3.close()




if __name__ == "__main__":
    main()
