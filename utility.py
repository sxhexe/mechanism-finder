import pybel
import openbabel as ob
import logging
import numpy as np
from collections import deque
ALLOWED_COORD = {(1,-1):[],
                 (1,0):[1],
                 (1,1):[0],
                 (3,0):[3],
                 (3,1):[],
                 (4,0):[2],
                 (5,-1):[2],
                 (5,0):[3],
                 (5,1):[4],
                 (6,-1):[3],
                 (6,0):[4], # there is a bug in smiles about carbene, so we are not allowing carbene here.
                 (6,1):[3],
                 (7,-1):[2],
                 (7,0):[3],
                 (7,1):[4],
                 (8,-1):[1],
                 (8,0):[2],
                 (8,1):[3],
                 (9,-1):[0],
                 (9,0):[1],
                 (11,1):[],
                 (12,0):[2],
                 (13,0):[3],
                 (14,0):[4],
                 (15,0):[3,5],
                 (15,1):[4],
                 (16,0):[2,3],
                 (17,-1):[0],
                 (17,0):[1],
                 (17,1):[],
                 (27,0):[7,8,9],
                 (35,0):[1],
                 (35,-1):[0],
                 (35,1):[2],
                 (58,0):[1,2,3],
                 (58,-1):[1,2,3]}

class EnergyReadingError(Exception):
    def __init__(self, value):
        self.message = value

    def __str__(self):
        return repr(self.message)


class SmilesError(Exception):
    def __init__(self, value):
        self.message = value

    def __str__(self):
        return repr(self.message)


def printMol(mol,fileFormat = "gjf", keywords = None, printOut = False):
    conv = ob.OBConversion()
    conv.SetOutFormat(fileFormat)
    if printOut:
        logging.info("printing the molecule")
        logging.info(conv.WriteString(mol, True))
    if fileFormat == 'gjf':
        if keywords is not None:
            conv.AddOption("k", ob.OBConversion.OUTOPTIONS, keywords)
        else:
            conv.AddOption("b", ob.OBConversion.OUTOPTIONS)
    elif fileFormat == 'svg':
        conv.AddOption("C", ob.OBConversion.OUTOPTIONS)
        tmpMol = ob.OBMol(mol)
        # tmpMol.DeleteHydrogens()
        return conv.WriteString(tmpMol, True)
    else:
        raise Exception("file format not supported in printMol")
    return conv.WriteString(mol, True)

def molToPngFile(mol, filename):
    mol = pybel.Molecule(mol)
    mol.draw(False, filename)

def printAtom(atom):
    logging.debug('atom index {}, atomic number {}'.format(atom.GetIdx(), atom.GetAtomicNum()))

def printBond(bond):
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()
    logging.debug('bond index {} - {}, atomic number {} - {}'.format(atom1.GetIdx(),
        atom2.GetIdx(), atom1.GetAtomicNum(), atom2.GetAtomicNum()))

def getCanonicalSmiles(mol, noChirality=True):
    conv = ob.OBConversion()
    conv.SetOutFormat("can")
    if noChirality:
        mol.SetDimension(2)
    result = conv.WriteString(mol, True)
    return '.'.join(sorted(result.split('.')))

def strToMol(type, s):
    if s is None:
        print("string is None in strToMol")
    conv = ob.OBConversion()
    conv.SetInFormat(type)
    mol = ob.OBMol()
    success = conv.ReadString(mol,s)
    if success:
        return mol
    else:
        logging.error("converting failure from {} to molecule".format(type))
        raise SmilesError("Failed to convert {} to molecule".format(type))

def smilesToFilename(smiles):
    fileName = ''
    for c in smiles:
        if c == '/':
            fileName += 'z'
            continue
        if c == '\\':
            fileName += 'x'
            continue
        if c == '#':
            fileName += '^'
            continue
        # if c == '(' or c == ')':
        #     fileName += '\\'
        fileName += c
    return fileName

def smilesToSysCall(smiles):
    fileName = smilesToFilename(smiles)
    call = ''
    for c in fileName:
        if c == '(' or c == ')' or c == '$':
            call += '\\'
        call += c
    return call

def numValenceElectron(atomicNumber):
    if atomicNumber <= 2:
        return atomicNumber
    elif atomicNumber <= 10:
        return atomicNumber - 2
    elif atomicNumber <= 18:
        return atomicNumber - 10
    elif atomicNumber <= 30:
        return atomicNumber - 18
    elif atomicNumber <= 36:
        return atomicNumber - 28
    elif atomicNumber <= 48:
        return atomicNumber - 36
    elif atomicNumber <= 54:
        return atomicNumber - 46
    elif atomicNumber == 58:
        return 10
    else:
        print('Atomic number not supported in calculating the number of valence electrons. Either it is from the 6th row and below or it is an invalid number, setting it to 10')
        return 10

def atomTotalBondOrder(atom):
    nBonds = 0
    for bond in ob.OBAtomBondIter(atom):
        nBonds += bond.GetBondOrder()
    return nBonds

def molToMat(mol):
    n = mol.NumAtoms()
    mat = np.array([[0 for _i in range(n+1)] for _j in range(n+3)])
    for i in range(1, n+1):
        mat[i][0] = mol.GetAtom(i).GetAtomicNum()
        mat[0][i] = mat[i][0]
    for atom in ob.OBMolAtomIter(mol):
        i = atom.GetIdx()
        nBonds = 0
        for bond in ob.OBAtomBondIter(atom):
            nBonds += bond.GetBondOrder()
        nonBondingElecs = numValenceElectron(atom.GetAtomicNum()) - nBonds - atom.GetFormalCharge()
        mat[i][i] = nonBondingElecs
        mat[n+1][i] = nBonds
        mat[n+2][i] = atom.GetFormalCharge()
    for bond in ob.OBMolBondIter(mol):
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        mat[i][j] = bond.GetBondOrder()
        mat[j][i] = mat[i][j]
    return mat

def matToMol(mat):
    n = len(mat) - 3
    mol = ob.OBMol()
    mol.BeginModify()
    for i in range(1, n+1):
        mol.NewAtom(i)
        atom = mol.GetAtom(i)
        atom.SetAtomicNum(mat[i][0])
        atom.SetFormalCharge(mat[n+2][i])
    for i in range(1, n+1):
        for j in range(1, i):
            if mat[i][j] != 0:
                mol.AddBond(i, j, mat[i][j])
    return mol

def separateFragments(pymol, dist=3.0):
    mol = pymol.OBMol
    nAtom = len(pymol.atoms)
    unvisited = set(range(1, nAtom+1))
    q = deque([1])
    fragments = [set()]
    while q or unvisited:
        if q:
            curr = q.popleft()
            unvisited.remove(curr)
            fragments[-1].add(curr)
        else:
            curr = unvisited.pop()
            fragments.append({curr})
        atom = mol.GetAtom(curr)
        for nbr in ob.OBAtomAtomIter(atom):
            nbrNum = nbr.GetIdx()
            if nbrNum in unvisited:
                q.append(nbrNum)
    coords = []
    for atom in pymol:
        coords.append(list(atom.coords))
    nFragments = len(fragments)
    delta = -(nFragments-1) * dist/2.0
    for fragment in fragments:
        for atomIdx in fragment:
            x, y, z = coords[atomIdx-1]
            coords[atomIdx-1] = [x + delta, y + delta, z + delta]
        delta += dist
    coords = [item for sublist in coords for item in sublist ]
    c_coords = ob.double_array(coords)
    mol.SetCoordinates(c_coords)
