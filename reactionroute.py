import argparse
import itertools
import os
import sqlite3
import subprocess
import sys
import time
from collections import deque, defaultdict

from GaussianHelper import *
from utility import *

# from seam_ts_search import *

np.set_printoptions(threshold=np.inf, linewidth=150)


class ReactionGraphEdge:
    def __init__(self, fromNode, node, eSources, eTargets):
        self.fromNode = fromNode
        self.node = node
        self.eSources = list(eSources)
        self.eTargets = list(eTargets)
        self.ts = None
        self.tsEnergy = 0.0
        self.onPath = False


class ReactionGraphNode:
    def __init__(self, mol=None, smiles=None, depth=None):
        if mol is not None:
            self.mol = ob.OBMol(mol)
            self.smiles = getCanonicalSmiles(mol)
        elif smiles is not None:
            self.mol = strToMol('smi', smiles)
            self.smiles = smiles
        else:
            logging.warning("a molecule is needed to create a ReactionGraphNode")
            sys.exit()
        self.neighbors = {}
        self.depths = []
        self.energy = 0.0
        self.onPath = False
        if depth is not None:
            self.depths.append(depth)


class ReactionRoute:
    def __init__(self, reactantString=None, productString=None, inputJson=None):
        # An equality holds for all atom: total bond order + # of non-bonding electrons = # of valence electrons + formal charge
        # This gives rise to the following rules for each atom. Atomic number determines the number of valence electrons, formal charge is formal charge, then there is total bond order. Once these three is fixed, the Luis structure is determined.
        self._allowedCoordNum = ALLOWED_COORD
        self._minFc = {}
        for pair, tboList in self._allowedCoordNum.items():
            if tboList != []:
                if pair[0] not in self._minFc:
                    self._minFc[pair[0]] = pair[1]
                else:
                    self._minFc[pair[0]] = min(self._minFc[pair[0]], pair[1])

        self._outputLevel = 2
        self._maxStep = 3
        self._maxExtraStep = 1
        self._doCalculation = False
        self._structureScreen = True
        self._energyScreen = True
        self._intermediateThresh = 200.0
        self._gaussianKeywords = "# pm6 3-21g opt"
        self._doTsSearch = False
        self._tsThresh = 200.0
        self._gaussianTsKeywords = '# pm6 3-21g opt=(ts,noeigen,calcfc,maxcyc=100)'
        self._energyBaseLine = 0.0
        self.ignoreList = set()
        self.activeList = set()
        self._invalidStructures = []
        self._reactantString = reactantString
        self._productString = productString
        self._targetLeastStep = 100
        self._targetFound = False
        self._reactionMap = {}
        self._energyMap = {}
        self._fragmentEnergyMap = {}
        self._gsub = False
        self._save = True
        self._pathOnly = True
        self._preEnergyScreen = False
        self._matrixForm = True
        self._filterFc = True
        self._noProduct = False
        self._allowedPairs = 2

        self.picDir = "pics"
        self.molPicFormat = 'svg'
        self.jobName = "tmpName"
        if inputJson is not None:
            self.inputJson(inputJson)

    def inputJson(self, inputJson):
        import json
        params = json.loads(inputJson)
        if 'reactant' in params:
            self._reactantString = params['reactant']
        if 'product' in params:
            self._productString = params['product']
        if 'maxStep' in params:
            self._maxStep = params['maxStep']
        if 'maxExtraStep' in params:
            self._maxExtraStep = params['maxExtraStep']
        if 'doCalculation' in params:
            self._doCalculation = params['doCalculation']
        if 'structureScreen' in params:
            self._structureScreen = params['structureScreen']
        if 'energyScreen' in params:
            self._energyScreen = params['energyScreen']
        if 'intermediateThresh' in params:
            self._intermediateThresh = params['intermediateThresh']
        if 'gaussianKeywords' in params:
            self._gaussianKeywords = params['gaussianKeywords']
        if 'doTsSearch' in params:
            self._doTsSearch = params['doTsSearch']
        if 'tsThresh' in params:
            self._tsThresh = params['tsThresh']
        if 'gaussianTsKeywords' in params:
            self._gaussianTsKeywords = params['gaussianTsKeywords']
        if 'ignoreList' in params:
            self.ignoreList = set(params['ignoreList'])
        if 'activeList' in params:
            self.activeList = set(params['activeList'])
        if 'gsub' in params:
            self._gsub = params['gsub']
        if 'outputLevel' in params:
            self._outputLevel = params['outputLevel']
        if 'save' in params:
            self._save = params['save']
        if 'pathOnly' in params:
            self._pathOnly = params['pathOnly']
        if 'preEnergyScreen' in params:
            self._preEnergyScreen = params['preEnergyScreen']
        if 'matrixForm' in params:
            self._matrixForm = params['matrixForm']
        if 'invalidStructures' in params:
            self._invalidStructures = params['invalidStructures']
        if 'allowedPairs' in params:
            self._allowedPairs = params['allowedPairs']
        if 'jobName' in params:
            self.jobName = params['jobName']

    def canBreakOrFormBond(self, atom, breakOrForm, nElec):
        # Decide if an atom can break or form bond in a certain way (get or lose certain number of electrons)
        formalCharge = atom.GetFormalCharge()
        atomicNum = atom.GetAtomicNum()
        nBonds = 0
        for bond in ob.OBAtomBondIter(atom):
            nBonds += bond.GetBondOrder()
        if breakOrForm.lower() == "break":
            nBondChange = -1
        elif breakOrForm.lower() == "form":
            nBondChange = 1
        else:
            raise Exception("breakOrForm can either be break or form")
        try:
            if nBonds + nBondChange in self._allowedCoordNum[(atomicNum, formalCharge + nBondChange * (nElec - 1))]:
                return 1
            else:
                return 0
        except KeyError:
            return 0

    def checkLuisRule(self, *args, **kwargs):
        for arg in args:
            if type(arg) is int:
                atom = kwargs['mol'].GetAtom(arg)
                pair = (atom.GetAtomicNum(), atom.GetFormalCharge())
                if pair not in self._allowedCoordNum or atomTotalBondOrder(atom) not in self._allowedCoordNum[pair]:
                    return False
            elif type(arg) is tuple:
                if not self.checkLuisRule(*arg, mol=kwargs['mol']):
                    return False
            elif type(arg) is ob.OBMol:
                for atom in ob.OBMolAtomIter(arg):
                    if not self.checkLuisRule(atom, mol=arg):
                        return False
            elif type(arg) is ob.OBBond:
                if not self.checkLuisRule(arg.GetBeginAtom()) or not self.checkLuisRule(arg.GetEndAtom()):
                    return False
            elif type(arg) is ob.OBAtom:
                pair = (arg.GetAtomicNum(), arg.GetFormalCharge())
                if pair not in self._allowedCoordNum or atomTotalBondOrder(atom) not in self._allowedCoordNum[pair]:
                    return False

        return True

    def obeyLuisRule(self, atom, nBondChange, nElectronChange):
        # Decide if an atom can break or form bond in a certain way (get or lose certain number of electrons)
        formalCharge = atom.GetFormalCharge()
        atomicNum = atom.GetAtomicNum()
        nBonds = atomTotalBondOrder(atom)
        if abs(nElectronChange) == 2:
            formalChargeChange = -nElectronChange / 2
        else:
            raise Exception("one electron reaction not supported yet")
        try:
            if nBonds + nBondChange in self._allowedCoordNum[(atomicNum, formalCharge + formalChargeChange)]:
                return 1
            else:
                return 0
        except KeyError:
            return 0

    def moveElec(self, mol, atom1Idx, atom2Idx, atom3Idx, nElec):
        mol.BeginModify()
        atom1 = None if atom1Idx is None else mol.GetAtom(atom1Idx)
        atom2 = None if atom2Idx is None else mol.GetAtom(atom2Idx)
        atom3 = None if atom3Idx is None else mol.GetAtom(atom3Idx)
        if atom1 is None:  # lone pair (atom2) to bond (atom2 - atom3)
            atom2.SetFormalCharge(atom2.GetFormalCharge() + 1)
            # ob.OBPairData(atom2.GetData('nLonePair')).SetValue(str(int(ob.OBPairData(atom.GetData('nLonePair')).GetValue())-2))
            bond = mol.GetBond(atom2, atom3)
            if bond is None:
                mol.AddBond(atom2.GetIdx(), atom3.GetIdx(), 1)
            else:
                bond.SetBondOrder(bond.GetBondOrder() + 1)
            atom3.SetFormalCharge(atom3.GetFormalCharge() - 1)
        elif atom3 is None:  # bond (atom1 - atom2) to lone pair (atom2)
            bond = mol.GetBond(atom1, atom2)
            bondOrder = bond.GetBondOrder()
            print('bondorder is {}'.format(bondOrder))
            if bondOrder == 1:
                mol.DeleteBond(bond)
            else:
                bond.SetBondOrder(bondOrder - 1)
            atom1.SetFormalCharge(atom1.GetFormalCharge() + 1)
            atom2.SetFormalCharge(atom2.GetFormalCharge() - 1)
            # ob.OBPairData(atom2.GetData('nLonePair')).SetValue(str(int(ob.OBPairData(atom.GetData('nLonePair')).GetValue())+2))
        else:  # bond1 (atom1 - atom2) to bond2 (atom2 - atom3)
            bond1 = mol.GetBond(atom1, atom2)
            bond2 = mol.GetBond(atom2, atom3)
            atom1.SetFormalCharge(atom2.GetFormalCharge() + 1)
            atom3.SetFormalCharge(atom3.GetFormalCharge() - 1)
            bondOrder1 = bond1.GetBondOrder()
            if bondOrder1 == 1:
                mol.DeleteBond(bond1)
            else:
                bond1.SetBondOrder(bondOrder1 - 1)
            if bond2 is None:
                mol.AddBond(atom2.GetIdx(), atom3.GetIdx(), 1)
            else:
                bond2.SetBondOrder(bond2.GetBondOrder() + 1)
        mol.EndModify()
        printMol(mol, printOut=True)
        for bond in ob.OBMolBondIter(mol):
            printBond(bond)

    def changeFormalCharge(self, mol, idx, change):
        atom = mol.GetAtom(idx)
        atom.SetFormalCharge(atom.GetFormalCharge() + change)

    def changeBondOrder(self, mol, i, j, change):
        bond = mol.GetBond(i, j)
        if bond is None:
            mol.AddBond(i, j, 1)
        else:
            bondOrder = bond.GetBondOrder()
            if bondOrder == 1 and change == -1:
                mol.DeleteBond(bond)
            else:
                bond.SetBondOrder(bondOrder + change)

    def oxidize(self, mol, eSource):
        if type(eSource) is int:
            self.changeFormalCharge(mol, eSource, +2)
        elif type(eSource) is tuple:
            i, j = eSource
            self.changeBondOrder(mol, i, j, -1)
            self.changeFormalCharge(mol, i, +1)
            self.changeFormalCharge(mol, j, +1)

    def reduce(self, mol, eTarget):
        if type(eTarget) is int:
            self.changeFormalCharge(mol, eTarget, -2)
        elif type(eTarget) is tuple:
            i, j = eTarget
            self.changeBondOrder(mol, i, j, +1)
            self.changeFormalCharge(mol, i, -1)
            self.changeFormalCharge(mol, j, -1)

    def checkStructure(self, molMat):
        # return True
        for item in self._invalidStructures:
            if item == 'noDoubleSameCharge':
                anion, cation = 0, 0
                for charge in molMat[self.nAtom + 2][1:]:
                    if charge > 0:
                        cation += 1
                    elif charge < 0:
                        anion += 1
                    if anion > 1 or cation > 1:
                        return False
        return True

    def doGaussian(self, mol, fullFileName, smiles):
        conn = sqlite3.connect('reactionroute.db')
        cursor = conn.cursor()
        records = cursor.execute('select energy from jobArchive where smiles == ? and keywords == ?',
                                 (smiles, self._gaussianKeywords))
        logging.debug('database connection established')
        record = records.fetchone()
        if record:
            logging.debug('{} record found in the database'.format(smiles))
            print('{} record found in the database'.format(smiles))
            conn.close()
            return record[0]
        else:
            logging.debug('{} not found in the database, doing calculations...'.format(smiles))
        molCopy = ob.OBMol(mol)
        molCopy.SetTitle("ReactionRoute")
        inputFile = open("gaussian" + os.sep + fullFileName + ".gjf", 'w')
        # inputFile = open(os.path.join("gaussian", fullFileName+".gjf"), 'w')
        op3d = ob.OBOp.FindType("gen3d")
        op3d.Do(molCopy, '3')
        inputFile.write(printMol(molCopy, fileFormat="gjf", keywords=self._gaussianKeywords))
        inputFile.close()
        gaussianCall = ''
        for i, c in enumerate(fullFileName):
            if (c == '(' or c == ')' or c == '$') and i > 0 and fullFileName[i - 1] != '\\':
                gaussianCall += '\\'
            gaussianCall += c
        if self._gsub:
            print("gsub -fastq gaussian" + os.sep + gaussianCall + ".gjf")
            logging.info("gsub -fastq gaussian" + os.sep + gaussianCall + ".gjf")
            output = subprocess.check_output('cd gaussian; gsub -fastq ' + gaussianCall + '.gjf; cd ..', shell=True)
            # print output
            jobId = output.split()[7]
            while True:
                time.sleep(1)
                outputQstat = subprocess.check_output('qstat', shell=True)
                if jobId not in outputQstat:
                    break
        else:
            print("gdv gaussian" + os.sep + gaussianCall + ".gjf")
            logging.info("gdv gaussian" + os.sep + gaussianCall + ".gjf")
            os.system("gdv gaussian" + os.sep + gaussianCall + ".gjf")

        molDict = logParser('gaussian' + os.sep + fullFileName + '.log')
        if molDict['result'] == 'Normal':
            cursor.execute('insert into jobArchive (smiles, keywords, formula, energy) values (?, ?, ?, ?)',
                           (smiles, self._gaussianKeywords, molDict['formula'], molDict['energy']))
            molCopyEnergy = molDict['energy']
        else:
            logging.error("First gaussian run failed. Trying second time with op3d.do(frag, 'dist')")
            inputFile = open("gaussian" + os.sep + fullFileName + ".gjf", 'w')
            op3d = ob.OBOp.FindType("gen3d")
            op3d.Do(molCopy, 'dist')
            inputFile.write(printMol(molCopy, fileFormat="gjf", keywords=self._gaussianKeywords))
            inputFile.close()
            os.system("gdv gaussian" + os.sep + gaussianCall + ".gjf")

            molDict = logParser('gaussian' + os.sep + fullFileName + '.log')
            if molDict['result'] == 'Normal':
                cursor.execute('insert into jobArchive (smiles, keywords, formula, energy) values (?, ?, ?, ?)',
                               (smiles, self._gaussianKeywords, molDict['formula'], molDict['energy']))
                molCopyEnergy = molDict['energy']
            else:
                cursor.execute('insert into jobArchive (smiles, keywords, formula, energy) values (?, ?, ?, ?)',
                               (smiles, self._gaussianKeywords, 'error', -999999999.0))
                logging.error("Second gaussian run failed. ")
                molCopyEnergy = -999999999.0
        conn.commit()
        conn.close()
        return molCopyEnergy

    def computeQMEnergy(self, mol, software, method, fragmentEnergyMap=None):
        if not os.path.isdir(software):
            os.system("mkdir " + software)
        molCopy = strToMol('smi', printMol(mol, 'smi'))
        molCopy.AddHydrogens()
        smiles = getCanonicalSmiles(molCopy)
        fileName = smilesToFilename(smiles)
        if software.lower() == "gaussian" or software.lower() == "gauss":
            logging.debug('using Gaussian to calculate...')
            logging.debug('the molecule that is about to be separated is {}'.format(smiles))
            fragments = molCopy.Separate()
            logging.debug('after separate')
            logging.debug('there are {} fragments'.format(len(fragments)))
            if len(fragments) >= 2:
                energySum = 0.0
                for i, frag in enumerate(fragments):
                    fragSmiles = getCanonicalSmiles(frag)
                    logging.debug('fragment{} = {}'.format(i, fragSmiles))
                    if fragSmiles in fragmentEnergyMap:
                        frag.SetEnergy(fragmentEnergyMap[fragSmiles])
                        logging.debug("this fragment's energy has been calculated. It's %d kcal/mol" % (
                            fragmentEnergyMap[fragSmiles]))
                    else:
                        fragmentEnergy = self.doGaussian(frag, fileName + str(i), fragSmiles)
                        logging.debug("the energy of this fragment is %d kcal/mol" % (fragmentEnergy))
                        frag.SetEnergy(fragmentEnergy)
                        fragmentEnergyMap[fragSmiles] = fragmentEnergy
                    energySum += fragmentEnergyMap[fragSmiles]
                logging.info("The energy of the molecule is %d kcal/mol" % (energySum))
                return energySum
            else:
                return self.doGaussian(molCopy, fileName, smiles)

    def getGaussianEnergy(self, fileName):
        f = open(fileName, 'r')
        outputChunk = ''
        for line in f:
            if '\\' in line:
                outputChunk += line.strip()
        outputList = outputChunk.split('\\')
        for word in outputList:
            if "HF=" in word:
                logging.debug("HF=" + word[3:])
                return 627.5 * float(word[3:])
        raise EnergyReadingError("Can't read energy from gaussian output %s" % (fileName))

    def checkChangeTable(self, molMat, changeTable, tboChange, fcChange):
        for atom in tboChange.keys() + fcChange.keys():
            newTbo = molMat[self.nAtom + 1][atom] + tboChange[atom]
            newFc = molMat[self.nAtom + 2][atom] + fcChange[atom]
            if newTbo not in self._allowedCoordNum.get((molMat[0][atom], newFc), []):
                return False
        return True

    def applyChanges(self, molMat, changeTable, tboChange, fcChange):
        for item in changeTable.items():
            molMat[item[0][0]][item[0][1]] += item[1]
        for item in tboChange.items():
            molMat[self.nAtom + 1][item[0]] += item[1]
        for item in fcChange.items():
            molMat[self.nAtom + 2][item[0]] += item[1]

    def isomerSearch(self):
        reactantMol = strToMol('smi', self._reactantString)
        self._reactantString = getCanonicalSmiles(reactantMol)
        logging.info("reactant = {}".format(self._reactantString))
        reactantMol.AddHydrogens()
        printMol(reactantMol, fileFormat="gjf", printOut=True)
        if not self._noProduct:
            productMol = strToMol('smi', self._productString)
            self._productString = getCanonicalSmiles(productMol)
            logging.info("product = {}".format(self._productString))

        self.nAtom = reactantMol.NumAtoms()

        # List indexes start from 1.
        if self.activeList and not self.ignoreList:
            allset = set(range(1, self.nAtom + 1))
            self.ignoreList = allset - self.activeList
        elif self.ignoreList and not self.activeList:
            allset = set(range(1, self.nAtom + 1))
            self.activeList = allset - self.ignoreList
        elif not self.activeList and not self.ignoreList:
            self.activeList = set(range(1, self.nAtom + 1))

        logging.info("ignoreList = {}".format(self.ignoreList))
        q = deque()
        head = ReactionGraphNode(mol=reactantMol)
        q.append(head)
        self._reactionMap[self._reactantString] = head
        self._energyMap = {self._reactantString: 0.0}
        if self._doCalculation:
            self._energyBaseLine = self.computeQMEnergy(reactantMol, "gaussian", self._gaussianKeywords,
                                                        self._fragmentEnergyMap)
        else:
            self._energyBaseLine = 0.0
        head.energy = 0.0
        nStep = -1
        while q:  # start Breadth-First-Search
            qSize = len(q)
            nStep += 1
            logging.info("=========================================================")
            logging.info("                     nStep = " + str(nStep))
            if nStep >= self._maxStep or nStep > self._targetLeastStep + self._maxExtraStep:
                logging.info("step number {}, exceeding maximum step {}".format(nStep, min(self._maxStep,
                                                                                           self._targetLeastStep + self._maxExtraStep)))
                break
            for nNode in range(qSize):  # process intermediates one generation at a time
                logging.info("***************************************************")
                logging.info("             processing a new molecule")
                currNode = q.popleft()
                if currNode.smiles == self._productString:
                    continue
                currMol = ob.OBMol(currNode.mol)

                oxidations = {'bond': set(), 'atom': set()}
                for atom in ob.OBMolAtomIter(currMol):
                    nLonePair = numValenceElectron(atom.GetAtomicNum()) - atomTotalBondOrder(
                        atom) + atom.GetFormalCharge()
                    if nLonePair > 0:
                        oxidations['atom'].add(atom)
                for bond in ob.OBMolBondIter(currMol):
                    oxidations['bond'].add(bond)

                def addMol(oxidized, reduced, tempMat=None):
                    logging.debug('in addMol')
                    logging.debug('oxidized: {}\nreduced: {}'.format(oxidized, reduced))
                    if self._matrixForm:
                        # logging.debug('\n'+str(tempMat))
                        newMol = matToMol(tempMat)
                    newMolSmiles = getCanonicalSmiles(newMol)
                    logging.info("newSmiles = " + newMolSmiles)
                    if newMolSmiles == self._productString:
                        logging.info("target found!!!")
                        self._targetLeastStep = nStep
                        self._targetFound = True
                    if newMolSmiles not in self._reactionMap:
                        logging.info("new molecule found! Adding it to the map")
                        if self._doCalculation and self._preEnergyScreen:
                            absoluteEnergy = self.computeQMEnergy(newMol, "gaussian", self._gaussianKeywords,
                                                                  self._fragmentEnergyMap)
                            logging.debug("absoluteEnergy is %f kcal/mol" % (absoluteEnergy))
                            logging.debug("energy base line is " + str(self._energyBaseLine))
                            energy = absoluteEnergy - self._energyBaseLine
                            logging.info("relative energy is %f kcal/mol" % (energy))
                            self._energyMap[newMolSmiles] = energy

                            logging.info("Screening energy")
                            if energy - currNode.energy < self._intermediateThresh:
                                logging.info("low energy intermediate found, adding it to the map...")
                                newNode = ReactionGraphNode(mol=newMol, depth=nStep)
                                newNode.energy = energy
                                self._reactionMap[newMolSmiles] = newNode
                                if newMolSmiles not in currNode.neighbors:
                                    logging.info('adding the edge')
                                    currNode.neighbors[newMolSmiles] = ReactionGraphEdge(currNode, newNode, oxidized,
                                                                                         reduced)
                                    q.append(newNode)
                            else:
                                logging.info("energy too high, discarded")
                        else:
                            newNode = ReactionGraphNode(mol=newMol, depth=nStep)
                            self._reactionMap[newMolSmiles] = newNode
                            if newMolSmiles not in currNode.neighbors:
                                logging.info('adding the edge')
                                currNode.neighbors[newMolSmiles] = ReactionGraphEdge(currNode, newNode, oxidized,
                                                                                     reduced)
                                q.append(newNode)
                    else:
                        logging.info("This molecule has been processed")
                        if currNode.smiles != newMolSmiles:
                            # self._reactionMap[newMolSmiles].depths.append(nStep)
                            if newMolSmiles not in currNode.neighbors:
                                logging.debug("adding {} - {}".format(currNode.smiles, newMolSmiles))
                                logging.debug(
                                    "Although this molecule has been added to reactionMap, "
                                    "it reveals a new route. Adding only the edge...")
                                currNode.neighbors[newMolSmiles] = ReactionGraphEdge(currNode,
                                                                                     self._reactionMap[newMolSmiles],
                                                                                     oxidized, reduced)
                    logging.debug("finish adding this molecule, no matter added or not")
                    # ====================the end of addMol====================

                if self._matrixForm:
                    if self._filterFc:
                        molMat = molToMat(currMol)
                        logging.debug('\n' + str(molMat))
                        eSources = set()
                        for i in self.activeList:
                            if molMat[i][i] > 0:
                                eSources.add((i,))
                            for j in self.activeList:
                                if j < i and molMat[i][j] > 0:
                                    eSources.add((i, j))
                        logging.debug('eSources = {}'.format(eSources))

                        def countChanges(atoms, redox):
                            # keep track of changes to A-BE by the redox operation
                            # redox = -1 if oxidation else 1
                            if len(atoms) is 1:
                                i = atoms[0]
                                changeTable[(i, i)] += 2 * redox
                                fcChange[i] -= 2 * redox
                            else:
                                i, j = atoms
                                changeTable[(i, j)] += 1 * redox
                                changeTable[(j, i)] += 1 * redox
                                fcChange[i] -= 1 * redox
                                fcChange[j] -= 1 * redox
                                tboChange[i] += 1 * redox
                                tboChange[j] += 1 * redox

                        logging.debug('one pair')
                        for eSource1 in eSources:
                            canReduce = set()
                            eTargets = set()
                            for i in self.activeList:
                                if molMat[self.nAtom + 2][i] > self._minFc[molMat[0][i]]:
                                    canReduce.add(i)
                            for atom in eSource1:
                                canReduce.add(atom)
                            canReduce = list(canReduce)
                            for i in range(len(canReduce)):
                                eTargets.add((canReduce[i],))
                                for j in range(i):
                                    eTargets.add((canReduce[i], canReduce[j]))
                            logging.debug('eSource1: {}'.format(eSource1))
                            logging.debug('eTargets: {}'.format(eTargets))
                            for eTarget1 in eTargets:
                                if set(eTarget1) == set(eSource1):
                                    continue
                                changeTable = defaultdict(int)
                                tboChange = defaultdict(int)  # total bond order change
                                fcChange = defaultdict(int)  # formal charge change
                                countChanges(eSource1, -1)
                                countChanges(eTarget1, 1)
                                if self.checkChangeTable(molMat, changeTable, tboChange, fcChange):
                                    tempMat = np.array(molMat)
                                    self.applyChanges(tempMat, changeTable, tboChange, fcChange)
                                    if self.checkStructure(tempMat):
                                        # logging.debug('\n'+str(tempMat))
                                        addMol([eSource1], [eTarget1], tempMat)
                                    else:
                                        logging.debug('this structure did not pass self.checkStructure, not adding it')
                            logging.debug('finishing this eTargets')

                        logging.debug('two pairs')
                        for eSource1, eSource2 in itertools.combinations(eSources, 2):
                            eTargets = set()
                            canReduce = set()
                            for i in self.activeList:
                                if molMat[self.nAtom + 2][i] > self._minFc[molMat[0][i]]:
                                    canReduce.add(i)
                            for atom in eSource1:
                                canReduce.add(atom)
                            for atom in eSource2:
                                canReduce.add(atom)
                            canReduce = list(canReduce)
                            for i in range(len(canReduce)):
                                eTargets.add((canReduce[i],))
                                for j in range(i):
                                    eTargets.add((canReduce[i], canReduce[j]))
                            logging.debug('eSource1 = {}, eSource2 = {}'.format(eSource1, eSource2))
                            logging.debug('eTargets: {}'.format(eTargets))
                            for eTarget1, eTarget2 in itertools.combinations(eTargets - set(eSource1) - set(eSource2),
                                                                             2):
                                changeTable = defaultdict(int)
                                tboChange = defaultdict(int)  # total bond order change
                                fcChange = defaultdict(int)  # formal charge change
                                countChanges(eSource1, -1)
                                countChanges(eTarget1, 1)
                                countChanges(eSource2, -1)
                                countChanges(eTarget2, 1)
                                if self.checkChangeTable(molMat, changeTable, tboChange, fcChange):
                                    tempMat = np.array(molMat)
                                    self.applyChanges(tempMat, changeTable, tboChange, fcChange)
                                    if self.checkStructure(tempMat):
                                        addMol([eSource1, eSource2], [eTarget1, eTarget2], tempMat)
                                    else:
                                        logging.debug('this structure did not pass self.checkStructure, not adding it')

                            logging.debug('finishing this eTargets')

                        if self._allowedPairs >= 3:
                            logging.debug('three pairs')
                            for eSource1, eSource2, eSource3 in itertools.combinations(eSources, 3):
                                eTargets = set()
                                canReduce = set()
                                for i in self.activeList:
                                    if molMat[self.nAtom + 2][i] > self._minFc[molMat[0][i]]:
                                        canReduce.add(i)
                                selectedESources = set(eSource1) | set(eSource2) | set(eSource3)
                                for atom in selectedESources:
                                    canReduce.add(atom)
                                canReduce = list(canReduce)
                                for i in range(len(canReduce)):
                                    eTargets.add((canReduce[i],))
                                    for j in range(i):
                                        eTargets.add((canReduce[i], canReduce[j]))
                                logging.debug(
                                    'eSource1 = {}, eSource2 = {}, eSource3 = {}'.format(eSource1, eSource2, eSource3))
                                logging.debug('eTargets: {}'.format(eTargets))
                                for eTarget1, eTarget2, eTarget3 in itertools.combinations(eTargets - selectedESources,
                                                                                           3):

                                    changeTable = defaultdict(int)
                                    tboChange = defaultdict(int)  # total bond order change
                                    fcChange = defaultdict(int)  # formal charge change
                                    countChanges(eSource1, -1)
                                    countChanges(eTarget1, 1)
                                    countChanges(eSource2, -1)
                                    countChanges(eTarget2, 1)
                                    countChanges(eSource3, -1)
                                    countChanges(eTarget3, 1)
                                    if self.checkChangeTable(molMat, changeTable, tboChange, fcChange):
                                        tempMat = np.array(molMat)
                                        self.applyChanges(tempMat, changeTable, tboChange, fcChange)
                                        if self.checkStructure(tempMat):
                                            addMol([eSource1, eSource2, eSource3], [eTarget1, eTarget2, eTarget3],
                                                   tempMat)
                                        else:
                                            logging.debug(
                                                'this structure did not pass self.checkStructure, not adding it')



                    else:  # no filter at ox/red level
                        molMat = molToMat(currMol)
                        logging.debug('\n' + str(molMat))
                        eSources, eTargets = set(), set()
                        for i in range(1, self.nAtom + 1):
                            if molMat[i][i] > 0:
                                eSources.add((i,))
                            eTargets.add((i,))
                            for j in range(1, i):
                                if molMat[i][j] > 0:
                                    eSources.add((i, j))
                                eTargets.add((i, j))
                        logging.debug(eTargets)

                        def countChanges(atoms, redox):  # redox = -1 if oxidation else 1
                            if len(atoms) is 1:
                                changeTable[(atoms[0], atoms[0])] += 2 * redox
                                fcChange[atoms[0]] -= 2 * redox
                            else:
                                changeTable[atoms[0], atoms[1]] += 1 * redox
                                changeTable[atoms[1], atoms[0]] += 1 * redox
                                fcChange[atoms[0]] -= 1 * redox
                                fcChange[atoms[1]] -= 1 * redox
                                tboChange[atoms[0]] += 1 * redox
                                tboChange[atoms[1]] += 1 * redox

                        for eSource1 in eSources:
                            for eTarget1 in eTargets:
                                changeTable = defaultdict(int)
                                tboChange = defaultdict(int)  # total bond order change
                                fcChange = defaultdict(int)  # formal charge change
                                countChanges(eSource1, -1)
                                countChanges(eTarget1, 1)
                                if self.checkChangeTable(molMat, changeTable, tboChange, fcChange):
                                    tempMat = np.array(molMat)
                                    self.applyChanges(tempMat, changeTable, tboChange, fcChange)
                                    addMol([eSource1], [eTarget1], tempMat)

                        for eSource1 in eSources:
                            for eSource2 in eSources:
                                for eTarget1 in eTargets:
                                    for eTarget2 in eTargets:
                                        changeTable = defaultdict(int)
                                        tboChange = defaultdict(int)  # total bond order change
                                        fcChange = defaultdict(int)  # formal charge change
                                        countChanges(eSource1, -1)
                                        countChanges(eTarget1, 1)
                                        countChanges(eSource2, -1)
                                        countChanges(eTarget2, 1)
                                        if self.checkChangeTable(molMat, changeTable, tboChange, fcChange):
                                            tempMat = np.array(molMat)
                                            self.applyChanges(tempMat, changeTable, tboChange, fcChange)
                                            addMol([eSource1, eSource2], [eTarget1, eTarget2], tempMat)

                else: # not using matrixform
                    eSources = set()
                    for atom in oxidations['atom']:
                        eSources.add(atom.GetIdx())
                    for bond in oxidations['bond']:
                        eSources.add((bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()))
                    eTargets = set()
                    for i in range(1, self.nAtom + 1):
                        eTargets.add(i)
                        for j in range(i + 1, self.nAtom + 1):
                            eTargets.add((i, j))

                    for eSource1 in eSources:
                        for eTarget1 in eTargets:
                            newMol = ob.OBMol(currMol)
                            self.oxidize(newMol, eSource1)
                            self.reduce(newMol, eTarget1)
                            if self.checkLuisRule(eSource1, eTarget1, mol=newMol):
                                addMol([eSource1], [eTarget1])

                    for eSource1 in eSources:
                        for eSource2 in eSources:
                            for eTarget1 in eTargets:
                                for eTarget2 in eTargets:
                                    newMol = ob.OBMol(currMol)
                                    self.oxidize(newMol, eSource1)
                                    self.oxidize(newMol, eSource2)
                                    self.reduce(newMol, eTarget1)
                                    self.reduce(newMol, eTarget2)
                                    if self.checkLuisRule(eSource1, eSource2, eTarget1, eTarget2, mol=newMol):
                                        addMol([eSource1, eSource2], [eTarget1, eTarget2])

        if not self._noProduct:
            logging.info("targetSmiles = " + self._productString)
        else:
            logging.info('no target provided')
        logging.info("targetLeastStep = {}".format(self._targetLeastStep))
        logging.info("===============End of the isomer search===============")
        if self._productString in self._reactionMap:
            return head, self._reactionMap[self._productString]
        else:
            logging.info("target not found")
            return head, None

    def printTextReactionMap(self, head):
        q = deque()
        q.append(ReactionGraphEdge(None, head, [], []))
        visited = set()
        while len(q) > 0:
            qSize = len(q)
            print("\n------------------------")
            for nLevel in range(qSize):
                currEdge = q.popleft()
                # currNode, brokenBonds, createdBonds = q.popleft()
                print(currEdge.node.smiles, 'b ', currEdge.eSources, 'c ', currEdge.eTargets),
                if currEdge.node.smiles not in visited:
                    visited.add(currEdge.node.smiles)
                    for molSmiles, nextEdge in currEdge.node.neighbors.items():
                        q.append(nextEdge)
        print

    def printGraphicReactionMap(self, head):
        q = deque()
        q.append(ReactionGraphEdge(None, head, [], []))
        visited = set()
        if not os.path.isdir("dot"):
            os.system("mkdir dot")
        if not os.path.isdir(self.picDir):
            os.system("mkdir {}".format(self.picDir))
        pwd = os.path.dirname(os.path.realpath(__file__))
        with open("dot" + os.sep + "dot.gv", "w") as dotFile:
            dotFile.write('digraph G  {\nconcentrate = true\nimagepath="' + pwd + os.sep + self.picDir + '"\n')
            edges = []
            nNodes = 0
            while len(q) > 0:
                qSize = len(q)
                nNodes += qSize
                for nLevel in range(qSize):
                    currEdge = q.popleft()
                    if currEdge.node.smiles not in visited:
                        visited.add(currEdge.node.smiles)
                        fileString = smilesToFilename(currEdge.node.smiles)
                        formatString = self.molPicFormat
                        if formatString == 'png':
                            molToPngFile(strToMol('smi', currEdge.node.smiles),
                                         self.picDir + os.sep + fileString + '.' + formatString)
                        else:
                            with open(self.picDir + os.sep + fileString + '.' + formatString,
                                      'w') as picFile:
                                picFile.write(printMol(strToMol('smi', currEdge.node.smiles), formatString))
                        if self._doCalculation:
                            dotFile.write('"{}" [image="{}.{}", label="{} kcal/mol", shape=none, labelloc=b]'
                                          .format(currEdge.node.smiles, fileString, formatString,
                                                  str(currEdge.node.energy)))
                        else:
                            dotFile.write('"{}" [image="{}.{}", label="", shape=none, labelloc=b]'
                                          .format(currEdge.node.smiles, fileString, formatString))
                        for molSmiles, nextEdge in currEdge.node.neighbors.items():
                            if self._pathOnly and nextEdge.onPath or not self._pathOnly:
                                q.append(nextEdge)
                                edges.append((currEdge.node.smiles, nextEdge.node.smiles))
                                if self._doTsSearch:
                                    dotFile.write('   "{}" -> "{}" [ label="{:<8}" ];\n'.format(currEdge.node.smiles,
                                                                                                nextEdge.node.smiles,
                                                                                                str(nextEdge.tsEnergy)))
                                else:
                                    dotFile.write(
                                        '   "{}" -> "{}";\n'.format(currEdge.node.smiles, nextEdge.node.smiles))
            dotFile.write("}\n")
            dotFile.write('//nNodes = {}\n'.format(len(visited)))
            dotFile.write('//nEdges = {}\n'.format(len(edges)))
        return edges

    def findDfsPath(self, head, end, paths, targetLeastStep, path=None):
        if path is None:
            path = [head]
        else:
            path.append(head)
        if len(path) > targetLeastStep + self._maxExtraStep:
            return
        if head == end:
            paths.append(path)
            return
        for molSmiles, edge in head.neighbors.items():
            if edge.node not in path:
                self.findDfsPath(edge.node, end, paths, targetLeastStep, path=list(path))

    def labelPathItems(self, paths, head):
        head.onPath = True
        for path in paths:
            breakPath = False
            for i, node in enumerate(path):
                if breakPath:
                    break
                if i + 1 < len(path):
                    if self._doCalculation and self._energyScreen and not self._preEnergyScreen:
                        if path[i + 1].smiles not in self._energyMap:
                            absoluteEnergy = self.computeQMEnergy(path[i + 1].mol, "gaussian", self._gaussianKeywords,
                                                                  self._fragmentEnergyMap)
                            logging.debug("absoluteEnergy is %f kcal/mol" % (absoluteEnergy))
                            logging.debug("energy base line is " + str(self._energyBaseLine))
                            energy = absoluteEnergy - self._energyBaseLine
                            logging.info("relative energy is %f kcal/mol" % (energy))
                            logging.info("Screening energy")
                            self._energyMap[path[i + 1].smiles] = energy
                        else:
                            logging.debug("energy already calculated")
                            energy = self._energyMap[path[i + 1].smiles]
                        if energy - path[i].energy < self._intermediateThresh:
                            logging.info("low energy intermediate found, marking it as onPath")
                            path[i + 1].energy = energy
                            node.neighbors[path[i + 1].smiles].onPath = True
                            path[i + 1].onPath = True
                        else:
                            logging.info("energy too high, discarded")
                            breakPath = True
                    else:
                        node.neighbors[path[i + 1].smiles].onPath = True
                        path[i + 1].onPath = True

    # def printGraphicPathMap(self, paths):
    #     if not os.path.isdir("dot"):
    #         os.system("mkdir dot")
    #     if not os.path.isdir("static" + os.sep + "pics"):
    #         os.system("mkdir static" + os.sep + "pics")
    #     dotFile = open("dot" + os.sep + "paths.gv", 'w')
    #     dotFile.write("digraph paths {")
    #     visitedNode = set()
    #     visitedEdge = set()
    #     for path in paths:
    #         for i, node in enumerate(path):
    #             if node not in visitedNode:
    #                 visitedNode.add(node)
    #                 if self._doCalculation:
    #                     node.energy = self.computeQMEnergy(node.mol, "gaussian", self._gaussianKeywords,
    #                                                        self._fragmentEnergyMap)
    #                 dotFile.write(
    #                     "    \"" + node.smiles + "\" [image = \".." + os.sep + "static" + os.sep + "pics" + os.sep
    #                     + smilesToFilename(node.smiles) + ".svg\", label = \""
    #                     + str(node.energy) + " kcal/mol\", shape = none, labelloc = b]\n")
    #             if i < len(path) - 1:
    #                 if (node, path[i + 1]) not in visitedEdge:
    #                     visitedEdge.add((node, path[i + 1]))
    #                     dotFile.write("    \"" + node.smiles + "\" -> \"" + path[i + 1].smiles + "\";\n")
    #     dotFile.write("}\n")

    # def getTsEstim(self, node, edge):
    #     mol1 = pybel.readstring('sdf', pybel.Molecule(node.mol).write('sdf'))
    #     mol1.make3D('uff')
    #     for bondData in edge.eTargets:
    #         self.createNewBond(mol1.OBMol, mol1.atoms[bondData[0]-1].OBAtom, mol1.atoms[bondData[1]-1].OBAtom, bondData[2], bondData[3])
    #     mol1.localopt('uff')
    #     mol2 = pybel.readstring('sdf', mol1.write('sdf'))
    #     for bondData in edge.eTargets:
    #         self.breakBond(mol1.OBMol, mol1.atoms[bondData[0]-1].OBAtom, mol1.atoms[bondData[1]-1].OBAtom, bondData[2], bondData[3])
    #     for bondData in edge.eSources:
    #         self.breakBond(mol2.OBMol, mol2.atoms[bondData[0]-1].OBAtom, mol2.atoms[bondData[1]-1].OBAtom, bondData[2], bondData[3])
    #     try:
    #         return SeamTsSearch(mol1, mol2, 'uff')
    #     except TsEstimConvergeError:
    #         print("TS estimate convergence failure")
    #         logging.error("TS estimate convergence failure (SeamTsSearch fails)")
    #         return None
    #
    # def findTsOnPath(self, head):
    #     preQ = [edge for edge in head.neighbors.values() if edge.onPath]
    #     q = deque(preQ)
    #     visitedEdge = set()
    #     if not os.path.isdir('gaussian'):
    #         os.system('mkdir gaussian')
    #     if os.path.isdir('gaussian/ts'):
    #         os.system('rm -f gaussian/ts/*')
    #     else:
    #         os.system('mkdir gaussian/ts')
    #     while q:
    #         currEdge = q.popleft()
    #         visitedEdge.add((currEdge.fromNode.smiles, currEdge.node.smiles))
    #         print('\n========finding TS=======')
    #         print(currEdge.fromNode.smiles, '->', currEdge.node.smiles)
    #         print(visitedEdge)
    #         if (currEdge.node.smiles, currEdge.fromNode.smiles) in visitedEdge:
    #             print('reversed TS is calculated before')
    #             print(currEdge.node.neighbors)
    #             try:
    #                 reverseEdge = currEdge.node.neighbors[currEdge.fromNode.smiles]
    #             except KeyError:
    #                 import pdb; pdb.set_trace()
    #             currEdge.ts = reverseEdge.ts
    #             currEdge.tsEnergy = reverseEdge.tsEnergy
    #         else:
    #             print('calculating TS')
    #             if len(currEdge.eSources) == 0:
    #                 print('pure bond forming reaction, energy goes downhill only, no TS')
    #                 return
    #             if len(currEdge.eTargets) == 0:
    #                 print('pure bond breaking reaction, energy goes uphill only, no TS')
    #                 return
    #
    #             mol = currEdge.fromNode.mol
    #             currTs = self.getTsEstim(currEdge.fromNode, currEdge)
    #             if currTs is not None:
    #                 print('TS esitimate:')
    #                 print(currTs.write('mol'))
    #                 currEdge.ts = currTs
    #                 filename = smilesToFilename(currEdge.fromNode.smiles) + '-' + smilesToFilename(currEdge.node.smiles)
    #                 currTs.title = "ReactionRoute.findTsOnPath"
    #                 opt = {'k': self._gaussianTsKeywords}
    #                 currTs.write('gjf', 'gaussian/ts/'+filename+'.com', overwrite=True, opt=opt)
    #                 currTs.title = ''
    #                 gaussianCall = ''
    #                 for c in filename:
    #                     if c == '(' or c == ')' or c == '$':
    #                         gaussianCall += '\\'
    #                     gaussianCall += c
    #                 print("gdv gaussian/ts/"+gaussianCall+".com")
    #                 logging.info("gdv gaussian/ts/"+gaussianCall+".com")
    #                 success = os.system('gdv gaussian/ts/'+gaussianCall+'.com')
    #                 try:
    #                     absoluteTsEnergy = self.getGaussianEnergy('gaussian/ts/'+filename+'.log')
    #                     currEdge.tsEnergy = absoluteTsEnergy - self._energyBaseLine
    #                     print('TS successfully calculated. The energy is {}'.format(currEdge.tsEnergy))
    #                 except EnergyReadingError:
    #                     currEdge.tsEnergy = 'gauTS E'
    #             else:
    #                 print('TS calculation failed')
    #                 currEdge.ts = None
    #                 currEdge.tsEnergy = 'tsEstim'
    #         for molSmiles, nextEdge in currEdge.node.neighbors.items():
    #             if nextEdge.onPath and (currEdge.node.smiles, molSmiles) not in visitedEdge:
    #                 print('adding {} -> {} to the queue'.format(currEdge.node.smiles, molSmiles))
    #                 q.append(nextEdge)


if __name__ == "__main__":
    logging.basicConfig(filename="log", level=logging.DEBUG, filemode='w')
    rr = ReactionRoute()
    flags = {}
    inputName = None

    parser = argparse.ArgumentParser(description='Search for a possible reaction mechanism')
    parser.add_argument('-j', help='provide a json input file', metavar='filename')
    parser.add_argument('-r', help='reactant SMILES', metavar='SMILES')
    parser.add_argument('-p', help='product SMILES', metavar='SMILES')
    parser.add_argument('-e', action='store_true',
                        help='enables calculations of intermediate energies and energy screen '
                             '(default threshold is 200kcal/mol. Use --intermThresh [energy] to change it)')
    parser.add_argument('--intermThresh', type=float,
                        help='energy in kcal/mol. Intermediates that have higher energy than this threshold '
                             'will be excluded from the reaction network. Default is 200 kcal/mol')
    parser.add_argument('-q', action='store_true',
                        help='use gsub to submit calculation jobs to the queue instead of running '
                             'on the local machine. Only works on pinconning or MERCED. ')
    parser.add_argument('-n', action='store_true',
                        help='exploration mode. No product is specified. The program will start the search '
                             'from the reactant and stop at the maximum number of steps allowed. '
                             'All relevant routes will be in the network')
    parser.add_argument('-a', action='store_true',
                        help='selecting active atom mode. No search will be run. A gaussian input file '
                             'will be generated instead to help you select active atoms. ')
    args = parser.parse_args()

    if not args.j and not args.r:
        parser.error('at least a reactant should be provided, either in an input file or through -r option')
    if args.j:
        inputName = args.j[:-5]
        with open(inputName + '.json') as f:
            rr.inputJson(f.read())
    if args.r:
        rr._reactantString = args.r
    if args.p:
        rr._productString = args.p
    if args.e:
        rr._doCalculation = True
        rr._energyScreen = True
    if args.q:
        rr._gsub = True
    if args.n:
        rr._noProduct = True
    if args.a:
        pymol = pybel.readstring('smi', rr._reactantString)
        mol = pymol.OBMol
        builder = ob.OBBuilder()
        builder.Build(mol)
        separateFragments(pymol)
        pymol.title = 'for select active atoms'
        pymol.addh()
        pymol.localopt()
        with open('activeatoms.com', 'w') as f:
            f.write(pymol.write('gjf'))
        exit()

    # cProfile.run('head, target= rr.isomerSearch()')
    head, target = rr.isomerSearch()
    # rr.printTextReactionMap(head)
    if target is not None and not rr._noProduct:
        paths = []
        rr.findDfsPath(head, target, paths, rr._targetLeastStep)
        rr.labelPathItems(paths, head)
    else:
        rr._pathOnly = False

    # if rr._doTsSearch:
    #     rr.findTsOnPath(head)
    edges = rr.printGraphicReactionMap(head)
    print(edges)
    if rr.jobName == 'tmpName':
        rr.jobName = os.path.basename(inputName)
    with open('dot' + os.sep + '{}.gv'.format(rr.jobName), 'w') as dotF:
        with open('{}.json'.format(inputName)) as inputF:
            for line in inputF:
                dotF.write('//{}'.format(line))
        with open('dot' + os.sep + 'dot.gv') as dotF_origin:
            dotF.write(dotF_origin.read())
    print("dot -Tsvg dot" + os.sep + "dot.gv -o dot" + os.sep + "{}.svg".format(rr.jobName))
    os.system("dot -Tsvg dot{}dot.gv -o {}{}{}.svg".format(os.sep, rr.picDir, os.sep, rr.jobName))
