# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 11:58:36 2018
Maximum Common Substructures
Command line: python3 max_com_substr.py
@author: stokm006
"""
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.AtomPairs import Pairs
from rdkit import DataStructs

def max_common_substructures():
    mol1 = Chem.MolFromSmiles("c1ccccc1O")
    mol2 = Chem.MolFromSmiles("c1ccccc1")
    mol3 = Chem.MolFromSmiles("c1ccccc1N")

    mols = [mol1,mol2,mol3]
    res=rdFMCS.FindMCS(mols)
    print (res.numAtoms)
    print (res.numBonds)
    print (res.smartsString)
    
    mols = [mol1, mol3]
    
    mols = (Chem.MolFromSmiles('NCC'),Chem.MolFromSmiles('OC=C'))
    rdFMCS.FindMCS(mols).smartsString
    print (rdFMCS.FindMCS(mols, atomCompare=rdFMCS.AtomCompare.CompareAny).smartsString)
    print (rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareAny).smartsString)
    print (rdFMCS.FindMCS(mols,bondCompare=rdFMCS.BondCompare.CompareOrderExact).smartsString)
    print (rdFMCS.FindMCS(mols,bondCompare=rdFMCS.BondCompare.CompareOrder).smartsString)

    mols = (Chem.MolFromSmiles('c1ccccc1'),Chem.MolFromSmiles('C1CCCC=C1'))
    print (rdFMCS.FindMCS(mols,bondCompare=rdFMCS.BondCompare.CompareAny).smartsString)
    print (rdFMCS.FindMCS(mols,bondCompare=rdFMCS.BondCompare.CompareOrderExact).smartsString)
    print (rdFMCS.FindMCS(mols,bondCompare=rdFMCS.BondCompare.CompareOrder).smartsString)

    mols = [Chem.MolFromSmiles("Nc1ccccc1"*5), Chem.MolFromSmiles("Nc1ccccccccc1")]
    print (rdFMCS.FindMCS(mols).numAtoms)
    
def atom_pairs():
    ms = [Chem.MolFromSmiles('C1CCC1OCC'),Chem.MolFromSmiles('CC(C)OCC'),Chem.MolFromSmiles('CCOCC')]
    pairFps = [Pairs.GetAtomPairFingerprint(x) for x in ms]
    d = pairFps[-1].GetNonzeroElements()
    print (d)
    print (Pairs.ExplainPairScore(558115))
    #dice similarity
    print (DataStructs.DiceSimilarity(pairFps[0],pairFps[1]))

    pairFps = [Pairs.GetAtomPairFingerprintAsBitVect(x) for x in ms]
    
    print (pairFps)
    print (DataStructs.DiceSimilarity(pairFps[0],pairFps[1]))

if __name__ == '__main__':
   # max_common_substructures()
    atom_pairs()


