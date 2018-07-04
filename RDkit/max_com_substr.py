# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 11:58:36 2018
Maximum Common Substructures and Atom Pairs
Command line: python3 max_com_substr.py
@author: stokm006
"""
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.AtomPairs import Pairs
from rdkit import DataStructs

def max_common_substructures():
    """ The FindMCS functionality finds a maximum common substructure (MCS) of 
    two or more molecules
    
    """
    
    # Generate molecules
    mol1 = Chem.MolFromSmiles("c1ccccc1O")
    mol2 = Chem.MolFromSmiles("c1ccccc1")
    mol3 = Chem.MolFromSmiles("c1ccccc1N")

    mols = [mol1,mol2,mol3]
    res=rdFMCS.FindMCS(mols)
    
    # number of common atoms
    print (res.numAtoms)

    # number of common bonds    
    print (res.numBonds)
    
    # common smart structure
    print (res.smartsString)
    
    mols = (Chem.MolFromSmiles('c1ccccc1'),Chem.MolFromSmiles('C1CCCC=C1'))
    # compare atom difference
    print (rdFMCS.FindMCS(mols,atomCompare=rdFMCS.AtomCompare.CompareAny).smartsString)    
    
    # compare bond difference    
    print (rdFMCS.FindMCS(mols,bondCompare=rdFMCS.BondCompare.CompareAny).smartsString)
    
def atom_pairs():
    """ Atom pair fingerprints, atom descriptor
    
    """
    
    # Generate molecules
    ms = [Chem.MolFromSmiles('C1CCC1OCC'),Chem.MolFromSmiles('CC(C)OCC'),Chem.MolFromSmiles('CCOCC')]
    pairFps = [Pairs.GetAtomPairFingerprint(x) for x in ms]
    
    # Get the list of bits and their counts for each fingerprint as a dictionary
    d = pairFps[-1].GetNonzeroElements()
    print (d)
    
    # Explanation of the bitscore.
    print (Pairs.ExplainPairScore(558115))
    
    # Dice similarity; The usual metric for similarity between atom-pair fingerprints
    print (DataStructs.DiceSimilarity(pairFps[0],pairFps[1]))
    
    # Atom decriptor without count
    pairFps = [Pairs.GetAtomPairFingerprintAsBitVect(x) for x in ms]
    print (DataStructs.DiceSimilarity(pairFps[0],pairFps[1]))

if __name__ == '__main__':
  #  max_common_substructures()
    atom_pairs()

