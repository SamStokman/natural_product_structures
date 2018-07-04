# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:57:45 2018
Simple python script that generates similarity maps using RDkit.
@author: stokm006
"""

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import AllChem


def draw_map():
    """ Creates similarity maps using the provided molecules (SMILES).
    
    """
    m1 = Chem.MolFromSmiles('c1ccccc1O')
    m2 = Chem.MolFromSmiles('c1ccccc1N')
    
    # Morgan Fingerprint (with normalization)
    # Can also be used with APFingerprint or TTFingerprint
    fig1, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(m1, m2, SimilarityMaps.GetMorganFingerprint)
    fig1.savefig('/path/to/similaritymap.png',bbox_inches='tight')
   
    # TT Fingerprint (with normalization)
    fig2, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(m1, m2, SimilarityMaps.GetTTFingerprint)
    fig2.savefig('/path/to/similaritymap.png',bbox_inches='tight')

    # Morgan Fingerprint (without normalization)
    weights = SimilarityMaps.GetAtomicWeightsForFingerprint(m1, m2, SimilarityMaps.GetMorganFingerprint)
    fig3 = SimilarityMaps.GetSimilarityMapFromWeights(m2, weights, size=(150, 150))
    fig3.savefig('/path/to/similaritymap.png',bbox_inches='tight')    
    
    # the degree of partial charge by using atomic charge
    AllChem.ComputeGasteigerCharges(m1)
    charges = [float(atom.GetProp('_GasteigerCharge')) for atom in m1.GetAtoms()]
    fig4 = SimilarityMaps.GetSimilarityMapFromWeights(m2,charges, size=(150, 150),scale=10)
    fig4.savefig('/path/to/molcharge_similaritymap.png',bbox_inches='tight')
  
   
if __name__ == "__main__":
    draw_map()
        