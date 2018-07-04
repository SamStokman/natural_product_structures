# -*- coding: utf-8 -*-
"""
Created on Thu May 24 14:51:39 2018
Try RD kit
Command line: python3 rdkit_try_script.py
@author: stokm006
"""

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
def use_rdkit1():
    m = Chem.MolFromSmiles('C1=CC=CN=C1')
   
    # C1=CC=CN=C1 == c1cccnc1 == n1ccccc1, all smiles, not smarts!!
    print ('\n')
    print ('m (mol from smiles) -->',m)
    print ('\n')
    print ('mol to smiles -->',Chem.MolToSmiles(m))
  
    # the kekule form of smile
    Chem.Kekulize(m)
    print ('\n')    
    print ('mol to smiles and Kekulize -->', '\n',Chem.MolToSmiles(m,kekuleSmiles=True))
    print ('\n')     
    
    # create molfile    
    print ('mol to block (molfile) -->',Chem.MolToMolBlock(m)) 
    m.SetProp("_Name","cyclobutane")
    print ('mol to block (molfile) with name -->', '\n',Chem.MolToMolBlock(m))
    # include 2D coordinates
    AllChem.Compute2DCoords(m)
    print('mol to block (molfile) with 2D coordinates -->', '\n', Chem.MolToMolBlock(m))
    # include 3D coordinates
    # add hydrogens
    m2 = Chem.AddHs(m)  
    AllChem.EmbedMolecule(m2, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    print('mol to block (molfile) with 3D coordinates and hydrogen -->', '\n', Chem.MolToMolBlock(m2))
    # remove hydogens
    m3 = Chem.RemoveHs(m2)
    AllChem.EmbedMolecule(m3, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    print('mol to block (molfile) with 3D coordinates without hydrogen -->', '\n', Chem.MolToMolBlock(m3))

    # looping over atoms and bonds
    for atom in m.GetAtoms():
        print('Atomic nr -->', atom.GetAtomicNum()) # 0.5 * AtomicNum C=12 g/mol, N=14 g/mol

 
 
    print('Symbol -->', m.GetAtomWithIdx(0).GetSymbol())
    print('Explicit Valence -->', m.GetAtomWithIdx(0).GetExplicitValence())
    print('Begin Atom Idex -->', m.GetBondWithIdx(0).GetBeginAtomIdx()) 
    print('End Atom Idex --', m.GetBondWithIdx(0).GetEndAtomIdx())   
    
    # bond type
    print('Bond type -->', m.GetBonds()[0].GetBondType())  
    print('Bondtype -->', m.GetBondWithIdx(0).GetBondType())
    print('Bond type between atoms -->', m.GetBondBetweenAtoms(0,1).GetBondType())
    print ('\n')
    
    # ring information
    m4 = Chem.MolFromSmarts('OC1C2C1CC2')

    
    # m4.GetAtomWithIdx(0).IsInRing()  0 == 0th atom
    print ('Boolean if atom is in ring -->', m4.GetAtomWithIdx(0).IsInRing())
    # m4.GetAtomWithIdx(2).IsInRingSize(4)  2 == 2nd atom, ringsize can also be 3
    print ('Boolean if atom is in ring with ringsize -->', m4.GetAtomWithIdx(2).IsInRingSize(4))
        
    # More detail about the smallest set of smallest rings (SSSR) is available
    ssr = Chem.GetSymmSSSR(m4)
    print('Smallest set of smallest ring (SSR) -->', ssr)
    print('len(SSR) -->', len(ssr))
    print('list(SSR[0]) -->', list(ssr[0]))
    print('list(SSR[1]) -->', list(ssr[1]))
    print('nr of SSSR/ nr of rings -->', Chem.GetSSSR(m4))  
    ri = m4.GetRingInfo()   
    print('nr of rings the atom is involved in -->', ri.NumAtomRings(2))
    print('boolean is atom in ring of size? -->', ri.IsAtomInRingOfSize(1,4)) # 1 = atom 4 = atoms in ring

from rdkit.Chem import Draw       
def draw_2D_structures():
    m = Chem.MolFromSmiles('C1=CC=CN=C1')
    m2 = Chem.AddHs(m)
    m4 = Chem.MolFromSmiles('OC1C2C1CC2')
    m5 = Chem.MolFromSmiles('c1nccc2n1ccc2')
    m_list = m, m2, m4, m5   
  
    # 2D; creates a png file with molecule in it
    Draw.MolToFile(m,'/home/stokm006/thesis/RDkit/images/molecule_m.png')
    Draw.MolToFile(m2,'/home/stokm006/thesis/RDkit/images/molecule_m2.png')
    Draw.MolToFile(m4,'/home/stokm006/thesis/RDkit/images/molecule_m4.png')
    Draw.MolToFile(m5,'/home/stokm006/thesis/RDkit/images/molecule_m5.png')
    # multiple molecules
    multiple_molecules = Draw.MolsToGridImage(m_list,molsPerRow=2,subImgSize=(100,100))
    multiple_molecules.save('/home/stokm006/thesis/images/multiple_molecules.png')
           
   
def use_rdkit2():
    # substructure matching
    m6 = Chem.MolFromSmiles('c1ccccc1O')
    patt = Chem.MolFromSmarts('ccO')
    print('Boolean if substructure matches --> ', m6.HasSubstructMatch(patt))
    print('Atoms in structure which match --> ',m6.GetSubstructMatch(patt))
    print('if multiple structures match--> ',m6.GetSubstructMatches(patt))
    
    # chemical transformations
    # deleting a substructure
    rm = AllChem.DeleteSubstructs(m6, patt)
    print ('The removed atoms -->', Chem.MolToSmiles(rm))
#    Draw.MolToFile(rm,'/home/stokm006/thesis/images/molecule_rm.png')
    # replacing a substructure
    m7 = Chem.MolFromSmiles('CC(=O)N')
    patt2 = Chem.MolFromSmarts('[$(NC(=O))]')
    repl = Chem.MolFromSmiles('OC')
    rms = AllChem.ReplaceSubstructs(m7,patt2,repl)
    print('New SMILE -->', Chem.MolToSmiles(rms[0]))


from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Descriptors
from rdkit.Chem import MACCSkeys



def use_rdkit3():
    # Fingerprinting and Molecular Similarity, default; Tanimoto similarity
    m_list2 = [Chem.MolFromSmiles('CCOC'), Chem.MolFromSmiles('CCO'), \
    Chem.MolFromSmiles('COC')]
    fps = [FingerprintMols.FingerprintMol(x) for x in m_list2]
    print('Fingerprint Similarity -->', DataStructs.FingerprintSimilarity(fps[0],fps[1]))
    print('Fingerprint Similarity -->', DataStructs.FingerprintSimilarity(fps[0],fps[0]))
    print('Fingerprint Similarity -->', DataStructs.FingerprintSimilarity(fps[2],fps[1]))
    # MACCS keys    
    fps = [MACCSkeys.GenMACCSKeys(x) for x in m_list2]
    print('Fingerprint MACCS keys -->', DataStructs.FingerprintSimilarity(fps[0],fps[1]))
    print('Fingerprint MACCS keys -->', DataStructs.FingerprintSimilarity(fps[0],fps[0]))
    print('Fingerprint MACCS keys -->', DataStructs.FingerprintSimilarity(fps[2],fps[1]))
    # Morgan/ circular fingerprints
    m7 = Chem.MolFromSmiles('CCOC')
    fp1 = AllChem.GetMorganFingerprint(m7,2)
    m8 = Chem.MolFromSmiles('CCO')
    fp2 = AllChem.GetMorganFingerprint(m8,2)
    print('Fingerprint Morgan similarity -->', DataStructs.DiceSimilarity(fp1,fp2))
    fp1 = AllChem.GetMorganFingerprintAsBitVect(m7,2,nBits=1024)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(m8,2,nBits=1024)
    print(DataStructs.DiceSimilarity(fp1,fp2))
    
    ffp1 = AllChem.GetMorganFingerprint(m7,2,useFeatures=True)
    ffp2 = AllChem.GetMorganFingerprint(m8,2,useFeatures=True)
    print(DataStructs.DiceSimilarity(ffp1,ffp2))

    #explaining bit (Morgan)
    m9 = Chem.MolFromSmiles('c1cccnc1C')
    info={}
    AllChem.GetMorganFingerprint(m9,2,bitInfo=info)
    print (info)
    print (info[98513984])
    print (info[4048591891])
    
    env = Chem.FindAtomEnvironmentOfRadiusN(m9,2,5)
    amap={}
    submol=Chem.PathToSubmol(m9,env,atomMap=amap)
    print (submol.GetNumAtoms())
    print (amap)
    # bit to smile
    print (Chem.MolToSmiles(submol))


    # Descriptor Calculation; used in papers or coding languages
    m6 = Chem.MolFromSmiles('c1ccccc1O')
    print('Descriptor TPSA -->', Descriptors.TPSA(m6))
    print('Descriptor MolLogP -->', Descriptors.MolLogP(m6))
    
    # Chemical reactions
#    rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2].[C:3]=[*:4][*:5]=[C:6]>>[C:1]1[C:2][C:3][*:4]=[*:5][C:6]1')
    rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]')
    ps = rxn.RunReactants((Chem.MolFromSmiles('CC(=O)O'),Chem.MolFromSmiles('NC=C')))
    print ('Reaction product -->', Chem.MolToSmiles(ps[0][0]))
    
    
if __name__ == "__main__":
  #  use_rdkit1()
  #  draw_2D_structures()
  #   use_rdkit2()
    use_rdkit3()
