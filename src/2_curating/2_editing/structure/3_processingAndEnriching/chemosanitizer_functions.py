# generic modules
import concurrent.futures
import time
# import multiprocessing
import pandas as pd
import sys
import gzip

# RDkit specific modules
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import SanitizeMol
from rdkit.Chem.MolStandardize import rdMolStandardize


# MolVS specific modules
import molvs
from molvs import Standardizer
from molvs import Validator
# from molvs.fragment import LargestFragmentChooser
# from molvs.fragment import FragmentRemover
# or eventually load the local (modified) fragment.py
import fragment
from fragment import LargestFragmentChooser
from fragment import FragmentRemover
from molvs.charge import Uncharger


def printer(inchi):
    print(inchi)


def MolFromInchi_fun(inchi):
    m = Chem.MolFromInchi(inchi)
    if m:
        return m
    return None


def MolFromSmiles_fun(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m:
        return m
    return None


def MolToSmiles_fun(romol):
    m = Chem.MolToSmiles(romol)
    if m:
        return m
    return None


def MolToInchi_fun(romol):
    m = Chem.MolToInchi(romol)
    if m:
        return m
    return None

        
def MolToIK_fun(romol):
    m = Chem.MolToInchiKey(romol)
    if m:
        return m
    return None


def MolToFlatMol_fun(romol):
    Chem.RemoveStereochemistry(romol) # See MOLVS examples (https://programtalk.com/python-examples/rdkit.Chem.RemoveStereochemistry/)
    return romol


def MolToMF_fun(romol):
    m = rdMolDescriptors.CalcMolFormula(romol)
    if m:
        return m
    return None


def MolToEmass_fun(romol):
    m = Descriptors.ExactMolWt(romol)
    if m:
        return m
    return None


def MolToLogP_fun(romol):
    m = Chem.Crippen.MolLogP(romol)
    if m:
        return m
    return None


# defining the validator log output format
fmt = '%(levelname)s - %(validation)s - %(message)s'


# save the Standardizer and LargestFragmentChooser classes as variables


def validator_fun(romol):
    m = Validator(log_format=fmt).validate(romol)
    if m:
        return m
    return None


def standardizor_fun(romol):
    m = Standardizer().standardize(romol)
    if m:
        return m
    return None


def canonicalizor_fun(romol):
    m = rdMolStandardize.TautomerEnumerator().Canonicalize(romol)
    if m:
        return m
    return None


def fragremover_fun(romol):
    m = FragmentRemover().remove(romol)
    if m:
        return m
    return None


def fragchooser_fun(romol):
    m = LargestFragmentChooser().choose(romol)
    if m:
        return m
    return None


def uncharger_fun(romol):
    m = Uncharger().uncharge(romol)
    if m:
        return m
    return None


def sik_fun(ik):
    if ik:
        m = str.partition(ik,"-")[0]
        if m:
            return m
        return None    
    return None


#  Courtesy of Richy Leroy (via Jean-Marc Nuzillard)

# iminol
smarts1 = "[NH0]=C([OH])([!O])" # secondary iminol
target1 = Chem.MolFromSmarts(smarts1)
smarts11 = "[NH1]=C([OH])([!O])" # primary iminol
target11 = Chem.MolFromSmarts(smarts11)
smarts12 = "[A&!H][NH0]=C([OH])O" # secondary carbamate
target12 = Chem.MolFromSmarts(smarts12)
smarts13 = "[NH]=C([OH])O" # primary carbamate
target13 = Chem.MolFromSmarts(smarts13)
smarts14 = "[CH]([OH])=N" # N-formyl
target14 = Chem.MolFromSmarts(smarts14)

# enol
smarts2 = "[C&!c]([!CX3&!O])([!CX3&!O])=[C&!c]([OH])([!N])"
target2 = Chem.MolFromSmarts(smarts2)

# enethiol
smarts3 = "[C&!c]([!CX3&!SHX2])([!CX3&!SHX2])=[C&!c]([SHX2])"
target3 = Chem.MolFromSmarts(smarts3)
smarts31 = "[N&!n]([!CX3&!SHX2])([!CX3&!SHX2])=[C&!c]([SHX2])"
target31 = Chem.MolFromSmarts(smarts31)


def tautomerizor_fun(m):

    ps = m

    # N-formyl
    if ps.HasSubstructMatch(target14):
        nb14 = len(ps.GetSubstructMatches(target14))
        #print(nb14, "N-Formyl")
        rxn11 = AllChem.ReactionFromSmarts('[CH1:1]([OH:2])=[N:3]>>[CH1:1](=[OH0D1:2])[NH:3]')
        ps = rxn11.RunReactants((ps,))
        ps = ps[0][0]
        if nb14 != 1 :
            for i in range(1,nb14) :
                ps = rxn11.RunReactants((ps,))
                ps = ps[0][0]
        Chem.SanitizeMol(ps)
        Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
    # Fails
    # Primary carbamate
    # if ps.HasSubstructMatch(target13):
    #     nb13 = len(ps.GetSubstructMatches(target13))
    #     #print(nb13, "primary carbamate")
    #     rxn11 = AllChem.ReactionFromSmarts('[C:1][O:2][CD3:3]([OH1:4])=[NH:5]>>[C:1][O:2][CD3:3](=[OH0D1:4])[NH2:5]')
    #     ps = rxn11.RunReactants((ps,))
    #     ps = ps[0][0]
    #     if nb13 != 1 :
    #         for i in range(1,nb13) :
    #             ps = rxn11.RunReactants((ps,))
    #             ps = ps[0][0]
    #     Chem.SanitizeMol(ps)
    #     Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
    # Secondary carbamate
    if ps.HasSubstructMatch(target12):
        nb12 = len(ps.GetSubstructMatches(target12))
        #print(nb12, "secondary carbamate")
        rxn11 = AllChem.ReactionFromSmarts('[O:4][CD3:1]([OH:2])=[NH0:3]>>[O:4][CD3:1](=[OH0D1:2])[NH:3]')
        ps = rxn11.RunReactants((ps,))
        ps = ps[0][0]
        if nb12 != 1 :
            for i in range(1,nb12) :
                ps = rxn11.RunReactants((ps,))
                ps = ps[0][0]
        Chem.SanitizeMol(ps)
        Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
    # Primary iminol
    if ps.HasSubstructMatch(target11):
        nb11 = len(ps.GetSubstructMatches(target11))
        #print(nb11, "primary iminol")
        rxn11 = AllChem.ReactionFromSmarts('[CD3:1]([OH:2])=[NH:3]>>[CD3:1](=[OH0D1:2])[NH2:3]')
        ps = rxn11.RunReactants((ps,))
        ps = ps[0][0]
        if nb11 != 1 :
            for i in range(1,nb11) :
                ps = rxn11.RunReactants((ps,))
                ps = ps[0][0]
        Chem.SanitizeMol(ps)
        Chem.AssignStereochemistry(ps,force=True,cleanIt=True)          
        if ps.HasSubstructMatch(target11)==True :
            ps = rxn11.RunReactants((ps,))
            ps = ps[0][0]
            Chem.SanitizeMol(ps)
            Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
    # Secondary iminol
    if ps.HasSubstructMatch(target1)==True :
        nb1 = len(ps.GetSubstructMatches(target1))
        #print(nb1, "secondary iminol")
        rxn1 = AllChem.ReactionFromSmarts('[C:1]([OH:2])=[NH0:3]>>[C:1](=[OH0:2])[NH:3]')
        ps = rxn1.RunReactants((ps,))
        ps = ps[0][0]
        if nb1 != 1 :
            for i in range(1,nb1) :
                ps = rxn1.RunReactants((ps,))
                ps = ps[0][0]
        Chem.SanitizeMol(ps)
        Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
        if ps.HasSubstructMatch(target1)==True :
            ps = rxn1.RunReactants((ps,))
            ps = ps[0][0]
            Chem.SanitizeMol(ps)
            Chem.AssignStereochemistry(ps,force=True,cleanIt=True)      
    # enol
    if ps.HasSubstructMatch(target2):
        nb2 = len(ps.GetSubstructMatches(target2))
        #print(nb2, "enol")
        rxn2 = AllChem.ReactionFromSmarts('[!c&C:1]([OH:2])=[!c&C:3]>>[C:1](=[OH0:2])[CH:3]')
        ps = rxn2.RunReactants((ps,))
        ps = ps[0][0]
        if nb2 != 1 :
            for i in range(1,nb2) :
                ps = rxn2.RunReactants((ps,))
                ps = ps[0][0]
        Chem.SanitizeMol(ps)
        Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
    # enethiol  
    if ps.HasSubstructMatch(target3):
        nb3 = len(ps.GetSubstructMatches(target3))
        #print(nb3, "enethiol")
        rxn3 = AllChem.ReactionFromSmarts('[!c&C:1]([SH:2])=[!c&C:3]>>[C:1](=[SH0:2])[CH:3]')
        ps = rxn3.RunReactants((ps,))
        ps = ps[0][0]
        if nb3 != 1 :
            for i in range(1,nb3) :
                ps = rxn3.RunReactants((ps,))
                ps = ps[0][0]
        Chem.SanitizeMol(ps)
        Chem.AssignStereochemistry(ps,force=True,cleanIt=True)
        Chem.AddHs(ps)

    return ps


def long_cleaning_function(myslice, smiles_column_header):
    myslice['ROMol'] = myslice[smiles_column_header].apply(MolFromSmiles_fun)
    myslice = myslice[~myslice['ROMol'].isnull()]
    myslice['validatorLog'] = myslice['ROMol'].apply(validator_fun)
    myslice['ROMolSanitizedLargestFragmentUncharged'] = myslice['ROMol'].apply(standardizor_fun).apply(fragremover_fun).apply(uncharger_fun).apply(tautomerizor_fun)
    myslice['smilesSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToSmiles_fun)
    myslice['inchiSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToInchi_fun)
    myslice['inchikeySanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToIK_fun)
    myslice['flatROMol'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToFlatMol_fun)
    myslice['smilesSanitizedFlat'] = myslice['flatROMol'].apply(MolToSmiles_fun)
    myslice['inchiSanitizedFlat'] = myslice['flatROMol'].apply(MolToInchi_fun)
    myslice['shortikSanitized'] = myslice['inchikeySanitized'].apply(sik_fun)
    myslice['formulaSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToMF_fun)
    myslice['exactmassSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToEmass_fun)
    myslice['xlogpSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToLogP_fun)
    myslice = myslice.drop(myslice.loc[:, ['ROMol', 'flatROMol', 'ROMolSanitizedLargestFragmentUncharged']], axis=1)
    return myslice


class CleaningFunc:
    def __init__(self, smiles_column_header):
        self.header = smiles_column_header

    def f(self, myslice):
        return long_cleaning_function(myslice, self.header)