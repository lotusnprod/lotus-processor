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
fmt = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'


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


def fragremover_fun(romol):
    m = FragmentRemover().remove(romol)
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


def long_cleaning_function(myslice, smiles_column_header):
    myslice['ROMol'] = myslice[smiles_column_header].apply(MolFromSmiles_fun)
    myslice = myslice[~myslice['ROMol'].isnull()]
    myslice['validatorLog'] = myslice['ROMol'].apply(validator_fun)
    myslice['ROMolSanitized'] = myslice['ROMol'].apply(standardizor_fun)
    myslice['ROMolSanitizedLargestFragment'] = myslice['ROMolSanitized'].apply(fragremover_fun)
    myslice['ROMolSanitizedLargestFragmentUncharged'] = myslice['ROMolSanitizedLargestFragment'].apply(uncharger_fun)
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
    return myslice


class CleaningFunc:
    def __init__(self, smiles_column_header):
        self.header = smiles_column_header

    def f(self, myslice):
        return long_cleaning_function(myslice, self.header)