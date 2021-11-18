# generic modules
import time
import multiprocessing
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


def MolToInchi_fun_safe(row):
    if '[O]' not in row['smilesSanitized']:
        m = Chem.MolToInchi(row['ROMolSanitizedLargestFragmentUncharged'])
        if m:
            return m
        return None
    else:
        print('Sayonara Robocop !')
        return None


def MolToInchi_fun_safe_flat(row):
    if '[O]' not in row['smilesSanitizedFlat']:
        m = Chem.MolToInchi(row['flatROMol'])
        if m:
            return m
        return None
    else:
        print('Sayonara Robocop !')
        return None
        
def MolToIK_fun(romol):
    m = Chem.MolToInchiKey(romol)
    if m:
        return m
    return None


def MolToIK_fun_safe(row):
    if '[O]' not in row['smilesSanitized']:
        m = Chem.MolToInchiKey(row['ROMolSanitizedLargestFragmentUncharged'])
        if m:
            return m
        return None
    else:
        print('Sayonara Babeee !')
        return None


# It looks like the  Chem.RemoveStereochemistry should be defined this way and not assigned to a variable. See MOLVS examples (https://programtalk.com/python-examples/rdkit.Chem.RemoveStereochemistry/)

def MolToFlatMol_fun(romol):
    Chem.RemoveStereochemistry(romol)
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


def long_cleaning_function(myslice, smiles_column_header):
    myslice['ROMol'] = myslice[smiles_column_header].apply(MolFromSmiles_fun)
    myslice = myslice[~myslice['ROMol'].isnull()]
    myslice['validatorLog'] = myslice['ROMol'].apply(validator_fun)
    myslice['ROMolSanitized'] = myslice['ROMol'].apply(standardizor_fun)
    myslice['ROMolSanitizedLargestFragment'] = myslice['ROMolSanitized'].apply(fragremover_fun)
    myslice['ROMolSanitizedLargestFragmentUncharged'] = myslice['ROMolSanitizedLargestFragment'].apply(uncharger_fun)
    myslice['smilesSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToSmiles_fun)
    myslice['inchiSanitized'] = myslice.apply(MolToInchi_fun_safe, axis=1)

    myslice['flatROMol'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToFlatMol_fun)
    myslice['smilesSanitizedFlat'] = myslice['flatROMol'].apply(MolToSmiles_fun)
    myslice['inchiSanitizedFlat'] = myslice.apply(MolToInchi_fun_safe_flat, axis=1)
    myslice['inchikeySanitized'] = myslice.apply(MolToIK_fun_safe, axis=1)
    myslice['shortikSanitized'] = myslice['inchikeySanitized'].str.split("-", n=1, expand=True)[0]
    myslice['formulaSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToMF_fun)
    myslice['exactmassSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToEmass_fun)
    myslice['xlogpSanitized'] = myslice['ROMolSanitizedLargestFragmentUncharged'].apply(MolToLogP_fun)
    return myslice


class CleaningFunc:
    def __init__(self, smiles_column_header):
        self.header = smiles_column_header

    def f(self, myslice):
        return long_cleaning_function(myslice, self.header)