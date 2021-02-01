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

# def MolFromInchi_fun(inchi):
#     m = Chem.MolFromInchi(inchi)
#     if m:
#         return m
#     print(inchi)


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


def MolToInchi_fun_safe(smiles, romol):
    # print(smiles) (beware as too much print statement will crash \\
    # interactve python console such as the one in vscode)
    if '[O]' not in smiles:
        m = Chem.MolToInchi(romol)
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


def MolToIK_fun_safe(smiles, romol):
    # print(smiles) (beware as too much print statement will crash \\
    # interactve python console such as the one in vscode)
    if '[O]' not in smiles:
        m = Chem.MolToInchiKey(romol)
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
    print(romol.GetNumAtoms)
    m = Validator(log_format=fmt).validate(romol)
    if m:
        return m
    return None


def standardizor_fun(romol):
    print('standardizer ' + str(romol.GetNumAtoms))
    m = Standardizer().standardize(romol)
    if m:
        return m
    return None


def fragremover_fun(romol):
    print('fragremover ' + str(romol.GetNumAtoms))
    m = FragmentRemover().remove(romol)
    if m:
        return m
    return None


def uncharger_fun(romol):
    print('uncharger ' + str(romol.GetNumAtoms))
    m = Uncharger().uncharge(romol)
    if m:
        return m
    return None
