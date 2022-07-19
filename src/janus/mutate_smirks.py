#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 12:15:57 2021

@author: akshat
"""
from typing import Dict
import random
import multiprocessing
from xml.etree.ElementTree import canonicalize

import rdkit
from rdkit.Chem import AllChem

import selfies 
from selfies import encoder, decoder

from .utils import get_selfies_chars

def mutate_smiles(smi):
    """
    Given a smile, make random changes to the molecule.
    Operations to the molecules are:
    1) Adding a benzene ring
    2) Deleting a benzene ring
    3) Breaking apart a molecule into its constituent pieces
    """
    print("entering mutate_smiles in mutate_smirks.py")

    add_ring_smirks = "[c;H1:1][c;H1:2]>>[c:1]1[c:3][c:4][c:5][c:6][c:2]1"
    add_rxn = AllChem.ReactionFromSmarts(add_ring_smirks)

    reacts= [AllChem.MolFromSmiles(smi)]
    print(f"Reactant: {smi}")
    
    products = add_rxn.RunReactants(reacts)    # will be a 2D array

    # pick a random molecule as a result
    print(f"Products: {products}")
    result = random.choice(products)
    print(result)

    print("exiting mutate_smiles in mutate_smirks.py")
    return AllChem.MolToSmiles(result[0], canonical=True)


def mutate_sf(sf_chars, alphabet, num_sample_frags, base_alphabet = None):
    """
    Given a list of SELFIES alphabets, make random changes to the molecule using 
    alphabet. Opertations to molecules are character replacements, additions and deletions. 

    Parameters
    ----------
    sf_chars : (list of string alphabets)
        List of string alphabets for a SELFIE string.
    alphabet : (list of SELFIE strings)
        New SELFIES characters are added here and sampled.
    num_sample_frags: (int)
        Number of randomly sampled SELFIE strings.
    base_alphabet: (list of SELFIE strings)
        Main alphabet that will be appended with the introduced characters above.
        If none, use the semantic robust alphabet.

    Returns
    -------
    Mutated SELFIE string.

    """
    print("entering mutate_sf in mutate_smirks.py")
    if base_alphabet is None:
        base_alphabet = list(selfies.get_semantic_robust_alphabet())
    random_char_idx = random.choice(range(len(sf_chars)))
    choices_ls = [1, 2, 3]  # 1 = replacement; 2 = addition; 3=delete
    mutn_choice = choices_ls[
        random.choice(range(len(choices_ls)))
    ]  # Which mutation to do:

    if alphabet != []:
        alphabet = random.sample(alphabet, num_sample_frags) + base_alphabet
    else:
        alphabet = base_alphabet

    # Mutate character:
    if mutn_choice == 1:
        random_char = alphabet[random.choice(range(len(alphabet)))]
        change_sf = (
            sf_chars[0:random_char_idx]
            + [random_char]
            + sf_chars[random_char_idx + 1 :]
        )

    # add character:
    elif mutn_choice == 2:
        random_char = alphabet[random.choice(range(len(alphabet)))]
        change_sf = (
            sf_chars[0:random_char_idx] + [random_char] + sf_chars[random_char_idx:]
        )

    # delete character:
    elif mutn_choice == 3:
        if len(sf_chars) != 1:
            change_sf = sf_chars[0:random_char_idx] + sf_chars[random_char_idx + 1 :]
        else:
            change_sf = sf_chars

    return "".join(x for x in change_sf)






if __name__ == "__main__":
    molecules_here = [
        "CCC",
        "CCCC",
        "CCCCC",
        "CCCCCCCC",
        "CS",
        "CSSS",
        "CSSSSS",
        "CF",
        "CI",
        "CBr",
        "CSSSSSSSSSSSS",
        "CSSSSSSSSSC",
        "CSSSSCCSSSC",
        "CSSSSSSSSSF",
        "SSSSSC",
    ]
    # A = get_mutated_smiles(
    #     molecules_here, alphabet=["[C]"] * 500, num_sample_frags=200, space="Explore"
    # )

