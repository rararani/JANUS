from os.path import dirname, join
import random

from rdkit import Chem
from rdkit.Chem import rdChemReactions, Draw, AllChem
from sympy import product

current_dir = dirname(__file__)
file_path = join(current_dir, "./DATA/acene_smiles.txt")

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

benzene = 'C1=CC=CC=C1'
napthalene = 'c1ccc2ccccc2c1'

# add a benzene ring
smirks = "[c;H1:1][c;H1:2]>>[c:1]1[c:3][c:4][c:5][c:6][c:2]1"
reaction = AllChem.ReactionFromSmarts(smirks)

if __name__ == "__main__":
    # print(AllChem.MolToSmiles(products[0][0]))
    # print(len(products[0]))
    # print(len(products))
    # smi = Chem.MolToSmiles(products[0][0])
    # for prod in products:
    #     print(Chem.MolToSmiles(prod[0]))
    # print("True")

    # orig_mol = mol_with_atom_index(AllChem.MolFromSmiles(benzene))
    # img = Draw.MolToFile(orig_mol, 'molecule.png')
    # mol = mol_with_atom_index(products[0][0])
    # img = Draw.MolToFile(mol, "molecule1.png")

    # # let's make one more extension

    # new_reactants = [Chem.MolFromSmiles(smi)]
    # products = reaction.RunReactants(new_reactants)
    # img = Draw.MolToFile(products[0][0], "mol2.png")

    # # alright now let's just go all out and generate maybe 150 reaction smiles
    # smi = benzene
    # reactants = [AllChem.MolFromSmiles(smi)]
    # file = open(file_path, "w")
    # file.write(benzene + "\n")

    # # this one works all the time
    # for i in range(200):
    #     products = reaction.RunReactants(reactants)
    #     mol = products[0][0]
    #     smi = AllChem.MolToSmiles(mol)
    #     file.write(smi + "\n")

    #     reactants = [Chem.MolFromSmiles(smi)]
    #     print(i)

    # goal: write out the "delete a benzene ring" operation in SMIRKS

    anthracene = "C1=CC=C2C=C3C=CC=CC3=CC2=C1"
    
    # this one will stricly reduce whatever you have to one benzene ring
    smirks = "[c:1]1[c:3][c:4][c:5][c:6][c:2]1>>[c:1]1[c:2][c:3][c:8][c:9][c:4]1"

    # smirks = "[c:1]1[c:3][c:4][c:5][c:6][c:2]1>>[c:1]1[c:3][c:4][c:5][c:6][c:2]1"

    orig_mol = mol_with_atom_index(AllChem.MolFromSmiles(anthracene))

    ri = orig_mol.GetRingInfo()
    rings = ri.AtomRings()
    print(rings)
    del_ring = random.choice(rings)

    # new_rings = []
    # for ring in rings:
    #     if ring != del_ring:
    #         new_rings.append(ring)

    # # construct SMARTS string preserving all carbon atoms except for del_ring
    # del_smarts = ""
    # num_rings = len(new_rings)
    # for i in range(new_rings):
    #     atoms = new_rings[]
    #     del_smarts += f"[c:{atoms[i]}]"
    #     if i == 0:
    #         del_smarts += "1"
    # del_smarts += "1"

    # smirks = "[c:1]1[c:3][c:4][c:5][c:6][c:2]1>>" + del_smarts
    # print(smirks)

    rxn = AllChem.ReactionFromSmarts(smirks)
    reactants = [AllChem.MolFromSmiles(anthracene)]
    products = rxn.RunReactants(reactants)

    print(f"Original: {anthracene}\nResult: {AllChem.MolToSmiles(products[0][0])}")

    prod_mol = mol_with_atom_index(products[0][0])

    Draw.MolToFile(orig_mol, "molecule.png")
    Draw.MolToFile(prod_mol, "molecule1.png")



