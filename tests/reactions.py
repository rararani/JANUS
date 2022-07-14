from os.path import dirname, join
import random

from rdkit import Chem
from rdkit.Chem import rdChemReactions, Draw, AllChem

current_dir = dirname(__file__)
file_path = join(current_dir, "./DATA/acene_smiles.txt")

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

benzene = 'C1=CC=CC=C1'
# napthalene = 'c1ccc2ccccc2c1'

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

    # alright now let's just go all out and generate maybe 150 reaction smiles
    smi = benzene
    reactants = [AllChem.MolFromSmiles(smi)]
    file = open(file_path, "w")
    file.write(benzene + "\n")

    # for i in range(200):
    #     products = reaction.RunReactants(reactants)
    #     mol = random.choice(products)
    #     smi = AllChem.MolToSmiles(mol[0])
    #     file.write(smi + "\n")

    #     reactants = [Chem.MolFromSmiles(smi)]

    # this one works all the time
    for i in range(200):
        products = reaction.RunReactants(reactants)
        mol = products[0][0]
        smi = AllChem.MolToSmiles(mol)
        file.write(smi + "\n")

        reactants = [Chem.MolFromSmiles(smi)]
        print(i)



