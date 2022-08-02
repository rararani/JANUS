from os.path import dirname, join
import random

from rdkit import Chem
from rdkit.Chem import rdChemReactions, Draw, AllChem

current_dir = dirname(__file__)
file_path = join(current_dir, "./DATA/acene_smiles.txt")
just_benz = join(current_dir, "./DATA/just_benzene.txt")

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def fill_file(smi, file_path):
    file = open(file_path, "w")
    for i in range(200):
        file.write(smi + "\n")


benzene = 'C1=CC=CC=C1'
napthalene = 'c1ccc2ccccc2c1'

# add a benzene ring
smirks = "[c;H1:1][c;H1:2]>>[c:1]1[c:3][c:4][c:5][c:6][c:2]1"
reaction = AllChem.ReactionFromSmarts(smirks)

# fill out bay area
smirks = "[c:1]1([c:2]([c;H1:3][c:4][c:5][c:6]2)[c:7]2[c:8][c:9]3)[c:10]3[c:11][c:12][c:13][c;H1:14]1>>[c:1]1([c:2]([c:3]([c:15][c:16]4)[c:4][c:5][c:6]2)[c:7]2[c:8][c:9]3)[c:10]3[c:11][c:12][c:13][c:14]14"

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
    # smi = benzene
    # reactants = [AllChem.MolFromSmiles(smi)]
    # file = open(file_path, "w")
    # file.write(benzene + "\n")

    file = open(file_path, "w")
    # this one works all the time
    # for j in range(20):
    smi = benzene
    reactants = [AllChem.MolFromSmiles(smi)]
    # file.write(benzene + "\n")
    for i in range(10):
        products = reaction.RunReactants(reactants)
        mol = products[0][0]
        smi = AllChem.MolToSmiles(mol)
        file.write(smi + "\n")

        reactants = [Chem.MolFromSmiles(smi)]
        print(i)

    # # just benzene
    # fill_file(benzene, just_benz)


    # goal: write out the "delete a benzene ring" operation in SMIRKS

    # anthracene = "C1=CC=C2C=C3C=CC=CC3=CC2=C1"
    # benzanthracene = "C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32"

    # this one will stricly reduce whatever you have to one benzene ring
    # smirks = "[c:1]1[c:3][c:4][c:5][c:6][c:2]1>>[c:1]1[c:2][c:3][c:8][c:9][c:4]1"

    # smirks = "[c;R2:1][c;R2:2][c;R1:3][c;R1:4][c;R1:5][c;R1:6]>>[c:1][c:2]" # THIS WORKS FOR DELETION

    # orig_mol = mol_with_atom_index(AllChem.MolFromSmiles(benzanthracene))

    # ri = orig_mol.GetRingInfo()
    # rings = ri.AtomRings()
    # print(rings)
    # del_ring = random.choice(rings)

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

    # rxn = AllChem.ReactionFromSmarts(smirks)
    # reactants = [AllChem.MolFromSmiles(benzanthracene)]
    # products = rxn.RunReactants(reactants)
    
    # print(f"Original: {benzanthracene}\nResult: {AllChem.MolToSmiles(products[0][0])}")
    
    # prod_mol = mol_with_atom_index(products[0][0])
    # prod2_mol = mol_with_atom_index(products[1][0])
    
    # Draw.MolToFile(orig_mol, "molecule.png")
    # Draw.MolToFile(prod_mol, "molecule1.png")
    # Draw.MolToFile(prod2_mol, "molecule2.png")

    # Draw.MolToFile(AllChem.MolFromSmiles("c1ccc2ccccc2c1"), "molecule3.png")