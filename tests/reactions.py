from os.path import dirname, join
import random

from rdkit import Chem
from rdkit.Chem import rdChemReactions, Draw, AllChem
from sympy import product

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

def generate_molecules(smi, rxn, num_generations=10):
    new_mols = []
    mol = AllChem.MolFromSmiles(smi)
    reactants = [mol]
    for _ in range(num_generations):
        products = rxn.RunReactants(reactants)
        new_mols.extend(products)
        if num_generations > 1:
            for product in products:
                new_mols.extend(generate_molecules(AllChem.MolToSmiles(product[0]), rxn, num_generations=1))
    
    return new_mols

def load_file_with_new_molecules(file_path, smi, rxn):
    file = open(file_path, "w")
    molecules = generate_molecules(smi, rxn)
    for molecule in molecules:
        file.write(AllChem.MolToSmiles(molecule[0]) + "\n")


benzene = 'C1=CC=CC=C1'
napthalene = 'c1ccc2ccccc2c1'
anthracene = 'C1=CC=C2C=C3C=CC=CC3=CC2=C1'
phenanthrene = 'C1=CC=C2C(=C1)C=CC3=CC=CC=C32'
unknown_mol = 'c1ccc2cc3cc4cc5cc6cc7cc8cc9cc%10cc%11ccccc%11cc%10cc9cc8cc7cc6cc5cc4cc3cc2c1'
unknown_mol2 = 'c1ccc2c(c1)cc1cc3ccc4c5ccccc5c5c6c7ccccc7c7ccccc7c6c6cc7cccc8c7c7c6c5c4c3c7c1c28'

add_nitrogen_smirks = "[#6&x2:1]>>[#7:1]"  # to find nitrogen
nitrogen_rxn = AllChem.ReactionFromSmarts(add_nitrogen_smirks)

# add a benzene ring
add_smirks = "[c;H1:1][c;H1:2]>>[c:1]1[c:3][c:4][c:5][c:6][c:2]1"
add_reaction = AllChem.ReactionFromSmarts(add_smirks)

fuse_smirks = "[c:1][c;H1:2][c;H1:3][c:4].[c:5][c;H1:6][c;H1:7][c:8]>>[c:1][c:2]([c:5])[c:3]([c:4])[c:8]"
fuse_rxn = AllChem.ReactionFromSmarts(fuse_smirks)

break_smirks = "[c;R2:1]1([c:2][c:3][c:4][c:5]2)[c;R2:6]2[c:7][c:8][c:9][c:10]1>>[c:1]1[c:2][c:3][c:4][c:5][c:6]1.[c:11]1[c:12][c:7][c:8][c:9][c:10]1"
break_reaction = AllChem.ReactionFromSmarts(break_smirks)

# fill out bay area
bay_smirks = "[c:1]1([c:2]([c;H1:3][c:4][c:5][c:6]2)[c:7]2[c:8][c:9]3)[c:10]3[c:11][c:12][c:13][c;H1:14]1>>[c:1]1([c:2]([c:3]([c:15][c:16]4)[c:4][c:5][c:6]2)[c:7]2[c:8][c:9]3)[c:10]3[c:11][c:12][c:13][c:14]14"
bay_reaction = AllChem.ReactionFromSmarts(bay_smirks)

if __name__ == "__main__":

    # add_nitro_rxn = AllChem.ReactionFromSmarts(add_nitrogen_smirks)
    # reactants = [Chem.MolFromSmiles(napthalene)]
    # products = add_nitro_rxn.RunReactants(reactants)
    # Draw.MolToFile(products[0][0], "nitrogen_rxn.png")

    smile = napthalene
    mol = AllChem.MolFromSmiles(smile)
    Draw.MolToFile(mol, "naphthalene.png")

    reactants = [mol, mol]
    products = fuse_rxn.RunReactants(reactants)
    print(len(products))
    Draw.MolToFile(products[0][0], 'fuse.png')

    # load_file_with_new_molecules(file_path, napthalene, add_reaction)

    # mol = AllChem.MolFromSmiles(unknown_mol2)
    # Draw.MolToFile(mol, "unknown_mol.png")

    # anthracene_mol = AllChem.MolFromSmiles(anthracene)
    # reactants = [anthracene_mol]
    # products = break_reaction.RunReactants(reactants)
    # mol = products[1][0]
    # print(len(products))
    # Draw.MolToFile(mol, "molecule.png")
    
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

    # file = open(file_path, "w")
    # # this one works all the time
    # # for j in range(20):
    # smi = benzene
    # reactants = [AllChem.MolFromSmiles(smi)]
    # # file.write(benzene + "\n")
    # for i in range(10):
    #     products = reaction.RunReactants(reactants)
    #     for mol in products:
    #         smi = AllChem.MolToSmiles(mol[0])
    #         file.write(smi + "\n")

    #     reactants = [Chem.MolFromSmiles(smi)]
    #     print(i)

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

    # napthalene_mol = AllChem.MolFromSmiles(napthalene)
    # reactants = [napthalene_mol]
    # add_rxn = AllChem.ReactionFromSmarts(add_smirks)
    # products = add_rxn.RunReactants(reactants)

    # print(len(products))

    # bay_rxn = AllChem.ReactionFromSmarts(bay_smirks)
    # print(products[1][0])
    # smi = AllChem.MolToSmiles(products[1][0])
    # mol = AllChem.MolFromSmiles(smi)
    # Draw.MolToFile(mol, "bay_mol.png")

    # reactants = [mol]
    # products = bay_rxn.RunReactants(reactants)

    # print(f"The # of products after bay fill operation: {len(products)}")

    # for i in range(len(products)):
    #     Draw.MolToFile(products[i][0], f"molecule{i}.png")

    # mol = AllChem.MolFromSmiles(unknown_mol)
    # Draw.MolToFile(mol, "molecule.png")