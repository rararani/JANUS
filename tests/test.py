from rdkit.Chem import AllChem, Draw

if __name__ == "__main__":
    mol = AllChem.MolFromSmiles('c1ccc2cc3cc4c(cc3cc2c1)c1cc2c3ccccc3c3c5ccc6c(ccc7ccc8cc9ccccc9cc8c76)c5c5cc6c(ccc7c6c6ccc8c9ccc%10ccccc%10c9c9ccccc9c8c6c6c8ccccc8c8c(ccc9ccc%10ccccc%10c98)c76)cc5c3c2cc1c1c2ccccc2c2c3ccccc3c3c5ccccc5ccc3c2c41')
    Draw.MolToFile(mol, "weirdmol.png")