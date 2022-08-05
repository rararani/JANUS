from os.path import dirname, join
from coulson import coulson

from janus import JANUS, utils
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, RDConfig, Descriptors
RDLogger.DisableLog("rdApp.*")

import selfies

current_dir = dirname(__file__)

def fitness_function(smi: str) -> float:
    """ User-defined function that takes in individual smiles
    and outputs a fitness value.
    """
    # calculatig the % TRE (Topological Resonance Energy) of a molecule
    # currently defined for molecules containing ONLY C and H
    print("entering tre fitness function")
    print(f"Smiles:{smi}\n")
    mol = Chem.MolFromSmiles(smi)
    input_data, _ = coulson.interface.process_rdkit_mol(mol)
    huckel_matrix, electrons = coulson.huckel.prepare_huckel_matrix(
        input_data.atom_types, input_data.connectivity_matrix
    )
    tre, p_tre = coulson.graph_aromaticity.calculate_tre(input_data.connectivity_matrix, sum(electrons))
    print(f"exiting tre fitness for: {smi}\n")
    return -p_tre

def custom_filter(smi: str):
    """ Function that takes in a smile and returns a boolean.
    True indicates the smiles PASSES the filter.
    """
    # smiles length filter
    if len(smi) > 81 or len(smi) == 0:
        return False
    else:
        return True

if __name__ == "__main__":

    population_path = join(current_dir, "./DATA/acene_smiles.txt")
    params_path = join(current_dir, "./default_params.yml")

    print("Start!")

    # all parameters to be set, below are defaults
    params_dict = {
        # Number of iterations that JANUS runs for
        "generations": 50,

        # The number of molecules for which fitness calculations are done,
        # exploration and exploitation each have their own population
        "generation_size": 20,

        # Number of molecules that are exchanged between the exploration and exploitation
        "num_exchanges": 5,

        # Callable filtering function (None defaults to no filtering)
        "custom_filter": custom_filter,

        # Fragments from starting population used to extend alphabet for mutations
        "use_fragments": True,

        # An option to use a classifier as selection bias
        "use_classifier": True,
    }

    # Set your SELFIES constraints (below used for manuscript)
    default_constraints = selfies.get_semantic_constraints()
    new_constraints = default_constraints
    new_constraints['S'] = 2
    new_constraints['P'] = 3
    selfies.set_semantic_constraints(new_constraints)  # update constraints

    print("going to create agent now")

    # Create JANUS object.
    agent = JANUS(
        work_dir = 'RESULTS',                                   # where the results are saved
        fitness_function = fitness_function,                    # user-defined fitness for given smiles
        start_population = population_path,   # file with starting smiles population
        **params_dict
    )

    p_tre = fitness_function('C1=CC=CC=C1')
    print(p_tre)

    # Alternatively, you can get hyperparameters from a yaml file
    # Descriptions for all parameters are found in default_params.yml
    # params_dict = utils.from_yaml(
    #     work_dir = 'RESULTS',
    #     fitness_function = fitness_function,
    #     start_population = population_path,
    #     yaml_file = params_path,       # default yaml file with parameters
    #     **params_dict                           # overwrite yaml parameters with dictionary
    # )
    # agent = JANUS(**params_dict)

    # Run according to parameters
    print("running the agent now")
    agent.run()     # RUN IT!
