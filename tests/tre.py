from cProfile import label
from os.path import dirname, join
from tokenize import PlainToken
from coulson import coulson

from janus import JANUS, utils
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, RDConfig, Descriptors
import matplotlib.pyplot as plt
import numpy as np
RDLogger.DisableLog("rdApp.*")

import selfies

current_dir = dirname(__file__)

def fitness_function(smi: str) -> float:
    """ User-defined function that takes in individual smiles
    and outputs a fitness value.
    """
    # calculatig the % TRE (Topological Resonance Energy) of a molecule
    # currently defined for molecules containing ONLY C and H
    # print("entering tre fitness function")
    # print(f"Smiles:{smi}\n")
    mol = Chem.MolFromSmiles(smi)
    input_data, _ = coulson.interface.process_rdkit_mol(mol)
    huckel_matrix, electrons = coulson.huckel.prepare_huckel_matrix( 
        input_data.atom_types, input_data.connectivity_matrix
    )
    tre, p_tre = coulson.graph_aromaticity.calculate_tre(huckel_matrix, sum(electrons), multiplicity=3)
    # print(f"exiting tre fitness for: {smi}\n")
    return p_tre

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
        "generations": 20,

        # The number of molecules for which fitness calculations are done,
        # exploration and exploitation each have their own population
        "generation_size": 10,

        # Number of molecules that are exchanged between the exploration and exploitation
        "num_exchanges": 5,

        # Callable filtering function (None defaults to no filtering)
        "custom_filter": custom_filter,

        # Fragments from starting population used to extend alphabet for mutations
        "use_fragments": False,

        # An option to use a classifier as selection bias
        "use_classifier": False,
    }

    # params_dict = {
    #     "generations": 5,
    #     "generation_size": 2,
    #     "num_exchanges": 1,
    #     "custom_filter": custom_filter,
    #     "use_fragments": False,
    #     "use_classifier": False

    # }

    # Set your SELFIES constraints (below used for manuscript)
    default_constraints = selfies.get_semantic_constraints()
    new_constraints = default_constraints
    new_constraints['S'] = 2
    new_constraints['P'] = 3
    selfies.set_semantic_constraints(new_constraints)  # update constraints

    print("going to create agent now")

    x = []
    for i in range(20):
        x.append(i)
    
    array1 = []
    array2 = []
    for i in range(1):

        # Create JANUS object.
        agent = JANUS(
            work_dir = 'RESULTS',                                   # where the results are saved
            fitness_function = fitness_function,                    # user-defined fitness for given smiles
            start_population = population_path,   # file with starting smiles population
            **params_dict
        )

        # p_tre = fitness_function('C1=CC=CC=C1')
        # print(p_tre)

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
        top_explr_mols, top_local_mols = agent.run()     # RUN IT!
        
        array1.append(top_explr_mols)
        array2.append(top_local_mols)

    
    # top_explr_mols = np.array(array1)
    # top_local_mols = np.array(array2)

    # top_explr_mols = np.sum(top_explr_mols, axis=0) / 20
    # top_local_mols = np.sum(top_local_mols, axis=0) / 20

    # plt.plot(x, top_explr_mols, 'b', label="Explore Stage")
    # plt.plot(x, top_local_mols, 'r', label="Exploit Stage")
    # plt.xticks(x)
    # plt.xlabel("Generation")
    # plt.ylabel("% Topological Resonance Energy (TRE)")
    # plt.title("% TRE per Generation for Ground-state PAHs")
    # plt.legend()
    # plt.savefig("per_tre_ground_graph_min.png")

    top_explr_mols = np.array(top_explr_mols)
    top_local_mols = np.array(top_local_mols)

    # top_explr_mols = top_explr_mols * -1
    # top_local_mols = top_local_mols * -1


    plt.plot(x, top_explr_mols, 'b', label="Explore Stage")
    plt.plot(x, top_local_mols, 'r', label="Exploit Stage")
    plt.xticks(x)
    plt.xlabel("Generation")
    plt.ylabel("% Topological Resonance Energy (TRE)")
    plt.title("% TRE per Generation for Ground-state PAHs")
    plt.legend()
    plt.savefig("per_tre_ground_graph_max.png")


    # plt.plot(x, top_explr_mols, label="Explore Stage")
    # plt.plot(x, top_local_mols, label="Exploit Stage")
    # plt.xlabel("Generation")
    # plt.ylabel("% Topological Resonance Energy (TRE)")
    # plt.title("% TRE per Generation for Excited-state PAHs")
    # plt.legend("Explore Stage", "Exploit Stage")
    # plt.show()


