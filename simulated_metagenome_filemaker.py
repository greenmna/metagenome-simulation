#!/usr/bin/env python

# Libraries
import numpy as np
import pandas as pd
import argparse as ap
import os
import pathlib
import subprocess
import shutil
import zipfile

# Silence warnings about chain assignment from Pandas
pd.options.mode.chained_assignment = None

# Arguments
parser = ap.ArgumentParser(description="""Generate files needed for simulation of short- and/or long-read metagenomic sequence data.""", 
                           formatter_class=ap.RawTextHelpFormatter)
# Arguments
parser.add_argument("-rf", "--reference-file", dest="ref", type=str, required=False, help="PATH to NCBI Dataset reference file. "
                    "This should be in the main directory where all simulated metagenomes will exist.")
parser.add_argument("-ng", "--metagenome-size", dest="ngenomes", required=False, 
                    default=10, type=int, help="Number of genomes to have within a metagenome. Must be an integer greater than 1. (default: 10)")
parser.add_argument("-ns", "--sample-number", dest="nsamples", required=False, 
                    default=1, type=int, help="Number of samples to simulate for a given metagenome. Must be an integer of at least 1. (default: 1)")
parser.add_argument("-nm", "--number-metagenomes", dest="nmetagenomes", required=False, 
                    default=1, type=int, help="Number of metagenomes to simulate. Must be an integer of at least 1 (default: 1)")
parser.add_argument("-us", "--unique-species", dest="unique", required=False, action='store_true', 
                    help='Limits metagenomes to one instance of each species. If not possible (reference file is too limited), ' 
                    'the warning "Generation of a metagenome with unique species is not possible. Please check your reference file for "'
                    'sufficient, unique species or decrease the number of genomes used in your simulated metagenome(s)')
parser.add_argument("-ce", "--conda-environment", dest="conda_env", required=False, 
                    type=str, default="ncbi_datasets", 
                    help="Specify the conda environment that ncbi_datasets_cli is accessed in. (default: ncbi_datasets)")
# set the type of DNA (if bacterial, assumed all are circular; otherwise, assumed linear (eukaryotes,viruses,etc.))
# set abundance variability (even/variable)
# set the quantity of data to generate (in gigabase-pairs)
args = parser.parse_args()  # The "" is temporary for Jupyter Notebook. Remove later
ref = args.ref
ngenomes = args.ngenomes
nsamples = args.nsamples
unique = args.unique
nmetagenomes=args.nmetagenomes
conda_env = args.conda_env

parent_dir = os.path.abspath('./')  # Parent Directory
ref_path = os.path.abspath(ref)  # Must be in the Parent Directory

# Create DataFrame from --reference-file input
reference = pd.read_csv(ref_path, sep='\t')
translation_table = dict.fromkeys(map(ord, "']["), None)
reference["Organism Name"] = reference["Organism Name"].apply(lambda x: x.translate(translation_table))
reference.rename(columns={"Assembly BioSample Strain ":"Assembly BioSample Strain"}, inplace=True)  # Should be fixed by NCBI...
reference['genus-species'] = reference['Organism Name'].apply(lambda x: ' '.join(x.split(' ')[0:2]))  # Create genus-species column

# Functions
def strain_exclusive(metagenome):
    """Checks whether the user wants to avoid multiple of the same species within a metagenome.

    Args:
        metagenome (DataFrame): metagenome created during process of input file creation
        r (DataFrame): this is the reference created before defining function
        ng (int): the number of genomes to have within the metagenome
    
    Returns:
        unique_metagenome (DataFrame): updated metagenome where all organisms differ at the species level
    """
    # Find entries with duplicate values in the genus-species column
    
    # Drop n-1 of the rows with duplicate values (e.g. if 2, then drop 1)
    
    # Attempt to add in a new organism from temp_metagenome
    
    # Check if entries are unique, and repeat steps 1-3 until all genus-species are unique
    max_attempts = 1000
    attempts = 0
    while metagenome.duplicated(subset=["genus-species"]).any() and attempts < max_attempts:
        metagenome.drop_duplicates(subset=["genus-species"], inplace=True)
        num_organisms_to_replace = ngenomes - len(metagenome)
        if num_organisms_to_replace <= 0:
            break
        
        new_entries = reference.sample(n=num_organisms_to_replace)
        
        metagenome = pd.concat([metagenome, new_entries])  # Replace temp_report with ref from args.ref
        attempts += 1
    
    if metagenome['genus-species'].duplicated().any() and attempts == max_attempts:
        print("Generation of a metagenome with unique species is not possible. "
              "Please check your reference file for a sufficient number of unique species "
              "or decrease the number of genomes used in your simulated metagenome(s)")
        return None  # In the main script, if strain_exclusive is None, then exit the code
        
    return metagenome


def download_genomes(ref_df):
    synth_metagenome = ref_df.sample(n=ngenomes)
    if unique:
        strain_exclusive(synth_metagenome)  # Generate metagenome with unique species
        accessions = ','.join(list(set(list(synth_metagenome["Assembly Accession"]))))
        download_command = f"conda run -n {conda_env} datasets download genome accession {accessions} --assembly-level complete --assembly-source RefSeq --exclude-atypical"
        subprocess.run(download_command, shell=True)
    
    else:
        accessions = ','.join(list(set(list(synth_metagenome["Assembly Accession"]))))
        download_command = f"conda run -n {conda_env} datasets download genome accession {accessions} --assembly-level complete --assembly-source RefSeq --exclude-atypical"
        subprocess.run(download_command, shell=True)
    
    synth_metagenome['parse'] = synth_metagenome['genus-species'].apply(lambda x: x.split(' '))
    return synth_metagenome

def extract_genomes(metagenome_dir):
    if os.path.isdir(f"metagenome_{nmetagenomes}_genomes"):
        # Make clean directory
        shutil.rmtree(f"metagenome_{nmetagenomes}_genomes")  
        os.mkdir(f"metagenome_{nmetagenomes}_genomes")
    
    # Unzip and delete zipped ncbi_datasets
    with zipfile.ZipFile(f"{metagenome_dir}/ncbi_dataset.zip", 'r') as genome_zip:
        genome_zip.extractall(f"{os.path.abspath('./')}/genome_data")
    
    genome_dirs = [str(d) for d in pathlib.Path(f"{metagenome_dir}/genome_data/ncbi_dataset/data").iterdir() if d.is_dir()]
    genome_files = [str(next(pathlib.Path(d).glob('*'))) for d in genome_dirs]
    for file in genome_files:
        shutil.move(file, f"metagenome_{nmetagenomes}_genomes")
    
    # Remove intermediary folder(s)
    shutil.rmtree(f"{metagenome_dir}/genome_data")
    
    return os.path.abspath(f"metagenome_{nmetagenomes}_genomes")


def create_genome_list(sim_metagenome_df, genome_location):
    genome_paths = [os.path.abspath(filepath) for filepath in os.listdir(genome_location)]
    tmp_genome_list = sim_metagenome_df[['Organism Name']]
    tmp_genome_list['genome_path'] = genome_paths
    return tmp_genome_list

def create_dna_list(sim_metagenome_df):
    tmp_dna_list = sim_metagenome_df[['Organism Name', 'Assembly Stats Number of Contigs']]
    tmp_dna_list = tmp_dna_list.loc[np.repeat(tmp_dna_list.index.values, tmp_dna_list['Assembly Stats Number of Contigs'])].reset_index(drop=True)
    counter = {}
    tmp_dna_list['DNA molecules'] = ''
    
    for index, row in tmp_dna_list.iterrows():
        organism = row['Organism Name']
        if organism not in counter:
            counter[organism] = 1
        else:
            counter[organism] += 1
        tmp_dna_list.at[index, 'DNA molecules'] = f'{organism}_dna_molecule_{counter[organism]}'
    
    tmp_dna_list.drop('Assembly Stats Number of Contigs', axis=1, inplace=True)
    return tmp_dna_list

# Main

def main():
    for i in range(nmetagenomes):
        if not os.path.isdir(f"metagenome_{nmetagenomes}"):
            os.mkdir(f"metagenome_{nmetagenomes}")  # Create the directory for a simulated metagenome

        os.chdir(f"metagenome_{nmetagenomes}")
        working_dir = os.path.abspath("./")
        
        # Generate simulated metagenome metadata
        sim_metagenome = download_genomes(reference)
        sim_metagenome.sort_values(by='Assembly Accession', inplace=True)
        sim_metagenome.to_csv(f"{working_dir}/metagenome_simulation_reference_file.tsv",sep="\t")
        
        # Extract downloaded information from NCBI
        genome_dir = extract_genomes(working_dir)
        
        # Generate genome_list
        genome_list = create_genome_list(sim_metagenome, genome_dir)
        genome_list.to_csv(f"{working_dir}/metagenome_{nmetagenomes}_genome_list.tsv", sep="\t", header=None, index=False)
        # Generate DNA_list
        dna_list = create_dna_list(sim_metagenome)
        dna_list.to_csv(f"{working_dir}/metagenome_{nmetagenomes}_dna_list.tsv", sep="\t", header=None, index=False)
        
        os.chdir(parent_dir)


if __name__ == '__main__':
    main()