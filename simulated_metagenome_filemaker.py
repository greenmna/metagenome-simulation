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
import random
import sys

# Silence warnings about chain assignment from Pandas
pd.options.mode.chained_assignment = None

# Arguments
parser = ap.ArgumentParser(description="""Generate files needed for simulation of short- and/or long-read metagenomic sequence data.""", 
                           formatter_class=ap.RawTextHelpFormatter)
# Arguments
parser.add_argument("-rf", "--reference-file", dest="ref", type=str, required=True, help="PATH to NCBI Dataset reference file. "
                    "Preferably, this is in the directory where simulated metagenomes are created.")
parser.add_argument("-ng", "--metagenome-size", dest="ngenomes", required=False, 
                    default=10, type=int, help="Number of genomes to have within a metagenome. Must be an integer greater than 1. (default: 10)")
parser.add_argument("-ns", "--sample-number", dest="nsamples", required=False, 
                    default=1, type=int, help="Number of samples to simulate for a given metagenome. Must be an integer of at least 1. (default: 1)")
parser.add_argument("-nm", "--number-metagenomes", dest="nmetagenomes", required=False, 
                    default=1, type=int, help="Number of metagenomes to simulate. Must be an integer of at least 1 (default: 1)")
parser.add_argument("-us", "--unique-species", dest="unique", required=False, action='store_true', 
                    help='Limits metagenomes to one instance of each species.\n' \
                    'If not possible, the warning "Generation of a metagenome with unique species is not possible.\n' \
                    'Please check your reference file for a sufficient number of unique species or decrease the number of genomes used."')
parser.add_argument("-ce", "--conda-environment", dest="conda_env", required=False, 
                    type=str, default="ncbi_datasets", 
                    help="Specify the conda environment that ncbi_datasets_cli is accessed in. (default: ncbi_datasets)")
parser.add_argument("-av", "--abundance-variation", dest="vary_abundance", required=False, 
                    type=str, choices={"even", "random"}, default="random", 
                    help="Describes how abundances for organisms in a metagenome will be distributed.\n" \
                        "even: all organisms have the same number of reads\n" \
                        "random: abundance values are assigned at random\n" \
                        "(default: random)")
parser.add_argument("-rc", "--read-count", dest="read_count", required=False, default=650_000, type=int, 
                    help="Specify how many reads will be simulated (default: 650_000)")
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
vary_abundance = args.vary_abundance
read_count = args.read_count

parent_dir = os.path.abspath('./')  # Parent Directory

if not ref:
    os.path.join(parent_dir, "example_files", "bacteria_refseq_catalog.tsv")
else:
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
    while metagenome.duplicated(subset=["genus-species"]).any() and attempts < max_attempts and len(metagenome) < ngenomes:
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
    
    return synth_metagenome

def extract_genomes(metagenome_dir):
    if os.path.isdir(f"metagenome_{nmetagenomes}_genomes"):
        # Make clean directory
        shutil.rmtree(f"metagenome_{nmetagenomes}_genomes")  
        os.mkdir(f"metagenome_{nmetagenomes}_genomes")
    else:
        os.mkdir(f"metagenome_{nmetagenomes}_genomes")
    
    # Unzip and delete zipped ncbi_datasets
    with zipfile.ZipFile(f"{metagenome_dir}/ncbi_dataset.zip", 'r') as genome_zip:
        genome_zip.extractall(f"{os.path.abspath('./')}/genome_data")
    
    genome_dirs = [str(d) for d in pathlib.Path(f"{metagenome_dir}/genome_data/ncbi_dataset/data").iterdir() if d.is_dir()]
    genome_files = [str(next(pathlib.Path(d).glob('*'))) for d in genome_dirs]
    for file in genome_files:
        shutil.copy(file, f"metagenome_{nmetagenomes}_genomes")
    
    # Remove intermediary folder(s)
    shutil.rmtree(f"{metagenome_dir}/genome_data")
    
    return os.path.abspath(f"metagenome_{nmetagenomes}_genomes")


def create_genome_list(sim_metagenome_df, genome_location):
    genome_paths = [os.path.abspath(os.path.join(genome_location, filepath)) for filepath in os.listdir(genome_location)]
    tmp_genome_list = sim_metagenome_df[['Organism Name']]
    tmp_genome_list['genome_path'] = genome_paths
    return tmp_genome_list

def create_dna_list(sim_metagenome_df, current_metagenome):
    # Gather taxonomic information
    taxa_to_download = ",".join([t.astype(str) for t in list(sim_metagenome_df["Organism Taxonomic ID"].unique())])
    download_taxa_info = f"conda run -n {conda_env} datasets summary taxonomy taxon {taxa_to_download} --as-json-lines > " \
                         f"tmp_taxa_info.jsonl"
    subprocess.run(download_taxa_info, shell=True)
    jsonl_file = os.path.join(os.path.abspath("./"), "tmp_taxa_info.jsonl")
    
    # Current issues with dataformat make it fail to recogniez .jsonl files despite being able to as indicated on their documentation
    # Errors can be ignored if claiming format issues with input file (conda will also kick an error, but code will proceed without issue)
    format_taxa_info = f"conda run -n {conda_env} dataformat tsv taxonomy --inputfile {jsonl_file} --template tax-summary > " \
                       f"{parent_dir}/metagenome_{nmetagenomes}/metagenome_{nmetagenomes}_taxa_report.tsv"
    subprocess.run(format_taxa_info, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
    # Add taxonomic information to reference file
    taxa_report = pd.read_csv(f"{parent_dir}/metagenome_{current_metagenome}/metagenome_{current_metagenome}_taxa_report.tsv", sep="\t")
    
    tmp_dna_list = pd.merge(sim_metagenome_df, taxa_report, left_on=["Organism Taxonomic ID"], right_on=["Query"])   
    tmp_dna_list = tmp_dna_list[['Organism Name', 'Assembly Stats Number of Contigs', 'Superkingdom name']]

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

    tmp_dna_list['DNA molecule type'] = tmp_dna_list['Superkingdom name'].apply(lambda x: "circluar" if x=="Bacteria" else "linear")  # Not perfect...
    tmp_dna_list.drop(['Assembly Stats Number of Contigs', 'Superkingdom name'], axis=1, inplace=True)
    
    return tmp_dna_list

def create_abundance_list(sim_metagenome_df):
    # Col 1 is called "Size"
    # Col 2 - Col n is the number of reads to simulate (default 650_000)
    # NanoSim uses percentages, where BBMap uses whole numbers of reads
    # Given this, we need to adjust the percentages such that whole number are providable for BBMap automatically  
    cols = ["Size"]
    [cols.append(x) for x in range(0, nsamples)]
    abund = pd.DataFrame(columns=cols, dtype=str)
    abund["Size"] = sim_metagenome_df['Organism Name']
    
    if vary_abundance == "even":
        abund[cols[1:]] = 100 / len(abund.Size)
    elif vary_abundance == "random":
        for i in range(0, len(abund.columns) - 1):
            abund[i] = random.choices(range(1, read_count), k=ngenomes)
        abund = abund.apply(lambda x: (x / x.sum()) * 100 if x.dtype == "int64" else x)
    else:
        print('An improper string was given for the vary_abundance argument. Please specify either "even" or "random"')
        sys.exit()
        
    abund.rename(columns=dict([x, read_count] for x in range(0, nsamples)), inplace=True)
    return abund


# Main

def main():
    for i in range(1, nmetagenomes + 1):
        if not os.path.isdir(f"metagenome_{i}"):
            os.mkdir(f"metagenome_{i}")  # Create the directory for a simulated metagenome

        os.chdir(f"metagenome_{i}")
        working_dir = os.path.abspath("./")
        
        # Generate simulated metagenome metadata
        sim_metagenome = download_genomes(reference)
        sim_metagenome.sort_values(by='Assembly Accession', inplace=True)
        sim_metagenome.to_csv(f"{working_dir}/metagenome_simulation_reference_file.tsv",sep="\t")
        
        # Extract downloaded information from NCBI
        genome_dir = extract_genomes(working_dir)
        
        # Generate genome_list
        genome_list = create_genome_list(sim_metagenome, genome_dir)
        genome_list.to_csv(f"{working_dir}/metagenome_{i}_genome_list.tsv", sep="\t", header=None, index=False)
        # Generate DNA_list
        dna_list = create_dna_list(sim_metagenome, i)
        dna_list.to_csv(f"{working_dir}/metagenome_{i}_dna_list.tsv", sep="\t", header=None, index=False)
        # Generate abudance_list
        abundance_list = create_abundance_list(sim_metagenome)
        abundance_list.to_csv(f"{working_dir}/metagenome_{i}_abundance_list.tsv", sep="\t", index=False)
        
        os.chdir(parent_dir)


if __name__ == '__main__':
    main()