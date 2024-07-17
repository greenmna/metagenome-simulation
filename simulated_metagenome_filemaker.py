#!/usr/bin/env python

# Libraries
import pandas as pd
import argparse as ap
import os
import zipfile

# Arguments
parser = ap.ArgumentParser(description="""Generate files needed for simulation of short- and/or long-read metagenomic sequence data.""", 
                           formatter_class=ap.RawTextHelpFormatter)
# Arguments
parser.add_argument("-rf", "--reference-file", dest="ref", required=False, 
                    type=open, help="NCBI Dataset reference file.")
parser.add_argument("-ng", "--metagenome-size", dest="ngenomes", required=False, 
                    default=10, type=int, help="Number of genomes to have within a metagenome. Must be an integer greater than 1. (default: 10)")
parser.add_argument("-ns", "--sample-number", dest="nsamples", required=False, 
                    default=1, type=int, help="Number of samples to simulate for a given metagenome. Must be an integer of at least 1. (default: 1)")
parser.add_argument("-nm", "--number-metagenomes", dest="nmetagenomes", required=False, 
                    default=1, type=int, help="Number of metagenomes to simulate. Must be an integer of at least 1 (default: 1)")
parser.add_argument("-us", "--unique-species", dest="unique", required=False, 
                    type=open, help='Limits metagenomes to one instance of each species. If not possible (reference file is too limited), ' 
                    'the warning "Generation of a metagenome with unique species is not possible. Please check your reference file for "'
                    'sufficient, unique species or decrease the number of genomes used in your simulated metagenome(s)')
# set the type of DNA (if bacterial, assumed all are circular; otherwise, assumed linear (eukaryotes,viruses,etc.))
# set abundance variability (even/variable)
# set the quantity of data to generate (in gigabase-pairs)
args = parser.parse_args("")  # The "" is temporary for Jupyter Notebook. Remove later
ref = args.ref
ngenomes = args.ngenomes
nsamples = args.nsamples
unique = args.unique
nmetagenomes=args.nmetagenomes

# Functions
def strain_exclusive(metagenome):
    """Checks whether the user wants to avoid multiple of the same species within a metagenome.

    Args:
        sim_metagenome (DataFrame): metagenome created during process of input file creation
    
    Returns:
        unique_metagenome (DataFrame): updated metagenome where all organisms differ at the species level
    """
    # Find entries with duplicate values in the genus-species column
    
    # Drop n-1 of the rows with duplicate values (e.g. if 2, then drop 1)
    
    # Attempt to add in a new organism from temp_metagenome
    
    # Check if entries are unique, and repeat steps 1-3 until all genus-species are unique
    max_attempts = 100
    attempts = 0
    while metagenome['genus-species'].duplicated().any() and attempts < max_attempts:
        tmp_df = metagenome.drop_duplicates(subset=["genus-species"])
        num_organisms_to_replace = ngenomes - len(tmp_df)
        metagenome = pd.concat([tmp_df, ref.sample(n=num_organisms_to_replace)])  # Replace temp_report with ref from args.ref
        attempts += 1
    
    if attempts == max_attempts and metagenome['genus-species'].duplicated().any():
        print("Generation of a metagenome with unique species is not possible. "
              "Please check your reference file for a sufficient number of unique species "
              "or decrease the number of genomes used in your simulated metagenome(s)")
        return None  # In the main script, if strain_exclusive is None, then exit the code
        
    return metagenome

# Main
def main():
    parent_dir = os.path.abspath('./')  # Declare parent directory
    for i in range(0, nmetagenomes + 1):
        os.mkdir(f"metagenome_{nmetagenomes}")  # Create the directory for a simulated metagenome
        os.chdir(f"metagenome_{nmetagenomes}")  # Change to that directory

if __name__ == '__main__':
    main()