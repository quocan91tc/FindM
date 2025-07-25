import pandas as pd
import os
import argparse
from module import *


CURR_DIR = os.getcwd()
COLUMNS = ['scaffold', 'jgi', 'type', 'start', 'end', 'direction1', 'direction2', 'phase', 'infos']
GROUPBY_ATTR = 'name'
TEXT_WIDTH = 80

if __name__ == "__main__":

    # Create the argument parser
    parser = argparse.ArgumentParser(description="Argument parsing.")

    # Add arguments
    parser.add_argument('--gff', '-g', type=str, required=True, help="Path of the input .gff file")
    parser.add_argument('--species', '-s', type=str, required=True, help="Species name")
    parser.add_argument('--fasta', '-f', type=str, required=True, help="Path of the input fasta file")
    parser.add_argument('--outtype', '-t', choices=['protein', 'dna', 'both'], default='prot', help="Format of output fasta files")
    parser.add_argument('--output', '-o', type=str, default=CURR_DIR, help="Saving directory (current directory by default)")

    # Parse the arguments
    args = parser.parse_args()

    # Raise input file error 
    if (not os.path.exists(args.gff)) or (not args.gff.endswith('.gff')):
        raise FileNotFoundError(f"GFF file does not exist or the file is not .gff")

    if (not os.path.exists(args.fasta)) or (not args.fasta.endswith('.fasta')):
        raise FileNotFoundError(f"Fasta file does not exist or the file is not .fasta")
    
    # print(f"Input file: {args.input}")
    # print(f"Output file: {args.output}")
    
    df_gff = load_gff(file_path=args.gff, columns=COLUMNS)
    fasta_dict = load_fasta(fasta_path=args.fasta)

    # print(df)
    df_coordinates = get_coordinates(df_gff=df_gff, groupby_attr=GROUPBY_ATTR)
    df_corrected = find_M(df=df_coordinates, fasta_dict=fasta_dict)  

    df_corrected.to_csv(f'{args.species}_corrected.csv')
    # export to new gff
    df_corrected.to_csv(os.path.join(args.output, f'{args.species}_all_gene_corrected.gff'), sep='\t', header=False, index=False)
    # export to new fasta
    get_fastas(df=df_corrected, fasta_dict=fasta_dict, out_dir=args.output, species=args.species, text_width=TEXT_WIDTH, output_type=args.outtype)
    
    del df_gff, fasta_dict, df_coordinates, df_corrected
