import pandas as pd
import os
import argparse
from module import *
import shutil
import time


CURR_DIR = os.getcwd()
COLUMNS = ['scaffold', 'jgi', 'type', 'start', 'end', 'direction1', 'direction2', 'phase', 'infos']
GROUPBY_ATTR = 'name'
TEXT_WIDTH = 80

if __name__ == "__main__":
    t_start = time.time()

    # Create the argument parser
    parser = argparse.ArgumentParser(description="Argument parsing.")

    # Add arguments
    parser.add_argument('--gff', '-g', type=str, required=True, help="Path of the input .gff file")
    parser.add_argument('--species', '-s', type=str, required=True, help="Species name")
    parser.add_argument('--fasta', '-f', type=str, required=True, help="Path of the input fasta file")
    parser.add_argument('--outtype', '-t', choices=['protein', 'dna', 'both'], default='prot', help="Format of output fasta files")
    parser.add_argument('--outdir', '-o', type=str, default=CURR_DIR, help="Saving directory (current directory by default)")

    # Parse the arguments
    args = parser.parse_args()

    # Raise input file error 
    if (not os.path.exists(args.gff)) or (not args.gff.endswith('.gff')):
        raise FileNotFoundError(f"GFF file does not exist or the file is not .gff")

    if (not os.path.exists(args.fasta)) or (not args.fasta.endswith('.fasta')):
        raise FileNotFoundError(f"Fasta file does not exist or the file is not .fasta")
    
    # print(f"Input file: {args.input}")
    # print(f"Output file: {args.output}")
    
    # df_combined = pd.read_csv("df_combined.csv") 
    # df_gff = pd.read_csv("df_gff.csv")    

    # with open("fasta_dict.json", "r") as f:
    #     fasta_dict = json.load(f)

    print("=====LOADING GFF DATATABLE=======")
    df_gff = load_gff(file_path=args.gff, columns=COLUMNS)
    print("=====> DONE LOADING GFF DATATABLE")

    print("=====LOADING DNA FASTA FILE=======")
    fasta_dict = load_fasta(fasta_path=args.fasta)
    print("=====>DONE LOADING DNA FASTA FILE=======")

    # print(df)
    print("=====EXTRACTING CDS COORDINATES=======")
    df_coordinates = get_coordinates(df_gff=df_gff, groupby_attr=GROUPBY_ATTR)
    print("=====>DONE EXTRACTING CDS COORDINATES=======")

    print("=====FINDING THE TRUTH M POSITION=======")
    # df_combined = find_M(df=df_coordinates, fasta_dict=fasta_dict)  
    df_combined = findM_parallel(df=df_coordinates, fasta_dict=fasta_dict)  
    print("=====>DONE FINDING THE TRUTH M POSITION=======")

    print("=====EXPORTING TO CSV=======")
    df_combined.to_csv(f'{args.species}_corrected.csv')
    # export to new gff
    df_combined.to_csv(os.path.join(args.outdir, f'{args.species}_all_gene_corrected.gff'), sep='\t', header=False, index=False)
    print("=====>DONE EXPORT TO CSV=======")

    print("=====WRITTING FASTA FILES=======")
    # temporal directory for processing
    process_dir = os.path.join(args.outdir, 'processed_fasta')
    os.makedirs(process_dir, exist_ok=True)
    # writting out the fasta files
    writting_parallel(
        species=args.species,
        out_dir=args.outdir,
        processed_dir=process_dir, 
        df_combined=df_combined, 
        df_gff=df_gff,
        fasta_dict=fasta_dict, 
        outtype=args.outtype, 
        text_width=TEXT_WIDTH)
    # delete processed redundancies
    shutil.rmtree(process_dir)
    # os.rmdir(process_dir)
    print("=====>DONE WRITTING FASTA FILES=======")
    print('Processed time (s): ', int(time.time() - t_start))

    del df_gff, fasta_dict, df_coordinates, df_combined
