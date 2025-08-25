import pandas as pd
import os
import argparse
from module import *
import shutil
import time


CURR_DIR = os.getcwd()
COLUMNS_V1 = ['scaffold', 'jgi', 'feature_type', 'start', 'end', 'score', 'direction', 'phase', 'infos']
COLUMNS_V3 = ['scaffold', 'source', 'feature_type', 'start', 'end', 'score', 'direction', 'phase', 'infos']
GROUPBY_ATTR_V1 = 'name'
GROUPBY_ATTR_V3 = 'Parent'

FEATURE_TYPE_V1 = 'type'
TEXT_WIDTH = 80

if __name__ == "__main__":
    t_start = time.time()

    # Create the argument parser
    parser = argparse.ArgumentParser(description="Argument parsing.")

    # Add arguments
    parser.add_argument('--gff', '-g', type=str, required=True, help="the input .gff file path")
    parser.add_argument('--species', '-s', type=str, required=True, help="Species name")
    parser.add_argument('--version', '-v', type=str, choices=['1', '3'], default='v1', help="Input GFF version (v1 as default)")
    parser.add_argument('--fasta', '-f', type=str, required=True, help="the input fasta file path (.fasta and .fa accepted)")
    parser.add_argument('--outdir', '-o', type=str, default=CURR_DIR, help="Saving directory (current directory as default)")

    # Parse the arguments
    args = parser.parse_args()

    # Raise input file error 
    if (not os.path.exists(args.gff)) or (not args.gff.endswith('.gff')):
        raise FileNotFoundError(f"GFF file does not exist or the file is not .gff")

    if (not os.path.exists(args.fasta)) or ((not args.fasta.endswith('fasta')) and (not args.fasta.endswith('fa'))):
        raise FileNotFoundError(f"Fasta file does not exist or the file is not .fasta or .fa")
    
    # ========= GFF VERSION 1 =============
    if args.version == '1':
        print("=====LOADING GFF TABLE=======")
        df_gff = load_gff(file_path=args.gff, columns=COLUMNS_V1)
        print("=====> DONE LOADING GFF TABLE")

        print("=====LOADING DNA FASTA FILE=======")
        fasta_dict = load_fasta(fasta_path=args.fasta)
        print("=====>DONE LOADING DNA FASTA FILE=======")

        # print(df)
        print("=====EXTRACTING CDS COORDINATES=======")
        df_cds_1st = get_coordinates(df_gff=df_gff, groupby_attr=GROUPBY_ATTR_V1)
        print("=====>DONE EXTRACTING CDS COORDINATES=======")

        print("=====FINDING THE TRUTH M POSITION=======")
        # df_combined = find_M(df=df_coordinates, fasta_dict=fasta_dict)  
        df_combined = findM_parallel(df=df_cds_1st, fasta_dict=fasta_dict)  
        print("=====>DONE FINDING THE TRUTH M POSITION=======")

        print("=====EXPORTING TO CSV=======")
        df_combined.to_csv(f'{args.species}_corrected.csv')
        # export to new gff
        df_combined.to_csv(os.path.join(args.outdir, f'{args.species}_all_gene_corrected.gff'), sep='\t', index=False)
        print("=====>DONE EXPORT TO CSV=======")

        print("=====WRITTING FASTA FILES=======")
        # temporal directory for processing
        process_dir = os.path.join(args.outdir, 'processed_fasta')
        os.makedirs(process_dir, exist_ok=True)
        # writting out the fasta files
        writting_parallel(
            ver=args.version,
            species=args.species,
            out_dir=args.outdir,
            processed_dir=process_dir, 
            df_combined=df_combined, 
            df_gff=df_gff,
            fasta_dict=fasta_dict, 
            outtype=args.outtype, 
            text_width=TEXT_WIDTH,
            groupby_col=GROUPBY_ATTR_V1)
        # delete processed redundancies
        shutil.rmtree(process_dir)
        # os.rmdir(process_dir)
        print("=====>DONE WRITTING FASTA FILES=======")
        print('Processed time (s): ', int(time.time() - t_start))

        del df_gff, fasta_dict, df_cds_1st, df_combined
        
    # ========= GFF VERSION 3 =============
    elif args.version == '3':
        print("=====LOADING GFF3 TABLE=======")
        df_gff, df_cds = load_gff3(file_path=args.gff, columns=COLUMNS_V3)
        print("=====> DONE LOADING GFF3 TABLE")
        df_gff.to_csv(f'{args.species}_loaded.gff', sep='\t', index=False)
        # df_cds.to_csv(f'{args.species}_cds.gff')

        print("=====LOADING DNA FASTA FILE=======")
        fasta_dict = load_fasta(fasta_path=args.fasta)
        # print(fasta_dict.keys())
        # with open(f"{args.species}_fasta_dict.json", "w", encoding="utf-8") as f:
        #     json.dump(fasta_dict, f, indent=4)   # indent=4 makes it pretty
        print("=====>DONE LOADING DNA FASTA FILE=======")

        # print(df)
        print("=====EXTRACTING CDS COORDINATES=======")
        df_cds_1st = get_coordinates(df_gff=df_cds, groupby_attr=GROUPBY_ATTR_V3)
        # df_cds_1st.to_csv(f'{args.species}_cds_1st.gff')
        print("=====>DONE EXTRACTING CDS COORDINATES=======")

        print("=====FINDING THE TRUTH M POSITION=======")
        # df_combined = find_M(df=df_coordinates, fasta_dict=fasta_dict)  
        df_combined = findM_parallel(df=df_cds_1st, fasta_dict=fasta_dict)
        df_combined.to_csv(f'{args.species}_combined.gff', sep='\t', index=False)  
        print("=====>DONE FINDING THE TRUTH M POSITION=======")

        print("=====EXPORTING TO CSV=======")
        df_corrected = df_combined.copy()
        df_corrected.loc[:, 'start'] = df_corrected['start'].astype(int) - 1
        df_corrected.loc[:, 'end'] = df_corrected['end'].astype(int)
        # export to new gff
        df_corrected.to_csv(os.path.join(args.outdir, f'{args.species}_all_gene_corrected.gff'), sep='\t', index=False)
        print("=====>DONE EXPORT TO CSV=======")

        print("=====WRITTING FASTA FILES=======")
        # temporal directory for processing
        process_dir = os.path.join(args.outdir, 'processed_fasta')
        os.makedirs(process_dir, exist_ok=True)
        # writting out the fasta files
        writting_parallel(
            ver=str(args.version),
            species=args.species,
            out_dir=args.outdir,
            processed_dir=process_dir, 
            df_combined=df_combined, 
            df_gff=df_gff,
            fasta_dict=fasta_dict, 
            outtype="prot", 
            text_width=TEXT_WIDTH,
            groupby_col=GROUPBY_ATTR_V3)
        # delete processed redundancies
        shutil.rmtree(process_dir)
        # os.rmdir(process_dir)
        print("=====>DONE WRITTING FASTA FILES=======")
        print('Processed time (s): ', int(time.time() - t_start))

        del df_gff, fasta_dict, df_cds_1st, df_combined
       
    # ========= GFF OTHER VERSIONS ============= 
    else:
        print("Sorry, the current tool does not support your gff file version!")