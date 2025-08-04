import os 
import pandas as pd
from Bio.Seq import Seq
import textwrap
import re
import multiprocessing as mp
import json
import glob
import shutil
from multiprocessing import Pool, cpu_count



def load_gff(file_path:str ='./Emihu1_all_genes.gff', 
             columns:list=['scaffold', 'jgi', 'type', 'start', 'end', 'direction1', 'direction2', 'phase', 'infos']
             ) -> pd.DataFrame:
    """
    Read table from .GFF file and load in to Pandas.DataFrame table
    """
    # check path
    if not os.path.exists(file_path):
        print("Non existed input path, Please re-check!")
        return None
    
    # Load the GFF file
    df_gff = pd.read_csv(file_path, sep='\t', comment='#', names=columns)

    # Split the 'attributes' column by ';' into multiple columns
    df_gff_split = df_gff['infos'].str.split(';', expand=True)

    # Rename the new columns (optional)
    df_gff_split.columns = ['name', 'proteinId', 'n_exon']
    df_gff.columns = df_gff.columns.str.strip()

    # Add the split columns back to the original DataFrame
    df_gff = pd.concat([df_gff, df_gff_split], axis=1)
    df_gff.drop(['infos'], axis=1, inplace=True)
    df_gff = df_gff[~df_gff['proteinId'].isna()]

    # convert CDS position to integer and subtract by 1 because string start from 0 in python
    df_gff['start'] = df_gff['start'].astype(int) - 1
    
    return df_gff


def get_coordinates(
        df_gff:pd.DataFrame, 
        groupby_attr:str ='name'
        ) -> pd.DataFrame:
    """
    Filtering the first CDS coordinate infomation for each attribute ('Name' by default) group
    from the GFF dataframe 
    """
    # get the first cds of protein which is translated by 5-3
    result_53 = df_gff.groupby('name').apply(
        lambda group: group.loc[
            (group['type'] == 'CDS') & (group['direction2'] == '+')
        ].head(1)  # get first match
    ).reset_index(drop=True)

    # get the last cds of protein which is translated by 3-5
    result_35 = df_gff.groupby('name').apply(
        lambda group: group.loc[
            (group['type'] == 'CDS') & (group['direction2'] == '-')
        ].iloc[-1] if any((group['type'] == 'CDS') & (group['direction2'] == '-')) else None
    ).dropna().reset_index(drop=True)

    return pd.concat((result_53, result_35))



def load_fasta(fasta_path:str = './Emihu1_masked_scaffolds.fasta'):
    """
    Load fasta sequence from file .fasta into a dictionary with (key: scaffold, value: dna sequence)
    """
    if not os.path.exists(fasta_path):
        print("Non existed input path, Please re-check!")
        return None
    
    ret_dict = {}
    key = ''
    with open(fasta_path, 'r') as f:
        for line in f.readlines():
            if ">" in line:
                # translate acide amine sequence to protein before advance
                if key != '':
                    ret_dict[key] = ret_dict[key]
                # save scaffold number as dict key
                key = line[1:].strip()  # avoid '>'
                # create a value for storing the sequence
                ret_dict[key] = ''
            else:
                ret_dict[key] += line.strip()
    f.close()

    return ret_dict


def find_firstM(seq):
    posM = -1
    for i in range(len(seq)):
        if seq[i] == '*':
            return posM
        if seq[i] == 'M':
            posM = i
    return posM


def re_verifyM(
        fasta_seq: str, 
        cds_start: float,
        cds_end:float,
        direction:str,
        ) -> int:
    """
    Take an input of fasta acide amine sequence, CDS position from JGI 
    and re-find the M position in the given sequence
    """
    cds_start = int(cds_start)
    cds_end = int(cds_end)
    if cds_start < 0 or cds_start >= len(fasta_seq) or cds_start > cds_end or cds_end < 0 or cds_end >= len(fasta_seq):
        return -1

    # Forward
    if direction == "+":
        # translate DNA to AA
        # include the first translated aa
        sub_seq = fasta_seq[:cds_start+3]
        # identify the translation framework: +0,+1 or +2 depend on the length of sequence
        framework = len(sub_seq)%3
        sub_seq = str(Seq(sub_seq[framework:]).translate())
        # find the closest codon stop
        # consider the closest * near the cds start position so inverse make the search run faster
        inversed_seq = sub_seq[::-1]
        M_pos = find_firstM(inversed_seq)
        if M_pos >= 0:
            M_pos = (len(inversed_seq) - M_pos - 1)*3 + framework
        return M_pos
    # Backward
    else:
        # find the true (closest) M position by going backward the CDS sequence
        sub_seq = fasta_seq[cds_start : cds_end]
        offset = len(sub_seq)%3
        if offset > 0:
            sub_seq = str(Seq(sub_seq[offset:]).reverse_complement().translate())
        else:
            sub_seq = str(Seq(sub_seq).reverse_complement().translate())
        # find the closest M encoutered from the end of dna sequence
        M_pos = sub_seq.find('M')
        if M_pos >= 0:
            return cds_end - (M_pos)*3
        return int(M_pos)
    
    
def re_verifyM_v2(
        fasta_subseq: str, 
        cds_start: float,
        cds_end:float,
        direction:str,
        ) -> int:
    """
    Take an input of fasta acide amine sequence, CDS position from JGI 
    and re-find the M position in the given sequence
    """
    
    cds_start = int(cds_start)
    cds_end = int(cds_end)
    if cds_start < 0 or cds_start >= len(fasta_subseq) or cds_start > cds_end or cds_end < 0 or cds_end >= len(fasta_subseq):
        return -1
    
    if direction == '+':
        sub_seq = fasta_subseq[:cds_start+3]
        offset = len(sub_seq) % 3
        sub_seq = str(Seq(sub_seq[offset:]).translate()) if offset > 0 else str(Seq(sub_seq).translate())
        inversed_seq = sub_seq[::-1]    
        M_pos = find_firstM(inversed_seq)
        if M_pos >= 0:
            M_pos = (len(inversed_seq) - M_pos - 1)*3 + offset
    else:
        sub_seq = fasta_subseq[cds_end-3:]
        offset = len(sub_seq) % 3
        sub_seq = str(Seq(sub_seq[:-offset]).reverse_complement().translate()) if offset > 0 else str(Seq(sub_seq).reverse_complement().translate())
        inversed_seq = sub_seq[::-1]
        M_pos = find_firstM(inversed_seq)
        if M_pos >= 0:
            M_pos = cds_end + M_pos*3

    return M_pos


def find_M(
    df:pd.DataFrame,
    fasta_dict:dict
    ) -> pd.DataFrame:

    results = []
    for scaffold, group_df in df.groupby(['scaffold']):
        fasta_seq = fasta_dict[scaffold[0]]  # get the value
        starts = group_df['start'].astype(int).tolist()
        ends = group_df['end'].astype(int).tolist()
        directions = group_df['direction2'].tolist()

        # Vectorized by list comprehension per group
        true_M_positions = [re_verifyM(fasta_seq, cds_start, cds_end, direction) for (cds_start, cds_end, direction) in zip(starts, ends, directions)]

        group_result = group_df.copy()
        group_result['true_M_pos'] = true_M_positions
        group_result['no_M'] = False
        group_result.loc[group_result['true_M_pos'] == -1, 'no_M'] = True 

        # calculate the differences
        group_result['M_pos_diff'] = 0
        # only gene in which has M can calculate the difference
        mask = ~group_result['no_M']
        
        # forward because of direction +
        mask2 = group_result['direction2'] == '+'
        group_result.loc[mask & mask2, 'M_pos_diff'] = (
            group_result.loc[mask & mask2, 'start'] - group_result.loc[mask & mask2, 'true_M_pos']
        ).abs()

        # backward because of direction -
        mask3 = group_result['direction2'] == '-'
        group_result.loc[mask & mask3, 'M_pos_diff'] = (
            group_result.loc[mask & mask3, 'end'] - group_result.loc[mask & mask3, 'true_M_pos']
        ).abs()
        results.append(group_result)

    # Concatenate all groups
    return pd.concat(results, ignore_index=True)


def get_fastas(df, fasta_dict, out_dir, species, text_width, output_type):
    """
    Writting the new fasta files based on the coordinates corrected
    """
    os.makedirs(out_dir, exist_ok=True)
    unch_dna_path = os.path.join(out_dir, f'{species}_Unchanged.fasta') 
    ch_dna_path = os.path.join(out_dir, f'{species}_Changed.fasta') 
    err_dna_path = os.path.join(out_dir, f'{species}_Error.fasta') 

    unch_prot_path = os.path.join(out_dir, f'{species}_Unchanged_prot.fasta') 
    ch_prot_path = os.path.join(out_dir, f'{species}_Changed_prot.fasta') 
    err_prot_path = os.path.join(out_dir, f'{species}_Error_prot.fasta') 

    def writting_fasta(scaffold, df_cleaned, df_gff, path_dna, path_prot, is_backward=False):
        #  boolean flag to remove the empty file
        is_empty_prot = True
        is_empty_dna = True
        prev_name = ''
        count = 1
        with open(path_dna, 'a') as f_dna, open(path_prot, 'a') as f_prot:
            for idx1, row in df_cleaned.iterrows():
                if prev_name != row['name']:
                    count = 1
                else:
                    count += 1
                prev_name = row['name']
                # print(row)
                # get the full CDS data of that model name
                df_cds = df_gff.loc[(df_gff['name'] == row['name']) & (df_gff['type'] == 'CDS')]
                df_cds = df_cds.sort_values(by='start', ascending=True)
                # forward
                if row['direction2'] == "+":
                    # the first CDS that we modified
                    seq_dna = fasta_dict[scaffold][int(row['start']-1):int(row['end'])]
                    # get the rest part
                    for idx2, row2 in df_cds.iloc[1:].iterrows():
                        seq_dna += fasta_dict[scaffold][row2['start']-1:row2['end']]
                # backward
                else:
                    seq_dna = ''
                    # get the first part of sequence
                    for idx2, row2 in df_cds.iloc[:-1].iterrows():
                        seq_dna += fasta_dict[scaffold][row2['start']-1:row2['end']]
                    # the last part that we modified
                    seq_dna = fasta_dict[scaffold][int(row['start']-1):int(row['end'])]
                
                # writting with output conditions
                if output_type == 'dna':
                    if count == 1: # only print at the first time
                        f_dna.write(f">jgi | {scaffold} | {row['proteinId']} | {row['name']}\n")
                    is_empty_dna = False
                    wraps_dna = textwrap.wrap(seq_dna, width=text_width)
                    for w in wraps_dna:
                        f_dna.write(w + "\n")
                # if user want prot fasta file or bth prot and dna fasta files
                else:
                    if count == 1: # only print at the first time
                        f_prot.write(f">jgi | {scaffold} | {row['proteinId']} | {row['name']}\n")
                    is_empty_prot = False
                    offset = len(seq_dna)%3
                    # forward
                    if row['direction2'] == "+":
                        if offset > 0:
                            seq_prot = str(Seq(seq_dna[:-offset]).translate())
                        else:
                            seq_prot = str(Seq(seq_dna).translate())
                    # backward
                    else:
                        if offset > 0:
                            seq_prot = str(Seq(seq_dna[offset:]).reverse_complement().translate())
                        else:
                            seq_prot = str(Seq(seq_dna).reverse_complement().translate())
                    wraps_prot = textwrap.wrap(seq_prot, width=text_width)
                    for w in wraps_prot:
                        f_prot.write(w + "\n")
                    if output_type == 'both':
                        f_dna.write(f">jgi | {scaffold} | {row['proteinId']} | {row['name']}\n")
                        is_empty_dna = False
                        wraps_dna = textwrap.wrap(seq_dna, width=text_width)
                        for w in wraps_dna:
                            f_dna.write(w + "\n")

        f_dna.close()
        f_prot.close()
        # delete the empty content text files
        if is_empty_dna:
            os.remove(path_dna)
        if is_empty_prot:
            os.remove(path_prot)

    for scaffold, group_df in df.groupby(['scaffold']):
        # writting the unchanged genes
        scaffold = scaffold[0]
        df_unchanged = group_df[(group_df['M_pos_diff'] == 0) & (~group_df['no_M'])]
        writting_fasta(scaffold, df_unchanged, unch_dna_path, unch_prot_path)

        # writting the changed genes
        df_changed = group_df[group_df['M_pos_diff'] != 0]
        # forward
        df_changed.loc[df_changed['direction2'] == "+", 'start'] = df_changed.loc[df_changed['direction2'] == "+", 'true_M_pos']
        # backward
        df_changed.loc[df_changed['direction2'] == "-", 'end'] = df_changed.loc[df_changed['direction2'] == "-", 'true_M_pos']
        writting_fasta(scaffold, df_changed, ch_dna_path, ch_prot_path)

        # writting the error genes
        df_error = group_df[group_df['no_M']]
        writting_fasta(scaffold, df_error, err_dna_path, err_prot_path)



### Multiprocessing

def process_one_scaffold(args):
    scaffold, group_df, fasta_dict = args

    fasta_seq = fasta_dict[scaffold]  # if scaffold is a tuple, e.g. ('scaffold_1',)
    starts = group_df['start'].astype(int).tolist()
    ends = group_df['end'].astype(int).tolist()
    directions = group_df['direction2'].tolist()

    true_M_positions = [
        re_verifyM_v2(fasta_seq, s, e, d)
        for s, e, d in zip(starts, ends, directions)
    ]

    group_result = group_df.copy()
    group_result['true_M_pos'] = true_M_positions
    group_result['no_M'] = group_result['true_M_pos'] == -1

    group_result['M_pos_diff'] = 0
    mask = ~group_result['no_M']

    # + direction
    mask2 = group_result['direction2'] == '+'
    group_result.loc[mask & mask2, 'M_pos_diff'] = (
        group_result.loc[mask & mask2, 'start'] - group_result.loc[mask & mask2, 'true_M_pos']
    ).abs()

    # - direction
    mask3 = group_result['direction2'] == '-'
    group_result.loc[mask & mask3, 'M_pos_diff'] = (
        group_result.loc[mask & mask3, 'end'] - group_result.loc[mask & mask3, 'true_M_pos']
    ).abs()

    return group_result


def findM_parallel(df, fasta_dict):
    args = [(scaffold, group.copy(), fasta_dict)
            for scaffold, group in df.groupby('scaffold')]

    with Pool(cpu_count()) as pool:
        results = pool.map(process_one_scaffold, args)

    return pd.concat(results, ignore_index=True)



def writting_fasta(scaffold, df_combined, df_gff, fasta_dict, path_dna, path_prot, output_type, text_width):
    prev_name = ''
    count = 1
    with open(path_dna, 'a') as f_dna, open(path_prot, 'a') as f_prot:
        for idx1, row in df_combined.iterrows():
            if prev_name != row['name']:
                count = 1
            else:
                count += 1
            prev_name = row['name']
            # print(row)
            # get the full CDS data of that model name
            df_cds = df_gff.loc[(df_gff['name'] == row['name']) & (df_gff['type'] == 'CDS')]
            df_cds = df_cds.sort_values(by='start', ascending=True)

            # forward
            if row['direction2'] == "+":
                # the first CDS that we modified in df_combined
                seq_dna = fasta_dict[scaffold][int(row['start']):int(row['end'])]
                # get the rest part from df_gff
                for idx2, row2 in df_cds.iloc[1:].iterrows():
                    seq_dna += fasta_dict[scaffold][row2['start']:row2['end']]
            
            # backward
            else:
                seq_dna = ''
                # get the first part of sequence
                for idx2, row2 in df_cds.iloc[:-1].iterrows():
                    seq_dna += fasta_dict[scaffold][row2['start']:row2['end']]
                # the last part that we modified
                seq_dna = fasta_dict[scaffold][int(row['start']):int(row['end'])]
            
            # writting with output conditions
            if output_type == 'dna':
                if count == 1: # only print at the first time
                    f_dna.write(f">jgi | {scaffold} | {row['proteinId']} | {row['name']}\n")
                # is_empty_dna = False
                wraps_dna = textwrap.wrap(seq_dna, width=text_width)
                for w in wraps_dna:
                    f_dna.write(w + "\n")
            # if user want prot fasta file or bth prot and dna fasta files
            else:
                if count == 1: # only print at the first time
                    f_prot.write(f">jgi | {scaffold} | {row['proteinId']} | {row['name']}\n")
                # is_empty_prot = False
                offset = len(seq_dna)%3
                # forward
                if row['direction2'] == "+":
                    if offset > 0:
                        seq_prot = str(Seq(seq_dna[:-offset]).translate())
                    else:
                        seq_prot = str(Seq(seq_dna).translate())
                # backward
                else:
                    if offset > 0:
                        seq_prot = str(Seq(seq_dna[offset:]).reverse_complement().translate())
                    else:
                        seq_prot = str(Seq(seq_dna).reverse_complement().translate())
                    # seq_prot = str(Seq(seq_dna).reverse_complement().translate())
                wraps_prot = textwrap.wrap(seq_prot, width=text_width)
                for w in wraps_prot:
                    f_prot.write(w + "\n")
                if output_type == 'both':
                    f_dna.write(f">jgi | {scaffold} | {row['proteinId']} | {row['name']}\n")
                    # is_empty_dna = False
                    wraps_dna = textwrap.wrap(seq_dna, width=text_width)
                    for w in wraps_dna:
                        f_dna.write(w + "\n")
    f_dna.close()
    f_prot.close()


def process_scaffold_v2(args):
    scaffold, \
    group_df, \
    df_gff, \
    fasta_dict, \
    process_dir, \
    outtype, \
    textwidth = args

    scaffold_id = scaffold
    # temporal fasta files for each scaffold
    unch_dna_path = os.path.join(process_dir, f"{scaffold}_unchanged_dna.fasta")
    ch_dna_path   = os.path.join(process_dir, f"{scaffold}_changed_dna.fasta")
    err_dna_path  = os.path.join(process_dir, f"{scaffold}_error_dna.fasta")

    unch_prot_path = os.path.join(process_dir, f"{scaffold}_unchanged_prot.fasta")
    ch_prot_path   = os.path.join(process_dir, f"{scaffold}_changed_prot.fasta")
    err_prot_path  = os.path.join(process_dir, f"{scaffold}_error_prot.fasta")


    df_unchanged = group_df[(group_df['M_pos_diff'] == 0) & (~group_df['no_M'])]
    if not df_unchanged.empty:
        writting_fasta(scaffold_id, df_unchanged, df_gff, fasta_dict, unch_dna_path, unch_prot_path, outtype, textwidth)

    df_changed = group_df[group_df['M_pos_diff'] != 0].copy()
    if not df_changed.empty:
        df_changed.loc[df_changed['direction2'] == "+", 'start'] = df_changed.loc[df_changed['direction2'] == "+", 'true_M_pos']
        df_changed.loc[df_changed['direction2'] == "-", 'end']   = df_changed.loc[df_changed['direction2'] == "-", 'true_M_pos']
        writting_fasta(scaffold_id, df_changed, df_gff, fasta_dict, ch_dna_path, ch_prot_path, outtype, textwidth)

    df_error = group_df[group_df['no_M']]
    if not df_error.empty:
        writting_fasta(scaffold_id, df_error, df_gff, fasta_dict, err_dna_path, err_prot_path, outtype, textwidth)


def merge_fasta(output_file, pattern):
    with open(output_file, "wb") as outfile:
        for fname in sorted(glob.glob(pattern)):
            with open(fname, "rb") as infile:
                shutil.copyfileobj(infile, outfile, length=1024 * 1024 * 10)  # 10MB buffer


def writting_parallel(species, out_dir, processed_dir, df_combined, df_gff, fasta_dict, outtype, text_width):
    # df_combined = pd.read_csv("df_combined.csv") 
    # df_gff = pd.read_csv("df_gff.csv")    

    # with open("fasta_dict.json", "r") as f:
    #     fasta_dict = json.load(f)

    args = [(scaffold,\
            group.copy(),\
            df_gff, \
            fasta_dict, \
            processed_dir,\
            outtype, \
            text_width) for scaffold, group in df_combined.groupby('scaffold')]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.map(process_scaffold_v2, args)

    unch_dna_path = os.path.join(out_dir, f'{species}_Unchanged_dna.fasta')
    ch_dna_path   = os.path.join(out_dir, f'{species}_Changed_dna.fasta')
    err_dna_path  = os.path.join(out_dir, f'{species}_Error_dna.fasta')
    unch_prot_path = os.path.join(out_dir, f'{species}_Unchanged_prot.fasta')
    ch_prot_path   = os.path.join(out_dir, f'{species}_Changed_prot.fasta')
    err_prot_path  = os.path.join(out_dir, f'{species}_Error_prot.fasta')

    merge_fasta(unch_dna_path, os.path.join(processed_dir, "*_unchanged_dna.fasta"))
    merge_fasta(ch_dna_path, os.path.join(processed_dir,"*_changed_dna.fasta"))
    merge_fasta(err_dna_path, os.path.join(processed_dir,"*_error_dna.fasta"))

    merge_fasta(unch_prot_path, os.path.join(processed_dir,"*_unchanged_prot.fasta"))
    merge_fasta(ch_prot_path, os.path.join(processed_dir,"*_changed_prot.fasta"))
    merge_fasta(err_prot_path, os.path.join(processed_dir,"*_error_prot.fasta"))