import os 
import pandas as pd
from Bio.Seq import Seq
import textwrap
import re


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
    
    # Group by 'group' and find the first row with a non-null value
    result = df_gff.groupby(groupby_attr, group_keys=False).apply(
        lambda group: group.loc[group[group['type'] == 'CDS'].index[0]] if any(group['type'] == 'CDS') else None,
        include_groups=False
    )

    return result


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


def re_verifyM(
        fasta_seq: str, 
        cds_pos: float,
        direction:str,
        ) -> int:
    """
    Take an input of fasta acide amine sequence, CDS position from JGI 
    and re-find the M position in the given sequence
    """
    cds_pos = int(cds_pos)
    if cds_pos < 0 or cds_pos >= len(fasta_seq):
        return cds_pos

    # include the first translated aa
    sub_seq = fasta_seq[:cds_pos+3]
    # identify the translation framework: +0,+1 or +2 depend on the length of sequence
    framework = len(sub_seq)%3
    # translate DNA to AA
    if direction == "+":
        sub_seq = str(Seq(sub_seq[framework:]).translate())
    else:
        sub_seq = str(Seq(sub_seq[framework:]).reverse_complement().translate())
    
    # find the closest codon stop
    # consider the closest * near the cds start position so inverse make the search run faster
    inversed_seq = sub_seq[::-1]
    match = re.search(r'\*', inversed_seq)
    # the first codon stop encoutered
    stop_pos = match.start() if match else len(sub_seq)

    # find the closest aa M to the first * by slicing again the inversed aa sequence
    tmp_seq = inversed_seq[:stop_pos]
    # find the first M encoutered
    # if no M, M_pos = -1
    M_pos = tmp_seq.find('M')
    if M_pos >= 0:
        return (len(inversed_seq) - M_pos - 1)*3 + framework

    return int(M_pos)


def find_M(
    df:pd.DataFrame,
    fasta_dict:dict
    ) -> pd.DataFrame:

    results = []
    for scaffold, group_df in df.groupby(['scaffold']):
        fasta_seq = fasta_dict[scaffold[0]]  # get once
        starts = group_df['start'].astype(int).tolist()
        directions = group_df['direction2'].tolist()

        # Vectorized by list comprehension per group
        true_M_positions = [re_verifyM(fasta_seq, cds_pos, direction) for (cds_pos, direction) in zip(starts, directions)]

        group_result = group_df.copy()
        group_result['true_M_pos'] = true_M_positions
        group_result['no_M'] = False
        group_result.loc[group_result['true_M_pos'] < 0, 'no_M'] = True 

        # calculate the differences
        group_result['M_pos_diff'] = 0
        # only gene in which has M can calculate the difference
        mask = ~group_result['no_M']
        group_result.loc[mask, 'M_pos_diff'] = (
            group_result.loc[mask, 'start'] - group_result.loc[mask, 'true_M_pos']
        ).abs()
        results.append(group_result)

    # Concatenate all groups
    result = pd.concat(results, ignore_index=True)


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

    def writting_fasta(scaffold, df, path_dna, path_prot):
        #  boolean flag to remove the empty file
        is_empty_prot = True
        is_empty_dna = True
        with open(path_dna, 'a') as f_dna, open(path_prot, 'a') as f_prot:
            f_dna.write(f'>{scaffold}\n')
            f_prot.write(f'>{scaffold}\n')

            for row in df.iterrows():
                seq_dna = fasta_dict[scaffold][row[start]:row[end]]
                
                if output_type == 'dna':
                    f_dna.write(f">jgi|{row['proteinId']} | {row['name']}\n")
                    is_empty_dna = False
                    wraps_dna = textwrap.wrap(seq_dna, width=text_width)
                    for w in wraps:
                        f_dna.write(w + "\n")
                # if user want prot fasta file or bth prot and dna fasta files
                else:
                    f_prot.write(f">jgi|{row['proteinId']} | {row['name']}\n")
                    is_empty_prot = False
                    seq_prot = str(Seq(seq_dna).translate())
                    wraps_prot = textwrap.wrap(seq_prot, width=text_width)
                    for w in wraps_prot:
                        f_prot.write(w + "\n")

                    if output_type == 'both':
                        f_dna.write(f">jgi|{row['proteinId']} | {row['name']}\n")
                        is_empty_dna = False
                        wraps_dna = textwrap.wrap(seq_dna, width=text_width)
                        for w in wraps:
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
        df_changed = group_df[group_df['M_pos_diff'] > 0]
        writting_fasta(scaffold, df_changed, ch_dna_path, ch_prot_path)

        # writting the unchanged genes
        df_error = group_df[group_df['no_M']]
        writting_fasta(scaffold, df_error, err_dna_path, err_prot_path)
