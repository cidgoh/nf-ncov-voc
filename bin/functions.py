import pandas as pd
import numpy as np
import logging

# standard variables used by all scripts
empty_attributes = 'ID=;original_mutation_description=;alias=;gene=;gene_name=;gene_symbol=; \
    strand_orientation=;gene_orientation=; \
    product=;protein_alias=;protein_name=;protein_symbol=;\
    protein_id=;alias_protein_id=;locus_tag=;ps_filter=;ps_exc=; \
    mat_pep=;mat_pep_desc=;mat_pep_acc=;ro=;ao=;dp=;sample_size=; \
    Reference_seq=;Variant_seq=;nt_name=;aa_name=;hgvs_nt=;hgvs_aa=;hgvs_alias=;heterozygosity=; \
    vcf_gene=;mutation_type=;viral_lineage=;clade=;multi_aa_name=; \
    multiaa_comb_mutation=;alternate_frequency=;measured_variant_functional_effect=;inferred_variant_functional_effect=; \
    viral_life_cycle_functional_effect=;author=;publication_year=;DOI=;PMID=;URL=;peer_review_status=; \
    citation=;comb_mutation=;variant_functional_effect_description=; \
    mutation_functional_annotation_resource=;functional_annotation_resource=;curator=; \
    clade_defining=;variant=;variant_type=;voi_designation_date=; \
    voc_designation_date=;vum_designation_date=;status=; \
    organism=;reference_accession=;reference_database_name=; \
    CVX_code=;DrugBank_Accession_Number=;Antibody_Registry_ID=;'

empty_attributes = empty_attributes.replace(" ", "")

gvf_columns = ['#seqid', '#source', '#type', '#start', '#end', '#score', '#strand', '#phase', '#attributes']

vcf_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                'FILTER', 'INFO', 'FORMAT', 'unknown']

# pragmas are in column 0
pragmas = pd.DataFrame([['##gff-version 3'],
                        ['##gvf-version 1.10'],
                        ['##sample-description'],
                        ['##species']])

# source: https://hgvs-nomenclature.org/stable/background/standards/#amino-acid-descriptions
amino_acid_codes_dict = {'*': 'Ter', 'A': 'Ala', 'B': 'Asx', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'U': 'Sec', 'V': 'Val', 'W': 'Trp', 'X': 'Xaa', 'Y': 'Tyr', 'Z': 'Glx'}

def convert_amino_acid_codes(one_letter_mutation_name):
    three_letter_code_result = "".join([
        amino_acid_codes_dict[ch]
        if (ch.isupper() and ch in amino_acid_codes_dict.keys()) 
        else ch for ch in one_letter_mutation_name])

    return(three_letter_code_result)

def rewrite_nt_snps_as_hgvs(original_nt_name):
    ## NOTE: this function doesn't add the ref sequence; this must be done separately

    pos = "".join([i for i in original_nt_name if i.isdigit()])
    nts_list = [ch for ch in original_nt_name if ch.isupper()]
    hgvs_name = "g." + pos + nts_list[0] + ">" + nts_list[1] ##may change "g." in future

    return(hgvs_name)

def remove_nts_from_nt_name(original_nt_name):
    # sufficient for nt names containing "del" or "dup"
    # will not work for "ins" or "delins"
    hgvs_name = "".join([ch for ch in original_nt_name if not ch.isupper()])

    return(hgvs_name)

'''
def convert_nt_delins_to_hgvs():
    delins_regex = "[a-z]\\.[A-Z]+[0-9]+_[0-9]+[A-Z]+"
    # eg. input: "g.TATT430_433ACTTCTA"
    # eg. output: "g.430_433delinsACTTCTA"

    return(hgvs_name)
'''

def add_hgvs_names(new_gvf):
    
    new_gvf['hgvs_nt'] = 'n/a'
    new_gvf['hgvs_aa'] = 'n/a'
    new_gvf['hgvs_alias'] = 'n/a'
    
    # fill in 'hgvs_nt'
    
    # define nt regex patterns
    # for SNPs, eg. g.C45T, g.C-45T
    nt_snp_regex = "[a-z]\\.[A-Z][0-9\\-]+[A-Z]"
    # for dels and dups, eg. g.254_259delTGGTTG, g.361delA
    nt_del_dup_regex = "[a-z]\\.[0-9\\-_]+(?:del|dup)[A-Z]+"
    # for ins
    nt_ins_regex = "[a-z]\\.[0-9\\-_]+ins[A-Z]+"
    # for delins, eg. g.GCC10182_10184ACA
    nt_delins_regex = "[a-z]\.[A-Z]+-?[0-9]+_-?[0-9]+[A-Z]+" 
    
    #df.loc[mask, 'val'] = df.loc[mask, 'val'].apply(f)
    # add hgvs nt snp names
    nt_snp_mask = new_gvf['nt_name'].str.contains(nt_snp_regex, regex=True)
    new_gvf.loc[nt_snp_mask, 'hgvs_nt'] = new_gvf['#seqid'] + "(" + new_gvf['locus_tag'] + "):" + \
                 new_gvf.loc[nt_snp_mask, 'nt_name'].apply(rewrite_nt_snps_as_hgvs)
    # add hgvs nt dels and dups
    nt_del_dup_mask = new_gvf['nt_name'].str.contains(nt_del_dup_regex, regex=True)
    new_gvf.loc[nt_del_dup_mask,'hgvs_nt'] = new_gvf['#seqid'] + "(" + new_gvf['locus_tag'] + "):" + \
                new_gvf.loc[nt_del_dup_mask, 'nt_name'].apply(remove_nts_from_nt_name)
    # add hgvs nt ins
    nt_ins_mask = new_gvf['nt_name'].str.contains(nt_ins_regex, regex=True)
    new_gvf.loc[nt_ins_mask, 'hgvs_nt'] = new_gvf['#seqid'] + "(" + new_gvf['locus_tag'] + "):" + new_gvf['nt_name']  
    # add hgvs nt delins: change to eg. g.123_129delinsAC
    nt_delins_mask = new_gvf['nt_name'].str.contains(nt_delins_regex, regex=True)
    #new_gvf.loc[nt_delins_mask, 'hgvs_nt'] = new_gvf['#seqid'] + "(" + new_gvf['locus_tag'] + "):" + "TBA!" #new_gvf['nt_name']  

    # remove transcript id in parentheses for all HGVS nucleotide names where 'locus_tag'=='n/a'
    new_gvf['hgvs_nt'] = new_gvf['hgvs_nt'].str.replace("(n/a)", "", regex=False)

    # define aa regex patterns
    aa_snp_regex = "[A-Z*][0-9\\-]+[A-Z*]"
    aa_other_regex = "[A-Z*]+[0-9\\-]+(?:del|delins|ins|dup|fs|ext)[A-Z*]*"
    
    # fill in 'hgvs_aa' 
    # add hgvs aa snps to rows with protein_id!=n/a
    aa_snp_mask = (new_gvf['aa_name'].str.contains(aa_snp_regex, regex=True)) & (new_gvf['protein_id']!='n/a')
    new_gvf.loc[aa_snp_mask, 'hgvs_aa'] = \
                new_gvf["protein_id"] + ":" + \
                new_gvf.loc[aa_snp_mask, 'aa_name'].apply(convert_amino_acid_codes)
    # add hgvs aa names for non-snps
    aa_non_snp_mask = (new_gvf['aa_name'].str.contains(aa_other_regex, regex=True)) & (new_gvf['protein_id']!='n/a')
    new_gvf.loc[aa_non_snp_mask, 'hgvs_aa'] = \
                new_gvf["protein_id"] + ":" + new_gvf.loc[aa_non_snp_mask, 'aa_name'].apply(
                    convert_amino_acid_codes)

    # fill in 'hgvs_alias'
    # add hgvs alias snps to rows with alias_protein_id!=n/a
    alias_snp_mask = (new_gvf['alias'].str.contains(aa_snp_regex, regex=True)) & (new_gvf['alias_protein_id']!='n/a') & (new_gvf['alias']!='n/a')
    new_gvf.loc[alias_snp_mask, 'hgvs_alias'] = \
                new_gvf["alias_protein_id"] + ":" + new_gvf.loc[alias_snp_mask, 'alias'].apply(
                    convert_amino_acid_codes)
    # add hgvs alias names for non-snps
    alias_non_snp_mask = (new_gvf['alias'].str.contains(aa_other_regex, regex=True)) & (new_gvf['alias_protein_id']!='n/a') & (new_gvf['alias']!='n/a')
    new_gvf.loc[alias_non_snp_mask, 'hgvs_alias'] = \
                new_gvf["alias_protein_id"] + ":" + new_gvf.loc[alias_non_snp_mask, 'alias'].apply(
                    convert_amino_acid_codes)
                
    return(new_gvf)


def separate_attributes(df):
    # expand #attributes column into multiple columns for each attribute,
    # keeping the original #attributes column
    
    # split #attributes column into separate columns for each tag
    # split at ;, form dataframe
    attributes = df['#attributes'].str.split(pat=';').apply(pd.Series)
    # last column is a copy of the index so drop it
    attributes = attributes.drop(labels=len(attributes.columns) - 1,
                                 axis=1)

    for column in attributes.columns:
        split = attributes[column].str.split(pat='=').apply(pd.Series)
        title = split[0].drop_duplicates().tolist()[0] #.lower()

        content = split[1]

        # ignore "tag=" in column content
        attributes[column] = content
        # make attribute tag as column label
        attributes.rename(columns={column: title}, inplace=True)

    # replace attributes column in the original df with the new
    # separated out attributes
    df = pd.concat((df, attributes), axis=1)

    return(df)


def rejoin_attributes(df, empty_attributes_str):
    # get column names as list
    columns_to_join = empty_attributes_str.split('=;')[:-1] #last one will be empty
    for col in columns_to_join:
        df[col] = col + "=" + df[col].astype(str) + ';'
    # replace #attributes column with filled attributes
    df['#attributes'] = df[columns_to_join].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    # make sure all empty key pairs are written as "attribute_key=n/a"
    df['#attributes'] = df['#attributes'].astype(str).str.replace('=;', '=n/a;', regex=False)
    # drop extra columns
    df = df.drop(columns=columns_to_join)
    
    return(df)
    

def get_unknown_labels(df):
# determines variant calling source (eg. iVar) based on pragmas
# returns GVF-relevant names for last column ("unknown") of vcf
    source = df['#CHROM'][df['#CHROM'].str.contains("##source=")].values[0].split("=")[1].split()[0]
    if source=="freeBayes":
        columns = [x.lower() for x in ["GT","DP","AD","RO","QR","AO","QA","GL"]]
    elif source=="iVar":
        # iVar names are: ["GT","REF_DP","REF_RV","REF_QUAL","ALT_DP","ALT_RV","ALT_QUAL","ALT_FREQ"]
        # use "RO" instead of "REF_DP" to match GVF standard
        # use "AO" instead of "ALT_DP" to match GVF standard
        columns = [x.lower() for x in ["GT","RO","REF_RV","REF_QUAL","AO","ALT_RV","ALT_QUAL","ALT_FREQ"]]
        
    return columns
        

def parse_pango_lineages(strain, dataframe):
    '''
    Expands the pango_lineage column in the clades file into
    a nested list called who_lineages, where each inside
    list is one expanded row of pango_lineages.
    
    Returns the nested list, as well as the row index for
    the input strain (if found).
    '''
    # initialize variables
    who_lineages = []
    var_to_match = 'None'
    
    # for each who_variant, expand the lineage names in
    # pango_lineages and save these to a list
    for var in dataframe["pango_lineage"]:
        var_pango_lineages = []
        if "," in var:
            for temp in var.split(","):
                if not "[" in var:
                    var_pango_lineages.append(temp)
                    var_pango_lineages.append(temp + ".*")
                    if strain.startswith(temp):
                        var_to_match = var
                else:
                    parent=temp[0]
                    child=temp[2:-1].split("|")
                    for c in child:
                        var_pango_lineages.append(parent + str(c))
                        var_pango_lineages.append(parent+str(c)+".*")
                        if strain.startswith(parent + str(c)):
                            var_to_match = var
        else:
            var_pango_lineages.append(var)
            if strain.startswith(var):
                var_to_match = var    
        
        # append the list for each variant to the larger
        # list for all variants
        who_lineages.append(var_pango_lineages)
        
        # get the row index for information on the input strain
        if var_to_match != 'None':
            strain_index = dataframe.index[dataframe['pango_lineage'] == \
                                 var_to_match].tolist()[0]
        else:
            strain_index = 'n/a'
        
    return who_lineages, strain_index


class get_variant_info:

    def __init__(self, strain, clades):

        # retrieve row number that matches the input strain
        who_lineages, var_index = parse_pango_lineages(strain, clades)

        # save status, WHO strain name, etc. from clades file
        if var_index != 'n/a':
            self.who_variant = clades.loc[var_index, 'variant']
            self.variant_type = clades.loc[var_index, 'variant_type']
            self.voi_designation_date = clades.loc[var_index, 'voi_designation_date']
            self.voc_designation_date = clades.loc[var_index, 'voc_designation_date']
            self.vum_designation_date = clades.loc[var_index, 'vum_designation_date']
            self.status = clades.loc[var_index, 'status']
            self.strain_in_cladefile = True # flag
        else:
            self.strain_in_cladefile = False
        

def unnest_multi(df, columns, reset_index=False):
# expands out columns of lists into 1d, as well as
# duplicating other non-specified rows as needed.
# all the lists must be the same length across columns in a given row, but
# can vary between rows
# adapted from https://stackoverflow.com/questions/21160134/flatten-a-column-with-value-of-type-list-while-duplicating-the-other-columns-va
    '''
    df_flat = pd.DataFrame(columns=columns)
    print("df_flat")
    print(df_flat[df_flat.index.duplicated()])
    print("")
    for col in columns:
        col_flat = pd.DataFrame([[i, x] 
                       for i, y in df[col].apply(list).iteritems() 
                           for x in y], columns=['I', col])
        print("whole col flat")
        print(col_flat)
        #col_flat = col_flat.set_index('I')
        print("")
        print("col_flat duplications")
        print(col_flat[col_flat.index.duplicated()])
        df_flat[col] = col_flat
    df = df.drop(labels=columns, axis=1)
    df = df.merge(df_flat, left_index=True, right_index=True)
        '''
    df = df.explode(columns, ignore_index=False)
    if reset_index:
        df = df.reset_index(drop=True)
    return df


def select_snpeff_records(eff_string, ao_count):
    
    eff_list = eff_string.split(",")

    # if any records in the row contain '|p.', take only those records
    EFF_records_list = [s for s in eff_list if '|p.' or 'LOF' in s]

    # if no records contain '|p.', take the "intergenic" record
    if len(EFF_records_list) == 0:
        EFF_records_list = [s for s in eff_list if 'intergenic_region' in s]

    # filter out annotations that include 'WARNING' or 'GU280_gp01.2'
    # this keeps 'GU280_gp01' annotations over 'GU280_gp01.2'
    EFF_records_list = [s for s in EFF_records_list if 'WARNING' not in s
                        and 'GU280_gp01.2' not in s]
    
    # of the filtered records, take only the first N, where N is the number
    # of comma-separated AO values given in the "unknown" column
    EFF_records_list = EFF_records_list[:ao_count]
    
    return EFF_records_list


def find_sample_size(table, lineage, vcf_file, wastewater):
    sample_size='n/a'
    if table != 'n/a':
        strain_tsv_df = pd.read_csv(table, delim_whitespace=True,
                                    usecols=['file', 'num_seqs'])

        # Reference mode
        if lineage != 'n/a':
            num_seqs = strain_tsv_df[strain_tsv_df['file'].str.startswith(
                lineage)]['num_seqs'].values
            sample_size = num_seqs[0]

        # wastewater data or user-uploaded fasta
        else:
            # num_seqs values should be identical, so take the
            # first value in num_seqs to be sample_size
            sample_size = strain_tsv_df['num_seqs'].values[0]
            # if forward and backward reads have different
            # num_seqs values, log an error
            if strain_tsv_df.num_seqs.nunique()!=1:
                err = "Different values in 'num_seqs' in " + table
                logging.info(err)

    # user-uploaded vcf
    elif table == 'n/a' and lineage == 'n/a':
        sample_size = 'n/a'

    if type(sample_size) == str:
        sample_size = sample_size.replace(",","")
        
    return sample_size



    
def parse_INFO(df, var_cols): # return INFO dataframe with named columns, including EFF split apart

    # extract these key-value pairs in INFO into their own columns
    cols_to_extract = ['DP', 'ps_filter', 'ps_exc', 'EFF', 'mat_pep', 'mat_pep_desc', 'mat_pep_acc']
    info = pd.DataFrame(columns=cols_to_extract)
    for col in cols_to_extract:
        pat = str(col) + "\=(.*?)(?=[\;]|$)" 
        info[col] = df['INFO'].str.extract(pat)
    # rename uppercase columns as lowercase
    info = info.rename(columns={'DP':'dp', 'EFF':'eff'})

    # concatenate info and df horizontally
    df = pd.concat([df, info], axis=1)
    df = df.drop(columns="INFO")

    # expand "unknown" column into multiple named columns
    unknown = df['unknown'].str.split(pat=':').apply(pd.Series)
    unknown.columns = var_cols
    #drop columns in df that have the same name as 'unknown' column names
    cols_to_drop = list(set(df.columns) & set(unknown.columns)) 
    df = df.drop(columns=cols_to_drop)

    df = pd.concat([df, unknown], axis=1)
    # make ALT, AO, type into lists
    for column in ["ao", "ALT"]:
        df[column] = df[column].str.split(",")
    if "type" in df.columns: # "type" is not an attribute of INFO for wastewater
        df["type"] = df["type"].str.split(",")   
    # get number of AO values given in "unknown" column
    df['ao_count'] = df["ao"].str.len()
    
    # parse EFF entry from INFO
    df["eff_result"] = [select_snpeff_records(x, y) for x, y in
                        zip(df['eff'], df["ao_count"])]
    #df.to_csv("eff_result_checking.tsv", sep="\t")
    # check how many "type" entries there are
    #df['eff_result_len'] = df["eff_result"].str.len()
    #df['type_len'] = df["type"].str.len()
    #print(df.query('eff_result_len != type_len'))
    #mismatch = df.query(df.query('eff_result_len != type_len'))[['POS', 'eff_result', 'eff_result_len', 'ao', 'type']]
    #mismatch.to_csv("mismatches.csv", sep='\t', header=True, index=True)
    # unnest list columns
    if "type" in df.columns: # "type" is not an attribute of INFO for wastewater
        df = unnest_multi(df, ["eff_result", "ao", "ALT", "type"], reset_index=True)
    else:
        df = unnest_multi(df, ["eff_result", "ao", "ALT"], reset_index=True)
    # calculate Alternate Frequency
    df['AF'] = df['ao'].astype(int) / df['dp'].astype(int)

    # expand the contents of eff_result into separate columns, named as in the 
    # VCF header
    eff_info = df['eff_result'].str.findall('\\((.*?)\\)').str[0]
    # split at pipe, form dataframe
    eff_info = eff_info.str.split(pat='|').apply(pd.Series)
    num_cols = len(eff_info.columns)
    eff_info_cols = ['Effect_Impact','Functional_Class','Codon_Change','Amino_Acid_Change','Amino_Acid_length','Gene_Name','Transcript_BioType','Gene_Coding','Transcript_ID','Exon_Rank','Genotype ERRORS', 'Genotype WARNINGS'][:num_cols]
    eff_info.columns = eff_info_cols

    df = pd.concat([df, eff_info], axis=1)
    df = df.drop(columns='eff_result')
    

    # split df['Amino_Acid_Change'] into two columns: one for HGVS amino acid
    # names, and the righthand column for nucleotide-level names
    name_mask = df['Amino_Acid_Change'].str.contains('/')
    df.loc[~name_mask, 'Amino_Acid_Change'] = '/' + df['Amino_Acid_Change']
    hgvs = df['Amino_Acid_Change'].str.rsplit(pat='/').apply(pd.Series)
    hgvs.columns = ["hgvs_protein", "hgvs_nucleotide"]
    df = pd.concat([df, hgvs], axis=1)
    
    # make adjustments to the nucleotide names
    # 1) change 'c.' to 'g.'convert_amino_acid_codes for nucleotide names ### double check that we want this
    df["hgvs_nucleotide"] = df["hgvs_nucleotide"].str.replace("c.", "g.", regex=False) 
    df["hgvs_nucleotide"] = df["hgvs_nucleotide"].str.replace("n.", "g.", regex=False) 
    # 2) change nucleotide names of the form "g.C*4378A" to g.C4378AN;
    asterisk_mask = df["hgvs_nucleotide"].str.contains('*', regex=False)
    df.loc[asterisk_mask, "hgvs_nucleotide"] = 'g.' + df['REF'] + df['POS'] + \
        df['ALT']
    # 3) change 'Gene_Name' to "intergenic" where names contain a "*"
    df.loc[asterisk_mask, 'Gene_Name'] = "intergenic"

    # create a "Names" column that holds the amino acid name (minus 'p.')
    # if there is one, or the nucleotide level name if not
    df["Names"] = df["hgvs_nucleotide"]
    protein_mask = df["hgvs_protein"].str.contains("p.")
    df.loc[protein_mask, "Names"] = df["hgvs_protein"].str.replace(
        "p.", "", regex=False)

    # rename some columns
    df = df.rename(columns={'Gene_Name': "vcf_gene", 'Functional_Class':
                            "mutation_type", 'hgvs_nucleotide': 'nt_name',
                            'hgvs_protein': 'aa_name', 'REF': 'Reference_seq',
                            'ALT': 'Variant_seq'})

    return(df)
    

def split_names(names_to_split, new_gvf, col_to_split):
    # separate multi-aa names noted in names_to_split into separate rows
    ### MZA: This needs immediate attention with Paul and his group. Need to update the notion of mutations
    # split multi-aa names from the vcf into single-aa names (multi-row)
    # load names_to_split spreadsheet
    names_to_split_df = pd.read_csv(names_to_split, sep='\t', header=0)
    # remove spaces and quotation marks from 'split_into' column
    for x in [' ', "'"]:
        names_to_split_df['split_into'] = names_to_split_df[
            'split_into'].str.replace(x, "")
    # merge "split_into" column into new_gvf, matching up by Names
    names_to_split_df = names_to_split_df.rename(columns={'name': col_to_split})
    new_gvf = new_gvf.merge(names_to_split_df, on=col_to_split, how = 'left')
    # add "multi_aa_name" column containing the original multi-aa names
    new_gvf["multi_aa_name"] = ''
    new_gvf.loc[new_gvf["split_into"].notna(), "multi_aa_name"] = new_gvf[col_to_split]
    # where "split_into" is notna, replace "Names" value with "split_into" value
    new_gvf.loc[new_gvf["split_into"].notna(), col_to_split] = new_gvf["split_into"]
    # make 'Names' into a column of lists
    new_gvf[col_to_split] = new_gvf[col_to_split].str.split(",")
    # unnest these lists (convert to 1d)                  
    new_gvf = unnest_multi(new_gvf, columns=[col_to_split], reset_index=True)   
    # 'multiaa_comb_mutation' attribute is "split_into" column left over from merge
    new_gvf['multiaa_comb_mutation'] = new_gvf['split_into']
    new_gvf = new_gvf.drop(columns=['split_into'])
    # make multiaa_comb_mutation contain everything in split_into
    # except for the name in Names
    new_gvf.multiaa_comb_mutation = new_gvf.multiaa_comb_mutation.fillna('')
    new_gvf["multiaa_comb_mutation"] = new_gvf.apply(lambda row : row["multiaa_comb_mutation"].replace(row[col_to_split], ''), axis=1)
    # strip extra commas
    new_gvf["multiaa_comb_mutation"] = new_gvf["multiaa_comb_mutation"].str.strip(',').str.replace(',,',',')

    return(new_gvf)


def map_pos_to_gene_protein(pos, GENE_PROTEIN_POSITIONS_DICT):
    """This function is inspired/lifted from Ivan's code.
    Map a series of nucleotide positions to SARS-CoV-2 genes.
    See https://www.ncbi.nlm.nih.gov/nuccore/MN908947.
    :param pos: Nucleotide position pandas series from VCF
    :param GENE_PROTEIN_POSITIONS_DICT: Dictionary of gene positions from cov_lineages
    :type pos: int
    :return: series containing SARS-CoV-2 chromosome region names at each
    nucleotide position in ``pos``
    """
    # make nucleotide positions series into a df
    df = pos.to_frame()
    pos_column = df.columns[0]
    # add new empty (NaN) columns
    cols_to_use = [pos_column, "gene", "protein_name", "protein_symbol", "protein_id"]
    df = df.reindex(columns = cols_to_use) 

    # loop through all CDS regions in dict to get attributes
    for entry in GENE_PROTEIN_POSITIONS_DICT.keys():

        if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="CDS" and ("protein_alias" in GENE_PROTEIN_POSITIONS_DICT[entry]):

            # extract values from JSON entry
            start = GENE_PROTEIN_POSITIONS_DICT[entry]["start"]
            end = GENE_PROTEIN_POSITIONS_DICT[entry]["end"]
            #aa_start = GENE_PROTEIN_POSITIONS_DICT[entry]["aa_start"]
            #aa_end = GENE_PROTEIN_POSITIONS_DICT[entry]["aa_end"]
            gene = GENE_PROTEIN_POSITIONS_DICT[entry]["gene"]
            product = GENE_PROTEIN_POSITIONS_DICT[entry]["product"]
            protein_alias = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_alias"]
            protein_name = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_name"]["label"]
            protein_symbol = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_symbol"]["label"]
            protein_id = GENE_PROTEIN_POSITIONS_DICT[entry]["protein_id"]
            locus_tag = GENE_PROTEIN_POSITIONS_DICT[entry]["locus_tag"]
    
            # fill in attributes for mutations in this CDS region
            cds_mask = df[pos_column].astype(int).between(start, end, inclusive="both")
            df.loc[cds_mask, "gene"] = gene
            df.loc[cds_mask, "product"] = product
            df.loc[cds_mask, "protein_alias"] = ','.join(protein_alias)
            df.loc[cds_mask, "protein_name"] = protein_name
            df.loc[cds_mask, "protein_symbol"] = protein_symbol
            df.loc[cds_mask, "protein_id"] = protein_id
            df.loc[cds_mask, "locus_tag"] = locus_tag

        if GENE_PROTEIN_POSITIONS_DICT[entry]["type"]=="gene":

            # extract values from JSON entry
            start = GENE_PROTEIN_POSITIONS_DICT[entry]["start"]
            end = GENE_PROTEIN_POSITIONS_DICT[entry]["end"]
            gene_name = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_name"]["label"]
            gene_symbol = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_symbol"]["label"]
            strand_orientation = GENE_PROTEIN_POSITIONS_DICT[entry]["strand_orientation"]["label"]
            gene_orientation = GENE_PROTEIN_POSITIONS_DICT[entry]["gene_orientation"]["label"]
            
            # fill in attributes for mutations in this gene region
            gene_mask = df[pos_column].astype(int).between(start, end, inclusive="both")
            df.loc[gene_mask, "gene_name"] = gene_name
            df.loc[gene_mask, "gene_symbol"] = gene_symbol
            df.loc[gene_mask, "strand_orientation"] = strand_orientation
            df.loc[gene_mask, "gene_orientation"] = gene_orientation

    # label all mutations that didn't belong to any gene as "intergenic"
    df.loc[df["gene"].isna(), "gene"] = "intergenic"
    # label all mutations that didn't belong to any protein as "n/a"
    df = df.fillna("n/a")

    return(df)


def clade_defining_threshold(threshold, df, sample_size):
    """Specifies the clade_defining attribute as True if AF >
    threshold, False if AF <= threshold, and n/a if the VCF is for a
    single genome """

    if sample_size == 1:
        df["clade_defining"] = "n/a"
    else:
        df.loc[df.alternate_frequency > threshold, "clade_defining"] = "True"
        df.loc[df.alternate_frequency <= threshold, "clade_defining"] = "False"
        
    return df


def add_alias_names(df, GENE_PROTEIN_POSITIONS_DICT):
    '''Creates alias names for Orf1ab mutations, reindexing the amino acid numbers.'''
    df.loc[:, 'alias'] = 'n/a'

    # get list of all NSP, 3CL, and PlPro proteins in the file:
    alias_mask = (df['gene_symbol'].str.contains("orf1ab")) & (df['mat_pep']!='n/a')
    nsps_list = sorted(list(set(df[alias_mask]['mat_pep'].tolist())))
    if len(nsps_list) > 0:
        
        ## note: gene and protein_name are based on our gene positions JSON

        # split up all names in alias_mask into letter-number-letter columns
        # hacky workaround to fix later: in rows that begin with a number, add "PLACEHOLDER" to the front before splitting them up, to stop NaNs
        df.loc[df['original_mutation_description'].str[0].str.isdigit(), 'original_mutation_description'] = "PLACEHOLDER" + df['original_mutation_description'].astype(str)

        # split at underscores
        if df['original_mutation_description'].str.contains("_").any():
            df[['mutation_1', 'mutation_2']] = df['original_mutation_description'].str.split('_', expand=True)
        else:
            df['mutation_1'] = df['original_mutation_description']
            df['mutation_2'] = None
        
        df[['1_start', '1_num', '1_end']] = df['mutation_1'].str.extract('([A-Za-z]+)(\\d+\\.?\\d*)([A-Za-z]*)', expand=True)
        df['1_num'] = df['1_num'].fillna(0)
        df['1_newnum'] = 0
        
        df[['2_start', '2_num', '2_end']] = df['mutation_2'].str.extract('([A-Za-z]+)(\\d+\\.?\\d*)([A-Za-z]*)', expand=True)
        df['2_num'] = df['2_num'].fillna(0)
        df['2_newnum'] = 0
        
        # for each nsp in nsps_list, operate on the number column based on the nsp start coordinates
        for nsp in nsps_list:
            if nsp=='PL_proPLpro':
                nsp='PL_pro'
            nsp_start_aa = int(GENE_PROTEIN_POSITIONS_DICT[nsp]["aa_start"])
            nsp_protein_id = GENE_PROTEIN_POSITIONS_DICT[nsp]["protein_id"]
            nsp_mask = df['mat_pep']==nsp
            # add alias_protein_id
            df.loc[nsp_mask, 'alias_protein_id'] = nsp_protein_id
            # for each half of the mutation name...
            # update the numeric part
            df.loc[nsp_mask, '1_newnum'] = df['1_num'].astype(int) - nsp_start_aa + 1
            df.loc[nsp_mask, '2_newnum'] = df['2_num'].astype(int) - nsp_start_aa + 1
            # put the three columns back together into a column called '1_alias' or '2_alias'
            df.loc[nsp_mask, '1_alias'] = df['1_start'] + df['1_newnum'].astype(str) + df['1_end']
            df.loc[nsp_mask, '2_alias'] = df['2_start'] + df['2_newnum'].astype(str) + df['2_end']
        
        # put both halves of the alias back together with a new underscore in the middle, for all ORF1ab rows
        df.loc[alias_mask, 'alias'] = df['1_alias'].astype(str) + '_' + df['2_alias'].astype(str)

        # remove the placeholder
        df['original_mutation_description'] = df['original_mutation_description'].str.replace("PLACEHOLDER","")
        df['alias'] = df['alias'].str.replace("PLACEHOLDER","")
        
        # remove nans
        df['alias'] = df['alias'].str.replace("nan_nan", "")
        df['alias'] = df['alias'].str.replace("_nan", "")
        
        cols_to_drop = [col for col in df.columns if (("1" in col) or ("2" in col))]
        df = df.drop(columns=cols_to_drop)
    
    return df
        

