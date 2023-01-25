import pandas as pd
import numpy as np

def parse_variant_file(dataframe):
    who_lineages = []
    for var in dataframe["pango_lineage"]:
        if "," in var:
            for temp in var.split(","):
                if not "[" in var:
                    who_lineages.append(temp)
                    who_lineages.append(temp+".*")
                else:
                    parent=temp[0]
                    child=temp[2:-1].split("|")
                    for c in child:
                        who_lineages.append(parent + str(c))
                        who_lineages.append(parent+str(c)+".*")
        else:
            who_lineages.append(var)
    return who_lineages
    
    

def find_sample_size(table, lineage, vcf_file):


    if table != 'n/a':
        strain_tsv_df = pd.read_csv(table, delim_whitespace=True,
                                    usecols=['file', 'num_seqs'])

        # Reference mode
        if lineage != 'n/a':
            num_seqs = strain_tsv_df[strain_tsv_df['file'].str.startswith(
                lineage + ".qc.")]['num_seqs'].values
            sample_size = num_seqs[0]

        # user-uploaded fasta
        else:
            filename_to_match = vcf_file.split(".sorted")[0] \
                # looks like "strain.qc"
            num_seqs = strain_tsv_df[strain_tsv_df['file'].str.startswith(
                filename_to_match)]['num_seqs'].values
            sample_size = num_seqs[0]

    # user-uploaded vcf
    elif table == 'n/a' and lineage == 'n/a':
        sample_size = 'n/a'

    return sample_size



    
def parse_INFO(df): # return INFO dataframe with named columns, including EFF split apart
    # parse INFO column
    
    # sort out problematic sites tag formats   
    df['INFO'] = df['INFO'].str.replace('ps_filter;', 'ps_filter=;')
    df['INFO'] = df['INFO'].str.replace('ps_exc;', 'ps_exc=;')
    df['INFO'] = df['INFO'].str.replace('=n/a', '')
  
    # make 'INFO' column easier to extract attributes from:
    # split at ;, form dataframe
    info = df['INFO'].str.split(pat=';').apply(pd.Series)

    for column in info.columns:
        split = info[column].str.split(pat='=').apply(pd.Series)
        title = split[0].drop_duplicates().tolist()[0]
        if isinstance(title, str):
            title = title.lower()
            content = split[1]
            # ignore "tag=" in column content
            info[column] = content
            # make attribute tag as column label
            info.rename(columns={column: title}, inplace=True)

    # parse EFF entry in INFO
    # series: extract everything between parentheses as elements of a
    # list
    eff_info = info['eff'].str.findall('\((.*?)\)')
    # take first element of list
    eff_info = eff_info.apply(pd.Series)[0]
    # split at pipe, form dataframe
    eff_info = eff_info.str.split(pat='|').apply(pd.Series)
    # hgvs names
    hgvs = eff_info[3].str.rsplit(pat='c.').apply(pd.Series)
    eff_info["hgvs_protein"] = hgvs[0].str[:-1]
    eff_info["hgvs_nucleotide"] = 'g.' + hgvs[1] # change 'c.' to 'g.' for nucleotide names

    # change nucleotide names of the form "g.C*4378A" to g.C4378AN;
    # change vcf_gene to "intergenic" here
    asterisk_mask = eff_info["hgvs_nucleotide"].str.contains('\*')
    eff_info["hgvs_nucleotide"][asterisk_mask] = 'g.' + df['REF'] + df['POS'] + \
                                     df['ALT']
    eff_info[5][asterisk_mask] = "intergenic"

    # use nucleotide name where protein name doesn't exist (for
    # 'Name' attribute)
    Names = hgvs[0].str[:-1]
    # fill in empty protein name spaces with nucleotide names ("c."...)
    Names[~Names.str.contains("p.")] = eff_info["hgvs_nucleotide"]
    # take "p." off the protein names
    Names = Names.str.replace("p.", "", regex=True)
    eff_info["Names"] = Names
    
    #rename some eff_info columns
    eff_info = eff_info.rename(columns={5: "vcf_gene", 1: "mutation_type"})

    #concatenate info and eff_info horizontally to return just one object
    info = pd.concat([info, eff_info], axis=1)

    return(info)
    
    
    
def add_variant_information(clade_file, merged_df, sample_size, strain):    
    # get clade_defining status, and then info from clades file
    # load clade-defining mutations file
    ### MZA: need to clean up this and add this into separate function "variant_info" 
    clades = pd.read_csv(clade_file, sep='\t', header=0)
    clades = clades.replace(np.nan, '', regex=True) #this is needed to append empty strings to attributes (otherwise datatype mismatch error)

    # find the relevant pango_lineage line in the clade file that
    # matches args.strain (call this line "var_to_match")
    ### replace this part with func "parse_variant_file" in functions.py
    #available_strains = parse_variant_file(clades)
    print("here", strain)
    cladefile_strain = 'None'
    available_strains = []
    for var in clades['pango_lineage'].tolist():
        if "," in var:
            for temp in var.split(","):
                if "[" not in var:
                    available_strains.append(temp)
                    if strain.startswith(temp):
                        var_to_match = var
                else:
                    parent = temp[0]
                    child = temp[2:-3].split("|")
                    for c in child:
                        available_strains.append(parent + str(c))
                        available_strains.append(parent + str(c) + ".*")
                        if strain.startswith(parent + str(c)):
                            var_to_match = var
        else:
            available_strains.append(var)
            if strain.startswith(var):
                var_to_match = var
                
    print("var_to_match", var_to_match)

    # if strain in available_strains:
    if var_to_match:
        # find the index of the relevant row
        var_index = clades.index[clades['pango_lineage'] == var_to_match].tolist()[0]
        print("var_index", var_index)
        print(clades.loc[var_index])
        # extract status, WHO strain name, etc. from clades file
        who_variant = clades.loc[var_index, 'variant']
        variant_type = clades.loc[var_index, 'variant_type']
        voi_designation_date = clades.loc[var_index, 'voi_designation_date']
        voc_designation_date = clades.loc[var_index, 'voc_designation_date']
        vum_designation_date = clades.loc[var_index, 'vum_designation_date']
        status = clades.loc[var_index, 'status']

        # add attributes from variant info file
        merged_df["#attributes"] = merged_df["#attributes"].astype(
            str) + "variant=" + who_variant + ';' + "variant_type=" + \
                                   variant_type + ';' + "voi_designation_date=" + \
                                   voi_designation_date + ';' + \
                                   "voc_designation_date=" + \
                                   voc_designation_date + ';' + \
                                   "vum_designation_date=" + \
                                   vum_designation_date + ';' + \
                                   "status=" + status + ';'
    else:
        merged_df["#attributes"] = merged_df["#attributes"].astype(
            str) + "clade_defining=n/a;" + "variant=n/a;" + \
                                   "variant_type=n/a;" + \
                                   "voi_designation_date=n/a;" + \
                                   "voc_designation_date=n/a;" + \
                                   "vum_designation_date=n/a;" + \
                                    "status=n/a;"
                                    
                                    
    return(merged_df)

