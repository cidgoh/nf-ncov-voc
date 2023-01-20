import pandas as pd

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
    
    

def find_sample_size(table, lineage):


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
            filename_to_match = args.vcffile.split(".sorted")[0] \
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
    info = pd.concat([info, eff_info], axis=0)

    return(info)
