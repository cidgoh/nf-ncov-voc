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