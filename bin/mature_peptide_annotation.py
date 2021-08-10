#!/usr/bin/env python3
import argparse
import gffutils
from cyvcf2 import VCF, Writer


def parse_args():
    parser = argparse.ArgumentParser(
        description='Adds mature peptide annotation to VCF files ')
    parser.add_argument('--vcf_file', type=str, default=None,
                        help='Variant calling output file in VCF '
                             'format')
    parser.add_argument('--annotation_file', type=str, default=None,
                        help='Annotation file (SARS-CoV-2) in GFF '
                             'format')
    parser.add_argument('--output_vcf', type=str, default=None,
                        help='Output VCF file')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    gene_protein = {
        "orf1ab": "cds-YP_009724389.1",
        "S": "cds-YP_009724390.1",
        "ORF3a": "cds-YP_009724391.1",
        "E": "cds-YP_009724392.1",
        "M": "cds-YP_009724393.1",
        "ORF6": "cds-YP_009724394.1",
        "ORF7a": "cds-YP_009724395.1",
        "ORF7b": "cds-YP_009725318.1",
        "ORF8": "cds-YP_009724396.1",
        "N": "cds-YP_009724397.2",
        "ORF10": "cds-YP_009725255.1",
    }
    db = gffutils.create_db(
        args.annotation_file, 'ncov_annotation.db',
        force=True, merge_strategy="warning")
    data_vcf = VCF(args.vcf_file)
    data_vcf.add_info_to_header(
        {'ID': 'mat_pep_id', 'Description': 'Mature Peptide ID',
         'Type': 'String', 'Number': '.'})
    data_vcf.add_info_to_header(
        {'ID': 'mat_pep_desc',
         'Description': 'Mature Peptide Description',
         'Type': 'String', 'Number': '.'})
    data_vcf.add_info_to_header(
        {'ID': 'mat_pep_acc',
         'Description': 'Mature Peptide Accession Number',
         'Type': 'String', 'Number': '.'})
    output_file_name = args.output_vcf
    w = Writer(output_file_name, data_vcf)
    gene = db[gene_protein["orf1ab"]]
    # gene = db[gene_protein[record.INFO.get('EFF').split("|")[5]]]

    for record in data_vcf:
        record.INFO["mat_pep_id"] = "n/a"
        record.INFO["mat_pep_desc"] = "n/a"
        record.INFO["mat_pep_acc"] = "n/a"
        for i in db.children(gene,
                             featuretype='mature_protein_region_of_CDS',
                             order_by='start'):
            if i.start <= record.POS <= i.end:
                record.INFO["mat_pep_id"] = ",".join(i.attributes[
                                                         'product'])
                record.INFO["mat_pep_desc"] = ",".join(i.attributes[
                                                           'Note'])
                record.INFO["mat_pep_acc"] = ",".join(i.attributes[
                                                          'protein_id'])
        w.write_record(record)
    w.close()
    data_vcf.close()
