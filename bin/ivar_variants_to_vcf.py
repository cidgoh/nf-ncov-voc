#!/usr/bin/env python
import os
import sys
import re
import errno
import argparse


def parse_args(args=None):
    description = 'Convert iVar variants tsv file to vcf format.'
    epilog = """
                Example usage: python ivar_variants_to_vcf.py <FILE_IN> 
                <FILE_OUT>
            """

    parser = argparse.ArgumentParser(description=description,
                                     epilog=epilog)
    parser.add_argument('FILE_IN',
                        help="Input tsv file.")
    parser.add_argument('FILE_OUT',
                        help="Full path to output vcf file.")
    parser.add_argument('-po', '--pass_only',
                        dest="PASS_ONLY",
                        help="Only output variants that PASS all "
                             "filters.",
                        action='store_true')
    parser.add_argument('-af', '--allele_freq_thresh',
                        type=float,
                        dest="ALLELE_FREQ_THRESH",
                        default=0,
                        help="Only output variants where allele "
                             "frequency greater than this number ("
                             "default: 0).")

    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def ivar_variants_to_vcf(file_in, file_out, pass_only=False,
                         min_alt_freq=0):
    filename = os.path.splitext(file_in)[0]
    header = ('##fileformat=VCFv4.2\n'
              '##source=iVar\n'
              '##contig=<ID=MN908947.3,length=29903>\n'
              '##INFO=<ID=DP,Number=1,Type=Integer,'
              'Description="Total read depth at the locus">\n'
              '##FILTER=<ID=PASS,Description="Result of p-value <= '
              '0.05">\n'
              '##FILTER=<ID=FAIL,Description="Result of p-value > '
              '0.05">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,'
              'Description="Genotype">\n'
              '##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Count'
              ' of full observations of the reference haplotype.">\n'
              '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,'
              'Description="Depth of reference base on reverse '
              'reads">\n'
              '##FORMAT=<ID=QR,Number=1,Type=Float,'
              'Description="Reference allele quality sum in phred">\n'
              '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Count'
              ' of full observations of this alternate haplotype."> \n'
              '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,'
              'Description="Depth of alternate base on reverse '
              'reads">\n'
              '##FORMAT=<ID=ALT_QUAL,Number=1,Type=String, '
              'Description="Mean quality of alternate base">\n'
              '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,'
              'Description="Mean quality of reference base">\n'
              '##FORMAT=<ID=QA,Number=A,Type=Integer,'
              'Description="Alternate allele quality sum in phred">\n'
              '##FORMAT=<ID=AF,Number=A,Type=Float,'
              'Description="Estimated allele frequency in the range '
              '(0,1]">\n')
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' \
              '\t' + filename + '\n'

    var_list = []
    var_count_dict = {'SNP': 0, 'INS': 0, 'DEL': 0}
    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)
    f_out = open(file_out, 'w')
    f_out.write(header)
    with open(file_in) as f:
        for line in f:
            print(line)
            if not re.match("REGION", line):
                line = re.split("\t", line)
                chromosome = line[0]
                position = line[1]
                id_var = '.'
                ref = line[2]
                alt = line[3]
                var_type = 'SNP'
                if alt[0] == '+':
                    alt = ref + alt[1:]
                    var_type = 'INS'
                elif alt[0] == '-':
                    ref += alt[1:]
                    alt = line[2]
                    var_type = 'DEL'
                quality = '.'
                pass_test = line[13]
                if pass_test == 'TRUE':
                    filter_flag = 'PASS'
                else:
                    filter_flag = 'FAIL'
                info = 'DP=' + line[11]
                info_format = \
                    'GT:RO:REF_RV:REF_QUAL:AO:ALT_RV' \
                    ':ALT_QUAL:AF'
                sample = '1:' + line[4] + ':' + line[5] + ':' + line[
                    6] + ':' + line[7] + ':' + line[8] + ':' + line[
                             9] + ':' + line[10]
                o_line = chromosome + '\t' + position + '\t' + id_var \
                    + '\t' + ref + '\t' + alt + '\t' + quality + '\t' \
                    + filter_flag + '\t' + info + '\t' + \
                    info_format + '\t' + sample + '\n'
                write_line = True
                if pass_only and filter_flag != 'PASS':
                    write_line = False
                if float(line[10]) < min_alt_freq:
                    write_line = False
                if (chromosome, position, ref, alt) in var_list:
                    write_line = False
                else:
                    var_list.append((chromosome, position, ref, alt))
                if write_line:
                    var_count_dict[var_type] += 1
                    f_out.write(o_line)
    f_out.close()

    # Print variant counts to pass to MultiQC

    var_count_list = [(k, str(v)) for k, v in
                      sorted(var_count_dict.items())]
    print('\t'.join(['sample'] + [x[0] for x in var_count_list]))
    print('\t'.join([filename] + [x[1] for x in var_count_list]))


def main(args=None):
    args = parse_args(args)
    ivar_variants_to_vcf(args.FILE_IN, args.FILE_OUT, args.PASS_ONLY,
                         args.ALLELE_FREQ_THRESH)


if __name__ == '__main__':
    sys.exit(main())
