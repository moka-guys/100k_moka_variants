#!/usr/bin/env python3
"""
v1.0 - AB 2019/05/30

usage: 100k_variant_validator.py [-h] -b BUILD -c CHR -p POS -r REF -a ALT

Queries VariantValidator API to extract details for a given variant for
entering to Moka

optional arguments:
  -h, --help            show this help message and exit
  -b BUILD, --build BUILD
                        Genome build
  -c CHR, --chr CHR     Chromosome
  -p POS, --pos POS     Position
  -r REF, --ref REF     Reference base(s)
  -a ALT, --alt ALT     Alternate base(s)
"""

import os
import sys
import argparse
import requests
import datetime
from configparser import ConfigParser

# Read config file (must be called config.ini and stored in same directory as script)
config = ConfigParser()
config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini"))

def process_arguments():
    """
    Uses argparse module to define and handle command line input arguments and help menu
    """
    # Create ArgumentParser object. Description message will be displayed as part of help message if script is run with -h flag
    parser = argparse.ArgumentParser(
        description='Queries VariantValidator API to extract details for a given variant for entering to Moka'
    )
    # Define the arguments that will be taken.
    parser.add_argument('-b', '--build', required=True, help='Genome build')
    parser.add_argument('-c', '--chr', required=True, help='Chromosome')
    parser.add_argument('-p', '--pos', required=True, help='Position')
    parser.add_argument('-r', '--ref', required=True, help='Reference base(s)')
    parser.add_argument('-a', '--alt', required=True, help='Alternate base(s)')
    # Return the arguments
    return parser.parse_args()


class Variant(object):
    """
    Represents a genetic variant, pulling details from VariantValidator
    """
    def __init__(self):
        """
        The fields below should be populated from VariantValidator using the get_info method.
        """
        self.chr37 = ''
        self.pos37 = ''
        self.ref37 = ''
        self.alt37 = ''
        self.chr38 = ''
        self.pos38 = ''
        self.ref38 = ''
        self.alt38 = ''
        self.transcripts = []
        self.var_val_version = ''

    def parse_intergenic(self, vv_json):
        """
        Parses VariantValidator JSON for intergenic variants
        """
        # For intergenic variants, we just need b37 + b38 coords, ref and alt
        # Can get these from field Intergenic_Variant_1
        # Author indicated the key will be made lower case in the future for consistency, so try both:
        # https://github.com/openvar/variantValidator/issues/57
        if 'Intergenic_Variant_1' in vv_json:
            loci_json = vv_json['Intergenic_Variant_1']['primary_assembly_loci']
        elif 'intergenic_variant_1' in vv_json:
            loci_json = vv_json['intergenic_variant_1']['primary_assembly_loci']
        # Capture b37 coordinates
        if 'grch37' in loci_json:
            grch37_vcf = loci_json['grch37']['vcf']
            self.chr37 = grch37_vcf['chr']
            self.pos37 = grch37_vcf['pos']
            self.ref37 = grch37_vcf['ref']
            self.alt37 = grch37_vcf['alt']
        # Capture b38 coordinates
        if 'grch38' in loci_json:
            grch38_vcf = loci_json['grch38']['vcf']
            self.chr38 = grch38_vcf['chr']
            self.pos38 = grch38_vcf['pos']
            self.ref38 = grch38_vcf['ref']
            self.alt38 = grch38_vcf['alt']

    def parse_gene_variant(self, vv_json):
        """
        Parses VariantValidator JSON for gene variants
        """
        ambiguous_37 = False
        ambiguous_38 = False
        # For variants that map to a transcript, store the json for each transcript in a list
        transcript_data = [vv_json[key] for key in vv_json.keys() if key not in ['flag', 'metadata']]
        # Loop through each transcript JSON, extract relevant fields and store in list of dictionaries
        for tx in transcript_data:
            # Capture b37 coordinates
            if 'grch37' in tx['primary_assembly_loci']:
                grch37_vcf = tx['primary_assembly_loci']['grch37']['vcf']
                # if we've already captured the coordinates from a previous transcript, check they match, if they don't set ambiguous flag to True
                if self.chr37 or self.pos37 or self.ref37 or self.alt37:
                    if not (
                        self.chr37 == grch37_vcf['chr'] and
                        self.pos37 == grch37_vcf['pos'] and
                        self.ref37 == grch37_vcf['ref'] and
                        self.alt37 == grch37_vcf['alt']
                            ):
                                ambiguous_37 = True
                # if we haven't captured the coordinates yet, capture them
                else:
                    self.chr37 = grch37_vcf['chr']
                    self.pos37 = grch37_vcf['pos']
                    self.ref37 = grch37_vcf['ref']
                    self.alt37 = grch37_vcf['alt']
            # Capture b38 coordinates
            if 'grch38' in tx['primary_assembly_loci']:
                grch38_vcf = tx['primary_assembly_loci']['grch38']['vcf']
                # if we've already captured the coordinates from a previous transcript, check they match, if they don't set ambiguous flag to True
                if self.chr38 or self.pos38 or self.ref38 or self.alt38:
                    if not (
                        self.chr38 == grch38_vcf['chr'] and
                        self.pos38 == grch38_vcf['pos'] and
                        self.ref38 == grch38_vcf['ref'] and
                        self.alt38 == grch38_vcf['alt']
                            ):
                                ambiguous_38 = True
                # if we haven't captured the coordinates yet, capture them
                else:
                    self.chr38 = grch38_vcf['chr']
                    self.pos38 = grch38_vcf['pos']
                    self.ref38 = grch38_vcf['ref']
                    self.alt38 = grch38_vcf['alt']
            # If there's a gene symbol, capture it
            gene = ''
            if tx['gene_symbol']:
                gene = tx['gene_symbol']
            # Capture the transcript name from start of HGVS
            transcript = tx['hgvs_transcript_variant'].split(':')[0]
            # Capture transcript HGVS minus the accession prefix e.g NM_153240.4:c.3762_3764dup would become c.3762_3764dup
            hgvst = tx['hgvs_transcript_variant'].split(':')[1]
            # If there's protein HGVS, capture the p. nomenclature single letter representation (slr).
            hgvsp = ''
            if tx['hgvs_predicted_protein_consequence']:
                hgvsp_full = tx['hgvs_predicted_protein_consequence']['slr']
                # Remove the NP protein accession prefix, and parentheses to make consistent with Ingenuity exported VCF
                # i.e. NP_000079.2:p.(G197C) would become p.G197C
                hgvsp = hgvsp_full.split(':')[1].replace('(', '').replace(')', '')
            # Store as dictionary in self.transcripts list
            self.transcripts.append({'gene': gene, 'transcript': transcript, 'hgvst': hgvst, 'hgvsp': hgvsp})
        # If either sets of coordinates are ambiguous, set back to empty strings
        if ambiguous_37:
            self.chr37 = ''
            self.pos37 = ''
            self.ref37 = ''
            self.alt37 = ''
        if ambiguous_38:
            self.chr38 = ''
            self.pos38 = ''
            self.ref38 = ''
            self.alt38 = ''

    def get_info(self, build, chr, pos, ref, alt):
        """
        Retrieves variant details required for Moka using VariantValidator API

        Args:
            build: Genome build (GRCh38, GRCh37)
            chr: Chromosome
            pos: Position
            ref: Reference bases(s)
            alt: Alternate base(s)
        """
        endpoint = f'{config.get("VariantValidator", "BASE_URL")}/{build}/{chr}-{pos}-{ref}-{alt}/all'
        r = requests.get(endpoint)
        # Error if response code not 200
        if r.status_code != 200:
            sys.exit(f'Response code {r.status_code} from VariantValidator')
        # Expecting flag to either be gene_variant or intergenic. Raise error for others.
        # https://github.com/openvar/variantValidator/issues/57
        if r.json()['flag'] not in ['gene_variant', 'intergenic']:
            sys.exit(f"{r.json()['flag']} flag returned from VaraintValidator")
        # Call the appropriate method depending on the variant type
        flag_parser = {'intergenic': self.parse_intergenic, 'gene_variant': self.parse_gene_variant}
        parser = flag_parser[r.json()['flag']]
        parser(r.json())
        # Some cases can't be lifted over, but check that there is at least a section for one of the two builds
        if not self.pos37 and not self.pos38:
            sys.exit('Missing coordinates for both GRCh37 and GRCh38')
        # Record the variant validator version
        self.var_val_version = r.json()['metadata']['variantvalidator_version']


def main():
    # Get arguments
    args = process_arguments()
    # Pull variant info from VariantValidator
    v = Variant()
    v.get_info(
        args.build,
        args.chr,
        args.pos,
        args.ref,
        args.alt
    )
    # Format transcript level information into a list of semi-colon separated strings
    tx_strs = [f"{tx['gene']};{tx['transcript']};{tx['hgvst']};{tx['hgvsp']}" for tx in v.transcripts]
    # Print the transcript info to stdout as a tab-separated string
    # Transcript strings are joined using commas and printed in the final field
    #Get datetime stamp to add to logifle name
    filename = datetime.datetime.now()
    datetime_stamp = filename.strftime("%y%m%d_%H%M%S")
    #Save stdout to logfile 
    sys.stdout = open('/home/mokaguys/Apps/100k_moka_variants/variant_validator/%s_variantvalidator.tsv' % datetime_stamp, 'w')
    print(f"{v.chr37}\t{v.pos37}\t{v.ref37}\t{v.alt37}\t{v.chr38}\t{v.pos38}\t{v.ref38}\t{v.alt38}\t{','.join(tx_strs)}\t{v.var_val_version}")


if __name__ == '__main__':
    main()
