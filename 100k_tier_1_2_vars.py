"""
v1.0 - AB 2019/05/30
Requirements:
    Python 3.7
    JellyPy
    GeL Report Models (v6 or higher)

usage: 100k_tier_1_2_vars.py [-h] -i IR_ID -p PROBAND_ID

Prints tier 1 and 2 variant details to stdout for a given 100k case

optional arguments:
  -h, --help            show this help message and exit
  -i IR_ID, --ir_id IR_ID
                        GeL Interpretation Request ID in format 12345-1
  -p PROBAND_ID, --proband_id PROBAND_ID
                        GeL participant ID for proband
"""

from configparser import ConfigParser
import os
import sys
import argparse
import glob
import re
import json
import datetime
# Append JellyPy to python path, needed when running via paramiko from Windows
sys.path.append('/home/mokaguys/Apps/JellyPy')
from pyCIPAPI.interpretation_requests import get_interpretation_request_list, get_interpretation_request_json
# Import InterpretedGenome from GeLReportModels v6.0
from protocols.reports_6_0_0 import InterpretedGenome

config = ConfigParser()
config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini"))

def process_arguments():
    """
    Uses argparse module to define and handle command line input arguments and help menu
    """
    # Create ArgumentParser object. Description message will be displayed as part of help message if script is run with -h flag
    parser = argparse.ArgumentParser(description='Prints tier 1 and 2 variant details to stdout for a given 100k case')
    # Define the arguments that will be taken.
    parser.add_argument('-i', '--ir_id', required=True, help='GeL Interpretation Request ID in format 12345-1')
    parser.add_argument('-p', '--proband_id', required=True, help='GeL participant ID for proband')
    # Return the arguments
    return parser.parse_args()


def group_vars_by_cip(interpreted_genomes_json):
    """
    Groups variants by CIP provider
    Args:
        interpreted_genomes_json: interpreted genome JSON from CIP API
    Returns:
        Nested dictionary of variants grouped by CIP and cip version. {key = CIP value = {key = cip_version, value = list_of_variants}}

    """
    # Note this does not pull out CNVs/SVs
    vars_by_cip = {}
    # possible values for cip at time of writing:
    # omicia, congenica, nextcode, genomics_england_tiering, illumina, exomiser
    # There will be a separate interepreted genome for each cip used
    for ig in interpreted_genomes_json:
        # Convert the interpreted genome JSON into InterpretedGenome object from GeL Report Models v6.0
        ig_obj = InterpretedGenome.fromJsonDict(ig['interpreted_genome_data'])
        # cip provider stored in the interpretationService field.
        # Store the list of reported variants for that cip
        cip = ig_obj.interpretationService.lower()
        cip_version = int(ig['cip_version'])
        # If cip not already in dictionary, add it in with an empty dictionary as value
        vars_by_cip[cip] = vars_by_cip.setdefault(cip, {})
        # If CIP is present multiple times each should have it's own version number
        # However do a quick test to make sure this is true and error out if not
        if cip_version in vars_by_cip[cip]:
            sys.exit(f"Multiple interpreted genomes with version number '{cip_version}' for interpretation service '{cip}'")
        # If there are variants, add the variant list for that cip/version to dictionary.
        if ig_obj.variants:
            vars_by_cip[cip][cip_version] = ig_obj.variants
        # If there aren't any variants just store empty list
        else:
            vars_by_cip[cip][cip_version] = []
    return vars_by_cip


def group_vars_by_tier(variants_json):
    """
    Groups variants according to their GeL tier
    Args:
        variants_json: list of GeL report model v6 variant objects
    Returns:
        Dictionary of variants grouped by tier {key = tier, value = list of variant objects}
    """
    # Takes list of variants and groups them by tier (Note this doesn't include SVs/CNVs)
    tiered_vars = {
        'TIER1': [],
        'TIER2': [],
        'TIER3': [],
        'OTHER': []
        }
    for variant in variants_json:
        # Log the tier for each report event for a variant
        tiers = []
        for reportevent in variant.reportEvents:
            # Possible values for tier in report model v6 are NONE, TIER1, TIER2, TIER3, TIER4, TIER5, TIERA, TIERB
            if reportevent.tier.upper() in ('TIER1', 'TIER2', 'TIER3'):
                # Record the tier number
                tiers.append(int(reportevent.tier.strip('TIER')))
            else:
                # Value of 4 used to represent any other or no tier.
                tiers.append(4)
        # The lowest number represents the highest ranked tier for that variant (e.g. tier 1 ranks higher than tier 2)
        top_rank_tier = min(tiers)
        # Record the variant in the dictionary based on it's highest ranked tier
        if top_rank_tier == 4:
            tiered_vars['OTHER'].append(variant)
        else:
            tiered_vars[f'TIER{top_rank_tier}'].append(variant)
    return tiered_vars


def get_tiered_vars(ir_json):
    """
    Returns GeL tiered variants from an interpretation request json
    Args:
        ir_json: GeL interpretation request JSON conforming to GeL report models v6
    Returns:
        Dictionary of GeL report model variant objects grouped by tier. Keys = 'TIER1', 'TIER2', 'TIER3' and 'OTHER'
    """
    # Create a dictionary of variant lists (value) grouped by CIP (key - lowercase)
    vars_by_cip = group_vars_by_cip(ir_json['interpreted_genome'])
    # There may be multiple interpreted genomes for a single CIP, so take the one with the highest cip_version number
    max_version = max(vars_by_cip['genomics_england_tiering'].keys())
    # Group variants from genomics_england_tiering by tier
    return group_vars_by_tier(vars_by_cip['genomics_england_tiering'][max_version])


def zygosity_to_vcf(gel_zygosity):
    """
    Convert GeL zygosity (as defined in v6 of GeL report models) to VCF format:
    https://gelreportmodels.genomicsengland.co.uk/html_schemas/org.gel.models.report.avro/6.0.1/InterpretedGenome.html#/schema/schemas%2FAVPRs%2Fbuild%2FInterpretedGenome.avpr/org.gel.models.report.avro.Zygosity
    Args:
        gel_zygosity: string representation of zygosity as stored in JSON e.g. 'heterozygous'
    Returns:
        VCF representation of zygosity e.g. '0/1'
    """
    vcf_lookup = {
        'reference_homozygous': '0/0',
        'heterozygous': '0/1',
        'alternate_homozygous': '1/1',
        'missing': './.',
        'half_missing_reference': './0',
        'half_missing_alternate': './1',
        'alternate_hemizigous': '1',
        'reference_hemizigous': '0',
        'na': 'na',
        'unk': 'unk'
    }
    return vcf_lookup[gel_zygosity]


def get_ir_json(ir_id, ir_version):
    """
    Get the interpretation request JSON for a case and return as a python json dictionary like object
    For speed this function will use a local cached json rather than downloading a fresh copy, provided
    the last_modified timestamp matches that found in the interpretation request list endpoint
    Args:

    """
    # Get list of local JSON file names (pulled using GeL2MDT)
    local_jsons = ';'.join(glob.glob('{local_cache}/{ir_id}-{ir_version}-*.json'.format(
            local_cache=config.get("GEL2MDT", "CIP-API-CACHE"),
            ir_id=ir_id,
            ir_version=ir_version
        )
    ))
    # GeL2MDT appends sequential version numbers to the end of the JSON each time a new copy is downloaded
    # Capture the GeL2MDT version number from end of each JSON file name to find latest file
    versions = re.findall(r'{ir_id}-{ir_version}-(\d+).json'.format(ir_id=ir_id, ir_version=ir_version), local_jsons)
    # If there's no local cache files, or filenames don't match expected format, versions will be an empty list so skip next block
    if versions:
        # Find highest version and construct filepath
        max_version = max(list(map(int, versions)))
        latest_cache_path = '{local_cache}/{ir_id}-{ir_version}-{max_version}.json'.format(
                local_cache=config.get("GEL2MDT", "CIP-API-CACHE"),
                ir_id=ir_id,
                ir_version=ir_version,
                max_version=max_version
            )
        # Read into a JSON object
        with open(latest_cache_path, 'r') as cache_json:
            ir_json = json.load(cache_json)
        # Get overview of case from interpretation request list endpoint (this is much faster than pulling the full JSON)
        ir_list = get_interpretation_request_list(interpretation_request_id=ir_id, version=ir_version)
        # Check that this only returns one result
        if len(ir_list) != 1:
            sys.exit(f"Length of interpretation request list was {list_length} which is different to expected length 1".format(list_length=len(ir_list)))
        # Get the last modified timestamp from cip api and local cache and convert to datetime objects
        last_modified_api = datetime.datetime.strptime(ir_list[0]['last_modified'], "%Y-%m-%dT%H:%M:%S.%fZ")
        last_modified_cache = datetime.datetime.strptime(ir_json['last_modified'], "%Y-%m-%dT%H:%M:%S.%fZ")
        # Check if they match. If they do, return the local cache JSON.
        if last_modified_api == last_modified_cache:
            return ir_json
    # If there is no local cache, or it's out of date, download a new JSON from CIP-API and return
    return get_interpretation_request_json(ir_id, ir_version, reports_v6=True)


def main():
    # Return arguments
    args = process_arguments()
    # Pull out interpretation request JSON
    ir_json = get_ir_json(args.ir_id.split('-')[0], args.ir_id.split('-')[1])
    # Group variants by tier
    tiered_vars = get_tiered_vars(ir_json)
    # Loop through tier 1 and tier 2 variants
    for var in tiered_vars['TIER1'] + tiered_vars['TIER2']:
        # Need to get genotype for the proband.
        # There is a separate 'variant_call' for each family member, so find proband using participant ID
        for variant_call in var.variantCalls:
            if variant_call.participantId == args.proband_id:
                # Capture proband genotype and convert to VCF format
                gt = zygosity_to_vcf(variant_call.zygosity)
        # The other details needed are stored in 'variantCoordinates'
        # Print variant details in tab separated list to stdout
        vc = var.variantCoordinates
        print(f"{vc.assembly}\t{vc.chromosome}\t{vc.position}\t{vc.reference}\t{vc.alternate}\t{gt}")

if __name__ == '__main__':
    main()
