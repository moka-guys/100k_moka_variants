#!/usr/bin/env python2
"""
v1.0 - AB 2019/05/30

usage: 100k_variant_import.py [-h] --ir_id IR_ID --proband_id PROBAND_ID
                              --ngstest_id NGSTEST_ID --internal_pat_id
                              INTERNAL_PAT_ID

Inserts tier 1 and 2 variants to Moka for a 100KGP case

optional arguments:
  -h, --help            show this help message and exit
  --ir_id IR_ID         GeL Interpretation Request ID in format 12345-1
  --proband_id PROBAND_ID
                        GeL participant ID for proband
  --ngstest_id NGSTEST_ID
                        Moka NGSTestID
  --internal_pat_id INTERNAL_PAT_ID
                        Moka InternalPatientID
"""
import os
import sys
import argparse
import subprocess
import datetime
from ConfigParser import ConfigParser
import Tkinter
import tkMessageBox
import paramiko
import pyodbc

# Read config file (must be called config.ini and stored in same directory as script)
config = ConfigParser()
config.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini"))


def process_arguments():
    """
    Uses argparse module to define and handle command line input arguments and help menu
    """
    # Create ArgumentParser object. Description message will be displayed as part of help message if script is run with -h flag
    parser = argparse.ArgumentParser(description='Inserts tier 1 and 2 variants to Moka for a 100KGP case')
    # Define the arguments that will be taken.
    parser.add_argument('--ir_id', required=True, help='GeL Interpretation Request ID in format 12345-1')
    parser.add_argument('--proband_id', required=True, help='GeL participant ID for proband')
    parser.add_argument('--ngstest_id', required=True, help='Moka NGSTestID')
    parser.add_argument('--internal_pat_id', required=True, help='Moka InternalPatientID')
    # Return the arguments
    return parser.parse_args()


class MokaConnector(object):
    """
    pyodbc connection to Moka database for use by other functions
    """
    def __init__(self):
        self.cnxn = pyodbc.connect('DRIVER={{SQL Server}}; SERVER={server}; DATABASE={database};'.format(
            server=config.get("MOKA", "SERVER"),
            database=config.get("MOKA", "DATABASE")
            ),
            autocommit=True
        )
        self.cursor = self.cnxn.cursor()

    def __del__(self):
        """
        Close connection when object destroyed
        """
        self.cnxn.close()

    def execute(self, sql):
        """
        Execute SQL, without catching return values (INSERT, UPDATE etc.)
        """
        self.cursor.execute(sql)

    def fetchall(self, sql):
        """
        Execute SQL catching all return records (SELECT etc.)
        """
        return self.cursor.execute(sql).fetchall()

    def fetchone(self, sql):
        """
        Execute SQL catching one returned record (SELECT etc.)
        """
        return self.cursor.execute(sql).fetchone()


class VariantAdder100KGP(object):
    """
    Methods for adding 100KGP variants to Moka
    Args:
        ngstest_id: Moka NGStestID from NGStest table for 100K case
        internal_pat_id: Moka InternalPatientID from Patients table
        moka_connection: MokaConnector object
    """
    def __init__(self, ngstest_id, internal_pat_id, moka_connection):
        self.internal_pat_id = internal_pat_id
        self.ngstest_id = ngstest_id
        self.mc = moka_connection
        self.prev_vars = self.lookup_prev_vars()
        self.moka_chr = self.lookup_chr()
        self.no_hgncid = []

    def lookup_prev_vars(self):
        """
        Find details of any variants already imported into Moka for this test so we can make sure we don't import them again.
        """
        sql = (
            "SELECT NGSVariant.ChrID, NGSVariant.Position_hg19, NGSVariant.ref, NGSVariant.alt, "
                    "NGSVariant.chr_id_38, NGSVariant.position_38, NGSVariant.ref_38, NGSVariant.alt_38 "
              "FROM NGSVariant "
             "WHERE NGSVariant.NGSTestID = {ngstest_id}"
            ).format(ngstest_id=self.ngstest_id)
        moka_prev_vars = self.mc.fetchall(sql)
        prev_vars = set([])
        for variant in moka_prev_vars:
            # If the variant has build 37 coordinates, store them
            if variant.Position_hg19:
                prev_vars.add(
                        ('GRCh37', str(variant.ChrID), str(variant.Position_hg19), variant.ref, variant.alt)
                    )
            # If the variant had build 38 coordinates, store them
            if variant.position_38:
                prev_vars.add(
                        ('GRCh38', str(variant.chr_id_38), str(variant.position_38), variant.ref_38, variant.alt_38)
                    )
        return prev_vars

    def lookup_chr(self):
        """
        Create a chromosome ID lookup dictionary from Moka Chromosome table
        """
        sql = "SELECT ChrID, Chr FROM Chromosome;"
        moka_chr_all = self.mc.fetchall(sql)
        moka_chr = {}
        for row in moka_chr_all:
            moka_chr[row.Chr] = str(row.ChrID)
        return moka_chr

    def lookup_hgncid(self, gene_symbol):
        """
        Return HGNCID from Moka for a given gene symbol
        """
        sql = (
            "SELECT HGNCID "
              "FROM GenesHGNC_current "
             "WHERE ApprovedSymbol = '{gene}'"
        ).format(gene=gene_symbol)
        results = self.mc.fetchall(sql)
        if len(results) != 1:
            raise ValueError("No HGNCID found for {gene_symbol}".format(gene_symbol=gene_symbol))
        return results[0].HGNCID

    def chr_to_id(self, variant):
        """
        Converts chromosome names in variant dictionary to Moka chromosome IDs
        """
        # Convert chromosomes to Moka chromosome IDs
        # Do this for both build 37 and build 38 chromosome fields (or whichever one is present if liftover failed)
        if variant['chr37']:
            variant['chr37'] = self.moka_chr[variant['chr37']]
        if variant['chr38']:
            variant['chr38'] = self.moka_chr[variant['chr38']]
        return variant

    def already_in_moka(self, variant):
        """
        Checks if a variant has already been imported to Moka
        """
        if variant['pos37']:
            if ('GRCh37', variant['chr37'], variant['pos37'], variant['ref37'], variant['alt37']) in self.prev_vars:
                return True
        if variant['pos38']:
            if ('GRCh38', variant['chr38'], variant['pos38'], variant['ref38'], variant['alt38']) in self.prev_vars:
                return True
        return False

    def create_hgncid_lookup(self, variant):
        """
        Creates an HGNCID lookup dictionary for any genes that overlap with the variant
        """
        # Capture list of gene symbols
        genes = set([tx['gene'] for tx in variant['transcript_annotations'] if tx['gene']])
        # Lookup HGNCIDs for each gene symbol in Moka and record result in dictionary
        # If no HGNCID can be found in Moka, add this to no_hgncid list so it can be reported back to user at end
        hgncid_lookup = {}
        for gene in genes:
            try:
                hgncid = self.lookup_hgncid(gene)
            except ValueError:
                self.no_hgncid.append(gene)
            else:
                hgncid_lookup[gene] = hgncid
        return hgncid_lookup

    def empty_to_null(self, variant):
        """
        Converts empty strings to Nulls for sql
        """
        for field in variant.keys():
            if variant[field] == '':
                variant[field] = 'Null'
        return variant

    def wrap_strings(self, variant, fields):
        """
        Wraps strings in supplied fields in quotes for SQL
        """
        for field in fields:
            if variant[field] != 'Null':
                variant[field] = "'{}'".format(variant[field])
        return variant

    def insert_variant(self, variant):
        """
        Insert variant to Moka
        Args:
            variant = dictionary in the format output from get_additional_info() function
        Returns:
            status string, indicating outcome of function (e.g. imported, skipped)
        """
        # Insert to NGSVariant table
        sql = (
            "INSERT INTO NGSVariant (NGSTestID, InternalPatientID, DateAdded, PanelType, PanelTypeName, Gene, "
                   "ChrID, Position_hg19, ref, alt, genotype, "
                   "chr_id_38, position_38, ref_38, alt_38, genotype_38) "
                   "VALUES ({NGSTestID}, {InternalPatientID}, '{DateAdded}', 5, '100k', {Gene}, "
                   "{ChrID}, {Position_hg19}, {ref}, {alt}, {genotype}, "
                   "{chr_id_38}, {position_38}, {ref_38}, {alt_38}, {genotype_38})"
        ).format(
            NGSTestID=self.ngstest_id,
            InternalPatientID=self.internal_pat_id,
            DateAdded=datetime.datetime.now().strftime(r'%Y%m%d %H:%M:%S %p'),
            Gene=variant['concat_genes'],
            ChrID=variant['chr37'],
            Position_hg19=variant['pos37'],
            ref=variant['ref37'],
            alt=variant['alt37'],
            genotype=variant['gt37'],
            chr_id_38=variant['chr38'],
            position_38=variant['pos38'],
            ref_38=variant['ref38'],
            alt_38=variant['alt38'],
            genotype_38=variant['gt38']
        )
        self.mc.execute(sql)
        # Return the NGSVariantID for the record just inserted
        return self.mc.fetchone("SELECT @@IDENTITY")[0]

    def insert_transcript(self, ngs_variant_id, transcript):
        # Insert to NGSVariantAnnotations table
        sql = (
            "INSERT INTO NGSVariantAnnotations (NGSVariantID, HGNCID, GENE_SYMBOL, TRANSCRIPT_ID, HGVS_Transcript, HGVS_Protein) "
                    "VALUES ({NGSVariantID}, {HGNCID}, {GENE_SYMBOL}, {TRANSCRIPT_ID}, {HGVS_Transcript}, {HGVS_Protein})"
        ).format(
            NGSVariantID=ngs_variant_id,
            HGNCID=transcript['hgncid'],
            GENE_SYMBOL=transcript['gene'],
            TRANSCRIPT_ID=transcript['transcript'],
            HGVS_Transcript=transcript['hgvst'],
            HGVS_Protein=transcript['hgvsp']
        )
        self.mc.execute(sql)


def get_tier_1_2_vars(ir_id, proband_id):
    """
    Calls script on linux server to get tier 1 and 2 variants from cip-api
    Args:
        ir_id: interpretation request with version suffix, but without cip prefix (i.e. 12345-1 for 100KGP case SAP-12345-1)
        proband_id: GeL 100KGP participant ID for the proband
    Returns:
        list of variants. Each variant is stored in a dictionary containing basic variant details (build, chromosome, position, ref, alt, genotype)
    """
    # Call the script on server
    p = subprocess.Popen(
        (
            r'\\gstt.local\shared\Genetics_Data2\Array\Software\Python\python.exe '
            r'\\gstt.local\apps\Moka\Files\Software\genapp_ssh_n_run\ssh_n_run.py -c '
            r'"/home/mokaguys/miniconda2/envs/jellypy_py3/bin/python /home/mokaguys/Apps/100k_moka_variants/100k_tier_1_2_vars.py -i {ir_id} -p {proband_id}"'
        ).format(ir_id=ir_id, proband_id=proband_id),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    # Capture stdout and stderr
    stdout, stderr = p.communicate()
    # If there was an error, raise exception and exit
    if stderr or p.returncode != 0:
        raise Exception("Following error encountered when retrieving variants from GeL for {ir_id}:\n{stderr}".format(ir_id=ir_id, stderr=stderr))
    # Variants are printed to stdout in tsv format with each variant on a separate line
    # Parse the stdout and store variant details in list of dictionaries
    vars_split = stdout.rstrip().split('\n')
    var_details = []
    for var in vars_split:
        # If there's no tier 1 or 2 variants, var will just be an empty string so skip
        if var:
            fields = var.rstrip().split('\t')
            var_details.append(
                {
                    'build': fields[0],
                    'chr': fields[1],
                    'pos': fields[2],
                    'ref': fields[3],
                    'alt': fields[4],
                    'gt': fields[5]
                }
            )
    return var_details


def get_additional_info(variant):
    """
    Calls script on linux server to get additional variant annotations from VariantValidator API
    Args:
        variant: variant dictionary as output from get_tier_1_2_vars() containing fields 'build', 'chr', 'pos', 'ref', 'alt', 'gt'
    """
    # Call the script on server
    p = subprocess.Popen(
        (
            r'\\gstt.local\shared\Genetics_Data2\Array\Software\Python\python.exe '
            r'\\gstt.local\apps\Moka\Files\Software\genapp_ssh_n_run\ssh_n_run.py -c '
            r'"/home/mokaguys/miniconda2/envs/jellypy_py3/bin/python /home/mokaguys/Apps/100k_moka_variants/100k_variant_validator.py -b {build} -c {chr} -p {pos} -r {ref} -a {alt}"'
        ).format(build=variant['build'], chr=variant['chr'], pos=variant['pos'], ref=variant['ref'], alt=variant['alt']),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    # Capture stdout and stderr
    stdout, stderr = p.communicate()
    # If there was an error, raise exception and exit
    if stderr or p.returncode != 0:
        raise Exception("Following error encountered when retrieving annotations for variant:\n{stderr}".format(stderr=stderr))
    # Variants are printed to stdout in tsv format
    fields = stdout.rstrip().split('\t')
    # Genotype could be different for b38 and b37 (e.g. if ref and alt bases swapped in reference, could go from 1/1 to 0/0)
    # Therefore have separate fields for b37 and b38 GTs in Moka
    # Use the build to work out which field the gt should be stored in
    if variant['build'] == 'GRCh38':
        gt38 = variant['gt']
        gt37 = ''
    elif variant['build'] == 'GRCh37':
        gt38 = ''
        gt37 = variant['gt']
    # The transcript field contains details for each transcript returned from API
    # Store transcript annotations in list of dictionaries
    transcript_annotations = []
    if fields[8]:
        # Separate the transcripts by splitting on comma, then capture the individual transcript annotations by splitting on semi-colon
        for tx in fields[8].split(','):
            tx_fields = tx.split(';')
            transcript_annotations.append(
                {
                    'gene': tx_fields[0],
                    'transcript': tx_fields[1],
                    'hgvst': tx_fields[2],
                    'hgvsp': tx_fields[3]
                }
            )
    # Return the above along with remaining annotations in a dictionary
    return {
        'submitted_variant': variant,
        'chr37': fields[0],
        'pos37': fields[1],
        'ref37': fields[2],
        'alt37': fields[3],
        'chr38': fields[4],
        'pos38': fields[5],
        'ref38': fields[6],
        'alt38': fields[7],
        'transcript_annotations': transcript_annotations,
        'gt38': gt38,
        'gt37': gt37,
        'var_val_version': fields[9]
    }


def patient_log(internal_pat_id, ir_id, var_val_version, num_imported, num_failed, num_skipped, moka_connection):
    """
    Record in patient log that variants have been imported. Record number imported, failed and skipped, as well as the variant validator version.
    Args:
        internal_pat_id: Moka InternalPatientID from Patients table
        ir_id: interpretation request with version suffix, but without cip prefix (i.e. 12345-1 for 100KGP case SAP-12345-1)
        var_val_version: version of VariantValidator used to annotate (returned in get_additional_info() output)
        num_imported: Number of variants imported
        num_failed: Number of variants that failed import
        num_skipped: Number of variants that skipped import
        moka_connection: MokaConnector object
    """
    log_message = (
        "Imported 100KGP variants for case {ir_id} from interp-API. Annotated with VariantValidator {var_val_version}. "
        "{num_imported} imported. {num_failed} failed import. {num_skipped} skipped (already imported)."
    ).format(
        ir_id=ir_id,
        var_val_version=var_val_version,
        num_imported=num_imported,
        num_failed=num_failed,
        num_skipped=num_skipped
    )
    sql = (
        "INSERT INTO PatientLog (InternalPatientID, LogEntry, Date, Login, PCName) "
        "VALUES ({internalPatientID}, '{log_message}', '{today_date}', '{username}', '{computer}');"
        ).format(
            internalPatientID=internal_pat_id,
            log_message=log_message,
            today_date=datetime.datetime.now().strftime(r'%Y%m%d %H:%M:%S %p'),
            username=os.getenv('username'),
            computer=os.getenv('computername')
            )
    moka_connection.execute(sql)


def message_box(message, type):
    """
    Displays a message box
    """
    # Create root window, hide it (because messagebox is separate window), and bring to top
    root = Tkinter.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    # Display messagebox
    if type == 'info':
        tkMessageBox.showinfo("", message)
    elif type == 'warning':
        tkMessageBox.showwarning("", message)
    elif type == 'error':
        tkMessageBox.showerror("", message)
    else:
        raise ValueError("Type must be 'showinfo', 'showwarning' or 'error'")


def add_to_moka(variants_annotated, ngstest_id, internal_pat_id, moka_connection):
    """
    Performs steps required to add variants to Moka
    """
    imported = []
    skipped = []
    v = VariantAdder100KGP(ngstest_id, internal_pat_id, moka_connection)
    for variant in variants_annotated:
        # Convert variant chromosome name to a Moka chromosome ID
        variant = v.chr_to_id(variant)
        # Check if this variant has already been added to avoid duplication. If it has, skip.
        if v.already_in_moka(variant):
            skipped.append(variant['submitted_variant'])
        else:
            # Create an HGNCID lookup dictionary for each gene in the transcript annotations
            hgncid_lookup = v.create_hgncid_lookup(variant)
            # Add the concatenated gene list to the variant dictionary
            variant['concat_genes'] = ';'.join(sorted(hgncid_lookup.keys()))
            # Convert empty strings to null for SQL
            variant = v.empty_to_null(variant)
            # Wrap strings in quotes for SQL
            variant = v.wrap_strings(
                variant,
                fields=['ref37', 'alt37', 'ref38', 'alt38', 'gt38', 'gt37', 'concat_genes']
                )
            # Insert variant to Moka. Capture the NGSVariantID
            ngs_variant_id = v.insert_variant(variant)
            # Next insert transcript annotations to the NGSVariantAnnotations table
            for tx in variant['transcript_annotations']:
                # Lookup HGNCID, if it's not there skip to next transcript
                try:
                    tx['hgncid'] = hgncid_lookup[tx['gene']]
                except KeyError:
                    continue
                # Convert empty strings to null for SQL
                tx = v.empty_to_null(tx)
                # Fields that are strings need surrounding quotes in the SQL unless they are Null
                tx = v.wrap_strings(
                    tx,
                    fields=['gene', 'transcript', 'hgvst', 'hgvsp', 'hgncid']
                    )
                # Insert transcript to Moka
                v.insert_transcript(ngs_variant_id, tx)
            # Record that variant has been imported
            imported.append(variant['submitted_variant'])
    return imported, skipped, v.no_hgncid


def summary_messages(skipped, failed, no_hgncid):
    """
    Displays summary messages for user at end of script
    """
    if skipped:
        message_box(
            "Skipped import for the following variants, they are already in Moka.\n\n{variants}".format(
                    variants='\n'.join(['{build} chr{chr} {pos} {ref} {alt}'.format(**variant) for variant in skipped])
            ),
            'info'
        )
    if failed:
        message_box(
            "Import failed for following variants. Please add manually.\n\n{variants}".format(
                variants='\n'.join(['{build} chr{chr} {pos} {ref} {alt}'.format(**variant) for variant in failed])
            ),
            'warning'
        )
    if no_hgncid:
        message_box(
            "Couldn't find HGNCIDs for the following genes, so didn't import annotations.\n\n{genes}".format(
                genes='\n'.join([gene for gene in no_hgncid])
            ),
            'warning'
        )


def main():
    variants_annotated = []
    failed = []
    # Return arguments
    args = process_arguments()
    # Get tier 1/2 variants
    variants = get_tier_1_2_vars(args.ir_id, args.proband_id)
    # If there aren't any variants returned, display message and exit
    if not variants:
        message_box("No tier 1 or 2 variants found for 100KGP case {ir_id}".format(ir_id=args.ir_id), 'info')
        return
    # For each variant, get variant details. If anything goes wrong, record in failed list.
    for variant in variants:
        try:
            variants_annotated.append(get_additional_info(variant))
        except Exception:
            failed.append(variant)
    # Capture the version of variant validator used to annotate
    var_val_version = ''
    if variants_annotated:
        var_val_version = variants_annotated[0]['var_val_version']
    # Insert each annotated variant to Moka
    mc = MokaConnector()
    imported, skipped, no_hgncid = add_to_moka(variants_annotated, args.ngstest_id, args.internal_pat_id, mc)
    # Record in patient log.
    patient_log(args.internal_pat_id, args.ir_id, var_val_version, len(imported), len(failed), len(skipped), mc)
    # Report any skipped/failed imports
    summary_messages(skipped, failed, no_hgncid)
    #Get datetime stamp
    filename = datetime.datetime.now() 
    datetime_stamp = filename.strftime("%y%m%d_%H%M%S")
    #Save details of imported variants in logfile
    with open('/home/mokaguys/Documents/100k_moka_variants_logfiles/moka_import/%s_mokaimport.tsv' % datetime_stamp, 'w') as file_obj:
        file_obj.write(f"{imported}\t{skipped}\t{no_hgncid}\n")


if __name__ == '__main__':
    main()
