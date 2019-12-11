# 100k_moka_variants v1.1
## 100k_tier_1_2_vars.py 

This script pulls out all tier 1 and 2 variants from the CIP-API for a 100k case.

The variant details are printed to stdout in a tab separated format with the following columns:

> Build&nbsp;&nbsp;&nbsp;&nbsp;Chromosome&nbsp;&nbsp;&nbsp;&nbsp;Position&nbsp;&nbsp;&nbsp;&nbsp;Ref&nbsp;&nbsp;&nbsp;&nbsp;Alt&nbsp;&nbsp;&nbsp;&nbsp;Genotype

e.g.

> GRCh38&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;132689263&nbsp;&nbsp;&nbsp;&nbsp;CCT&nbsp;&nbsp;&nbsp;&nbsp;C&nbsp;&nbsp;&nbsp;&nbsp;0/1<br>GRCh38&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;132682750&nbsp;&nbsp;&nbsp;&nbsp;C&nbsp;&nbsp;&nbsp;&nbsp;CCGA&nbsp;&nbsp;&nbsp;&nbsp;0/1


### Usage

This script requires access to the CIPAPI so must be run on our trust linux server.

Requirements:

* Python 3.7+
* Access to CIPAPI
* JellyPy (in PYTHONPATH)
* GelReportModels (v6 or higher)


On `SV-TE-GENAPP01` activate the `jellypy_py3` conda environment so that above requirements are met:

```
source activate jellypy_py3
```

Run the script:

```
python /home/mokaguys/Apps/100k_moka_variants/100k_tier_1_2_vars.py -i IR_ID -p PROBAND_ID
```
Arguments:
```
  -i IR_ID, --ir_id IR_ID
                        GeL Interpretation Request ID in format 12345-1
  -p PROBAND_ID, --proband_id PROBAND_ID
                        GeL participant ID for proband
```

## 100k_variant_validator.py

This script uses the variant details output from `100k_tier_1_2_vars.py` to query the [VariantValidator](https://variantvalidator.org/) API and returns variant details formatted for entry to Moka.

The variant details are printed to stdout in a tab separated format with the following columns:

> Chr37&nbsp;&nbsp;&nbsp;&nbsp;Pos37&nbsp;&nbsp;&nbsp;&nbsp;ref37&nbsp;&nbsp;&nbsp;&nbsp;Alt37&nbsp;&nbsp;&nbsp;&nbsp;Chr38&nbsp;&nbsp;&nbsp;&nbsp;Pos38&nbsp;&nbsp;&nbsp;&nbsp;ref38&nbsp;&nbsp;&nbsp;&nbsp;Alt38&nbsp;&nbsp;&nbsp;&nbsp;Transcript_Details&nbsp;&nbsp;&nbsp;&nbsp;Variant_Validator_Version


The Transcript_Details field is a string containing details for each transcript returned from the API. 

For each transcript this string contains the gene, transcript accession, transcript hgvs and protein hgvs separated by semi-colons. i.e.

> gene;transcript_accession;hgvs_trsncript;hgvs_protein

Multiple transcripts are separated by commas. So the final string is in the format:

> gene1;transcript_accession1;hgvs_t1;hgvs_p1,gene2;transcript_accession2;hgvs_t2;hgvs_p2 ...etc

An example output would be:

> 3&nbsp;&nbsp;&nbsp;&nbsp;132401594&nbsp;&nbsp;&nbsp;&nbsp;C&nbsp;&nbsp;&nbsp;&nbsp;CCGA&nbsp;&nbsp;&nbsp;&nbsp;3&nbsp;&nbsp;&nbsp;&nbsp;132682750&nbsp;&nbsp;&nbsp;&nbsp;C&nbsp;&nbsp;&nbsp;&nbsp;CCGA&nbsp;&nbsp;&nbsp;&nbsp;NPHP3;NM_153240.4;c.3762_3764dup;p.R1255dup,NPHP3-ACAD11;NR_037804.1;n.3768_3770dup;n.&nbsp;&nbsp;&nbsp;&nbsp;v0.2.4_post6


### Usage

This script requires access to the VariantValidator API so must be run on our linux server.

Requirements:

* Python 3.7+
* requests

On `SV-TE-GENAPP01` activate the `jellypy_py3` conda environment so that above requirements are met:

```
source activate jellypy_py3
```

Run the script:

```
python /home/mokaguys/Apps/100k_moka_variants/100k_variant_validator.py -b BUILD -c CHR -p POS -r REF -a ALT
```
Arguments:
```
  -b BUILD, --build BUILD
                        Genome build
  -c CHR, --chr CHR     Chromosome
  -p POS, --pos POS     Position
  -r REF, --ref REF     Reference base(s)
  -a ALT, --alt ALT     Alternate base(s)
```

## 100k_variant_import.py 

This script acts as a wrapper for `100k_tier_1_2_vars.py ` and `100k_variant_validator.py`, and inserts the variants into the Moka database.

It is designed to be called from the Moka UI.

### Usage

This script requires access to Moka so must be run from a trust Windows desktop with access to `\\gstt.local\apps\Moka\Files\Software\`. It uses the [`ssh_n_run.py`](https://github.com/moka-guys/ssh_to_genapp) script to run code remotely on the trust linux server.

Requirements:

* Python 2.7
* paramiko
* pyodbc

The connection details for Moka are in `config.ini` in the same directory of the script. See `example_config.ini` for format.

Run the script:

```
python 100k_variant_import.py [-h] --ir_id IR_ID --proband_id PROBAND_ID
                              --ngstest_id NGSTEST_ID --internal_pat_id
                              INTERNAL_PAT_ID
```

Arguments:
```
  --ir_id IR_ID         GeL Interpretation Request ID in format 12345-1
  --proband_id PROBAND_ID
                        GeL participant ID for proband
  --ngstest_id NGSTEST_ID
                        Moka NGSTestID
  --internal_pat_id INTERNAL_PAT_ID
                        Moka InternalPatientID
```
