For use with full length 16S microbiome data, raw output from a 16S sequencing run using the Sequel IIe Platform from PacBio. You will need to clone the mcsmrt repository as well as the report script located in this repository. Additionally you must have Pacbio SMRT Link software running version 11.


# Pipeline and report

Please note the standardardized file-naming convention that must be used in order to run this program.

Organization is as follows: 
project_dir/SMRTcell_subdir/CCS_or_DEMUX_subdir/outputs

example: 

m54257_12345_20211108/1_A01/ccs/outputs/

another example: 

m54257_67890_202111/2_B01/demux_no_peek/outputs/

It may be possible to edit file names and directory structure as per your own needs, but this is the convention we use.

# Required Files:

CCS Sub-report:
1) ccs_passes.tsv --> 
a tsv file with recorded number of passess associated with each read (required for MCSMRT pipeline, script can be found on MCSMRT repository)
2) ccs_report.zip OR ccs_statistics.csv --> 
a csv file output from pbcromwell pb_ccs 

DEMUX Sub-report:
1) barcode_to_sample.csv --> 
a csv file with a translation of barcode ID and bio sample ID. (barcode,sample)
2) barcode_ccs_summary.csv --> 
a csv file output from pbcromwell pb_demux_ccs

note that the script will check to make sure these files exist.
