#! /bin/env bash

# This script takes three parameters, the pacbio sequel IIe runID, the cell ID, and the collaborator last name
# For the Demultiplex json the expectation is that it exists in the required_ccs_files
# directory.  If it does not, link it or copy it there. You MUST update all required CCS and DEMUX files prior to
# starting this pipeline. CCS input file should stay the same, but you need to add the location of the barcode2sample CSV
# file to the DEMUX input file.

#example command:
#./ccs_demux_v9.sh r54257_20200728_141738 1_A01 czerski


todaysDate=$(date +%Y%m%d)
runID=$1
cell=$2
nickname=$3
run_location=/data/pacbio/sequel_IIe/userdata/data_root/$runID
cell_location=/data/pacbio/sequel_IIe/userdata/data_root/$runID/$cell

echo "This is the run location:"
echo $run_location
echo "This is the cell location:"
echo $cell_location

echo "Making project directory..."
#Check if the directory exists, if not make it
[[ -d $runID/$cell ]] || mkdir -p $runID/$cell

echo "Linking required files..."
#Check if there are already links to the data, and required files. Create if not.
[[ -L $runID/$cell/$cell ]] || ln -s /data/pacbio/sequel_IIe/userdata/data_root/$runID/$cell $runID/$cell
[[ -L $runID/$cell/pb_cromwell_files ]] || ln -s /data/shared/homes/sam/pb_cromwell_files $runID/$cell
[[ -L $runID/$cell/"$runID"_"$cell".csv ]] || ln -s /data/shared/homes/sam/16S_projects_sequel/$nickname/"$runID"_"$cell"/"$runID"_"$cell".csv $runID/$cell

#Make Subreadset Variable
subreadset=$runID/$cell/$cell/m64328e_*.subreadset.xml

echo "Executing CCS analysis..."
#Execute ccs analysis
pbcromwell run pb_ccs \
	   -e $subreadset \
	   --config $runID/$cell/pb_cromwell_files/2023_cromwell.conf \
	   --output-dir $runID/$cell/ccs \
	   #--overwrite \ #I don't need to run CCS twice.
	   -i $runID/$cell/pb_cromwell_files/2023_inputs_ccs.json

echo "Demultiplexing..."
#Execute demultiplexing
pbcromwell run pb_demux_ccs \
	   -e $runID/$cell/ccs/outputs/final.consensusreadset.xml \
	   -e /scratch/sam/required_ccs_files/all_sequel_bc_v1v2.xml \
	   --config $runID/$cell/pb_cromwell_files/2023_cromwell.conf \
	   --output-dir $runID/$cell/demux_no_peek \
	   #--overwrite \
	   -i $runID/$cell/pb_cromwell_files/2023_inputs_demux.json

echo "Preparing for microbiome classification..."
#Execute scripts required for MCSMRT
#1) Unzip FASTQ files
gunzip -f $runID/$cell/demux_no_peek/outputs/*.gz

#2) Rename FASTQ with Sample ID
sh /data/shared/homes/sam/apps/file_manipulation/rename_fqs_from_csv.sh \
    -f $runID/$cell/demux_no_peek/outputs/ \
    -c $runID/$cell/"$runID"_"$cell".csv \
    -o $runID/$cell/renamed_fqs

#3) Generate CCS Passes TSV
ln -s /data/shared/homes/sam/apps/file_manipulation/bam_to_passes_tsv.sh $runID/$cell/ccs/outputs/
sh $runID/$cell/ccs/outputs/bam_to_passes_tsv.sh $runID $cell

#4) Add Passes to Renamed FASTQs
ruby mcsmrt/add_samp_and_ccs_passes_qserr.rb \
     -f $runID/$cell/renamed_fqs \
     -p $runID/$cell/ccs/outputs/ccs_passes.tsv \
     -o $runID/$cell/renamed_fqs_w_passes

echo "Executing MCSMRT..."
#ln -s /data/shared/homes/sam/apps/mcsmrt
ruby mcsmrt/mcsmrt.rb -d 32 \
     -f $runID/$cell/renamed_fqs_w_passes \
     -c mcsmrt/all_required_mcsmrt_files/rdp_gold.fa \
     -t mcsmrt/all_required_mcsmrt_files/16sMicrobial_ncbi_lineage_reference_database.udb \
     -l mcsmrt/all_required_mcsmrt_files/16sMicrobial_ncbi_lineage_reference.fasta \
     -g mcsmrt/all_required_mcsmrt_files/human_genome_for_mcsmrt/human_g1k_v37.fasta \
     -p mcsmrt/all_required_mcsmrt_files/primers.fasta \
     -b mcsmrt/all_required_mcsmrt_files/ncbi_clustered_table.tsv \
     -o yes \
     -j after \
     -v 

echo "Organizing..."
#make dir for mcmsrt output
mcsmrt_outdir="$nickname"_"$todaysDate"
mkdir $mcsmrt_outdir

#move output to new dir
mv $runID/$cell/pre_* $mcsmrt_outdir/.
mv $runID/$cell/post_* $mcsmrt_outdir/.
mv $runID/$cell/split_otus $mcsmrt_outdir/.
mv $runID/$cell/all_bc_tax_info* $mcsmrt_outdir/.

#Make a new folder in 16S_projects_sequel if it doesn't exist
proj_dirs=/data/shared/homes/sam/16S_projects_sequel
if [[ -d $proj_dirs/$nickname ]] ; then
    mkdir $proj_dirs/$nickname/"$runID"_"$cell"
else
    mkdir $proj_dirs/$nickname ; mkdir $proj_dir/$nickname/"$run_ID"_"$cell" 
fi

#Make dir for report output
mkdir $runID/$cell/sequencing_analysis

echo "Executing report generator script..."
python3 /data/shared/homes/sam/apps/ccs_demux/ccs_and_demux_report.py -c $runID/$cell -o $runID/$cell/sequencing_analysis

echo 'See "${proj_dirs)"/"${nickname}"/"${runID}"_"${cell}" for all results and raw data'
ln -s -t $proj_dirs/$nickname/"$runID"_"$cell"/ $runID/$cell
mv $proj_dirs/$nickname/"$runID"_"$cell"/$cell/ $proj_dirs/$nickname/$runID/$cell/ccs_and_demultiplex_output

ln -s -t $proj_dirs/$nickname/"$runID"_"$cell"/ $runID/$cell/renamed_fqs

ln -s -t $proj_dirs/$nickname/"$runID"_"$cell"/ $runID/$cell/$mcsmrt_output
mv $proj_dirs/$nickname/"$runID"_"$cell"/$mcsmrt_output $proj_dirs/$nickname/"$runID"_"$cell"/mcsmrt_output

ln -s -t $proj_dirs/$nickname/"$runID"_"$cell"/ $runID/$cell/sequencing_analysis

echo 'Typical outputs to be emailed to the group are be found here: "$proj_dirs"/"$nickname"/"$runID"_"$cell"/results'
mkdir $proj_dirs/$nickname/$runID_$cell/results
res_dir=$proj_dirs/$nickname/$runID_$cell/results

ln -s -t $res_dir $proj_dirs/$nickname/$runID_$cell/sequencing_analysis/barcodes_plots.pdf
ln -s -t $res_dir $proj_dirs/$nickname/$runID_$cell/sequencing_analysis/ccs_plots.pdf
ln -s -t $res_dir $proj_dirs/$nickname/$runID_$cell/sequencing_analysis/barcodes_summary.pdf
ln -s -t $res_dir $proj_dirs/$nickname/$runID_$cell/mcsmrt_output/post_final_results.tsv

echo "fin"
