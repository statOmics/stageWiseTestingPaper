#!/bin/bash

## Define paths to software and reference files
BASEDIR=/Volumes/HDKoen2/data/dtu/diff_splice_paper_Kvdb #OK
REFERENCEDIR=/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files #OK
RSEM=/Applications/RSEM-1.2.30/ #OK
NULLSIMULATION_NODE=/Volumes/HDKoen/data/dtu/diff_splice_paper_Kvdb/drosophila/diffexpression/null_simulation #OK
NONNULLSIMULATION_NODE=/Volumes/HDKoen/data/dtu/diff_splice_paper_Kvdb/drosophila/diffexpression/non_null_simulation_dge #OK
ASTALAVISTA=/Applications/astalavista-3.2/bin/ #OK
FIGDIR=/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/figures #OK
RCODEGEN=/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/software/Rcode #OK
ROUT=/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/software/Rout #OK
DEXSEQ=/Users/koenvandenberge/Library/R/3.2/library/DEXSeq/python_scripts #OK
TOPHAT=/Applications/tophat-2.1.1.OSX_x86_64/ #OK
CUFFLINKS=/Applications/cufflinks-2.2.1.OSX_x86_64 #OK
KALLISTO=/Applications/kallisto/ #OK

## ------------------------- INPUT PREPARATION ----------------------------- ##

## The basis for the simulation is generated from a sample downloaded from ENA.
## This sample was sequenced with an Illumina HiSeq, with a paired-end protocol,
## and with a read length of 101 bp.

## Extract only lines corresponding to protein coding genes from gtf file
grep "protein_coding" $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.gtf \
> $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf

## Prepare the reference files for RSEM
$RSEM/rsem-prepare-reference \
--gtf $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf --bowtie2 \
$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa \
$REFERENCEDIR/rsem_reference/Drosophila_melanogaster.BDGP5.70

## Estimate the model file with RSEM
$RSEM/rsem-calculate-expression \
--paired-end --bowtie2 --seed 123 \
$REFERENCEDIR/SRR1501444_1.fastq $REFERENCEDIR/SRR1501444_2.fastq \
$REFERENCEDIR/rsem_reference/Drosophila_melanogaster.BDGP5.70 \
$REFERENCEDIR/rsem_model/SRR1501444


## Plot some characteristics of the estimated model
$RSEM/rsem-plot-model $REFERENCEDIR/rsem_model/SRR1501444 $FIGDIR/RSEM_model.pdf

## Modify the quality score distribution in the RSEM model file so that the probability
## of quality score 2 is 0. Otherwise, the quality scores of the simulated data may be very low.
## Also make sure that the transition probabilities into this state are 0.
## ->> $REFERENCEDIR/rsem_model/SRR1501444.stat/SRR1501444.highQ_soneson.model
## this can alternatively be done with the adaptModelFileRSEM_Drosophila.R script.

## Estimate and plot the isoform percentage distributions from the RSEM results
R CMD BATCH --no-restore --no-save "--args referencefile='$REFERENCEDIR/rsem_model/SRR1501444.isoforms.results' outdir='$FIGDIR'" $RCODEGEN/isopct_distribution.R $ROUT/isopct_distribution_drosophila.Rout

## Run ASTALAVISTA to classify splicing events
$ASTALAVISTA/astalavista -t asta \
-i $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf -e [ASE,ASI,DSP,VST]
gunzip $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding_sorted.gtf_astalavista.gtf.gz

## Duplicate the protein_coding gtf file into one with the same name as the genome
## fasta file, for index building (not really necessary)
scp $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf \
$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.dna.toplevel.gtf

## Build TopHat index
## make TopHatIndex folder before running
bowtie2-build -f $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa \
$REFERENCEDIR/TopHatIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel
ln -s $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa  \
$REFERENCEDIR/TopHatIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa

## Build TopHat transcriptome index
tophat -G $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.dna.toplevel.gtf \
--transcriptome-index=$REFERENCEDIR/TopHatTranscriptomeIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel \
$REFERENCEDIR/TopHatIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel

## Build kallisto index from the fasta file created in the TopHat transcriptome index
##make KallistoIndex folder
$KALLISTO/kallisto index \
--index=$REFERENCEDIR/KallistoIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel \
$REFERENCEDIR/TopHatTranscriptomeIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa

## Create conversion table from transcript index to transcript name (to interpret kallisto output)
grep "^>" $REFERENCEDIR/TopHatTranscriptomeIndex/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa | sed -e 's/>//' | cut -d" " -f1,2 > $REFERENCEDIR/KallistoIndex/TranscriptID_conversion.txt

## Prepare flattened annotations (for DEXSeq)
python $DEXSEQ/dexseq_prepare_annotation.py \
$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf \
$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.flattened.gff

python $DEXSEQ/dexseq_prepare_annotation.py --aggregate='no' \
$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf \
$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.flattened.nomerge.gff

## Manually create flattened gtf/gff file (for use with featureCounts and DEXSeq)
#R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf' output_gtf='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.flatman.ign.gtf' ignore_strand=TRUE" $RCODEGEN/generate_flattened_gtf.R $ROUT/generate_flattened_gtf.Rout
#sed -i -e 's/[*]/./' $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.flatman.ign.gff

## Manually create flattened gtf/gff file (for use with featureCounts and DEXSeq) - strand specific overlap
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf' output_gtf='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.flatman.strandspec.gtf' ignore_strand=FALSE" $RCODEGEN/generate_flattened_gtf.R $ROUT/generate_flattened_gtf_strandspec.Rout
sed -i -e 's/[*]/./' $REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.flatman.strandspec.gff

## Fix original gtf file to count on (real) exon level with featureCounts
R CMD BATCH --no-restore --no-save "--args input_gtf='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf' output_gtf='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.exongene.gtf'" $RCODEGEN/generate_renamed_gtf_for_featurecounts.R $ROUT/generate_renamed_gtf_for_featurecounts.Rout


## -------------------------- DATA SIMULATION ------------------------------ ##
## Estimate the simulation parameters for the individual samples
## note that you need internet connection for this command as it calls biomaRt.

## switching order for some number of most expresed transcripts within a gene
## changed to five samples per group
R CMD BATCH --no-restore --no-save "--args path_to_generate_rsem_files='$RCODEGEN/kvdb/generate_rsem_files_function_kvdb_flipMostExpressedTranscripts_v2_fiveSamples_alsoLowExpressedGenes.R' seed=123 isoform_results_file='$REFERENCEDIR/rsem_model/SRR1501444.isoforms.results' nbr_per_group=5 meandisp.file='$REFERENCEDIR/Pickrell.Cheung.Mu.Phi.Estimates.rds' outdirbase='$BASEDIR/drosophila/diffexpression/non_null_simulation_dge' librarysize=25000000 keepchr=NULL nbr_diff_spliced=1000 nbr_diff_expr=1000 fold_changes='expon'" $RCODEGEN/kvdb/generate_rsem_files_drosophila_run_kvdb.R $ROUT/generate_rsem_files_drosophila_run_node_kvdb.Rout

R CMD BATCH --no-restore --no-save "--args path_to_sim_details='$NONNULLSIMULATION_NODE/3_truth/simulation_details.txt' output_pdf='$FIGDIR/simulation_details_nodiffexp_nonnull.pdf'" $RCODEGEN/plot_simulation_details_2.R $ROUT/plot_simulation_details_2_node_nonnull_drosophila_kvdb.Rout

#truth file## Generate truth files
R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$NONNULLSIMULATION_NODE_FLIP/3_truth/simulation_details.txt' out.file='$NONNULLSIMULATION_NODE_FLIP/3_truth/truth_drosophila_non_null.txt' astalavista.file='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding_sorted.gtf_astalavista.gtf' gtf.file='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf' flattened.gtf.file='$REFERENCEDIR/Drosophila_melanogaster.BDGP5.70.protein_coding.flattened.nomerge.gff' missing.annot.file=NULL" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_drosophila_nonnull_node.Rout

## Simulate reads for 10 samples (25 million pairs/sample),

for n in 1 2 3 4 5 6 7 8 9 10
do
$RSEM/rsem-simulate-reads \
$REFERENCEDIR/rsem_reference/Drosophila_melanogaster.BDGP5.70 \
$REFERENCEDIR/rsem_model/SRR1501444.stat/SRR1501444.highQ_soneson.model \
$BASEDIR/drosophila/diffexpression/non_null_simulation_dge/non_null_simulation/1_reads/rsem_files/sample${n}.txt \
0.05 25000000 $BASEDIR/drosophila/diffexpression/non_null_simulation_dge/non_null_simulation/1_reads/reads/sample${n}/sample_${n} \
--seed 123
done



