#!/bin/bash

## Define paths to software and reference files
BASEDIR=/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb #OK
SRATOOLKIT=/Applications/sratoolkit.2.7.0-mac64/bin #OK
REFERENCEDIR=$BASEDIR/hsapiens/reference_files #OK
RSEM=/Applications/RSEM-1.2.30/ #OK
NULLSIMULATION_NODE=$BASEDIR/hsapiens/no_diffexpression/null_simulation #OK
NONNULLSIMULATION_NODE=$BASEDIR/hsapiens/no_diffexpression/non_null_simulation #OK
NULLSIMULATION_DE=$BASEDIR/hsapiens/with_diffexpression/null_simulation #OK
NONNULLSIMULATION_DE=$BASEDIR/hsapiens/with_diffexpression/non_null_simulation #OK
ASTALAVISTA=/Applications/astalavista-3.2/bin/ #OK
FIGDIR=/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/hsapiens/figures #OK
RCODEGEN=/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/software/Rcode #OK
ROUT=/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/software/Rout #OK
DEXSEQ=/Users/koenvandenberge/Library/R/3.2/library/DEXSeq/python_scripts #OK
TOPHAT=/Applications/tophat-2.1.1.OSX_x86_64/ #OK
CUFFLINKS=/Applications/cufflinks-2.2.1.OSX_x86_64 #OK

## ------------------------- INPUT PREPARATION ----------------------------- ##
## The basis for the simulation is generated from a sample downloaded from the SRA.
## This sample was sequenced with an Illumina HiSeq, with a paired-end protocol,
## and with a read length of 101 bp.
## Extract the fastq files from the sra archive
$SRATOOLKIT/fastq-dump --split-files \
-Q 33 -O $REFERENCEDIR $REFERENCEDIR/SRR493366.sra

## Extract only lines corresponding to primary assembly chromosomes from gtf file
grep -v "^H" $REFERENCEDIR/Homo_sapiens.GRCh37.71.gtf \
> $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.gtf

## Extract only lines corresponding to protein coding genes from gtf file
grep "protein_coding" $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.gtf \
> $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf

## Prepare the reference files for RSEM
$RSEM/rsem-prepare-reference \
--gtf $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf --bowtie2 \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71

## Estimate the model file with RSEM
$RSEM/rsem-calculate-expression \
--paired-end --bowtie2 --seed 123 \
$REFERENCEDIR/SRR493366_1.fastq $REFERENCEDIR/SRR493366_2.fastq \
$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71 \
$REFERENCEDIR/rsem_model/SRR493366

## Plot some characteristics of the estimated model
$RSEM/rsem-plot-model $REFERENCEDIR/rsem_model/SRR493366 $FIGDIR/RSEM_model.pdf

## Modify the quality score distribution in the RSEM model file so that the probability
## of quality score 2 is 0. Otherwise, the quality scores of the simulated data may be very low.
## Also make sure that the transition probabilities into this state are 0.
## ->> $REFERENCEDIR/rsem_model/SRR493366.stat/SRR493366.highQ_adapted.model
## kvdb: wrote adaptModelFileRSEM_Human.R script for this

## Estimate and plot the isoform percentage distributions from the RSEM results
R CMD BATCH --no-restore --no-save "--args referencefile='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' outdir='$FIGDIR'" $RCODEGEN/isopct_distribution.R $ROUT/isopct_distribution_human.Rout

## Run ASTALAVISTA to classify splicing events
$ASTALAVISTA/astalavista -t asta \
-i $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf -e [ASE,ASI,DSP,VST]
gunzip $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf.gz

## Duplicate the protein_coding gtf file into one with the same name as the genome
## fasta file, for index building (not really necessary)
scp $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf

## Build TopHat index
bowtie2-build -f $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly
ln -s $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa

## Build TopHat transcriptome index
tophat -G $REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.gtf \
--transcriptome-index=$REFERENCEDIR/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
$REFERENCEDIR/TopHatIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly

## Build kallisto index from the fasta file created in the TopHat transcriptome index
kallisto index \
--index=$REFERENCEDIR/KallistoIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly \
$REFERENCEDIR/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa

## Create conversion table from transcript index to transcript name (to interpret kallisto output)
grep "^>" $REFERENCEDIR/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa | sed -e 's/>//' | cut -d" " -f1,2 > $REFERENCEDIR/KallistoIndex/TranscriptID_conversion.txt

## -------------------------- DATA SIMULATION ------------------------------ ##
## Estimate the simulation parameters for the individual samples

R CMD BATCH --no-restore --no-save "--args path_to_generate_rsem_files='$RCODEGEN/kvdb/generate_rsem_files_function_kvdb_flipMostExpressedTranscripts_v2_fiveSamples_alsoLowExpressedGenes_human.R' seed=123 isoform_results_file='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' nbr_per_group=5 meandisp.file='$REFERENCEDIR/Pickrell.Cheung.Mu.Phi.Estimates.rds' outdirbase='/Volumes/HDKoen2/data/dtu/diff_splice_paper_Kvdb/hsapiens/diffexpression/non_null_simulation_dge' librarysize=40000000 keepchr=NULL nbr_diff_spliced=1000 nbr_diff_expr=1000 fold_changes='expon'" $RCODEGEN/kvdb/generate_rsem_files_human_run_kvdb.R $ROUT/generate_rsem_files_human_run_de.Rout

R CMD BATCH --no-restore --no-save "--args path_to_sim_details='$NONNULLSIMULATION_DE/3_truth/simulation_details.txt' output_pdf='$FIGDIR/simulation_details_withdiffexp_nonnull.pdf'" $RCODEGEN/plot_simulation_details.R $ROUT/plot_simulation_details_withde_nonnull.Rout

## Generate truth files
R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$NULLSIMULATION_DE/3_truth/simulation_details.txt' out.file='$NULLSIMULATION_DE/3_truth/truth_human_null.txt' astalavista.file=NULL gtf.file=NULL flattened.gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' missing.annot.file=NULL" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_human_null_de.Rout

R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$NONNULLSIMULATION_DE/3_truth/simulation_details.txt' out.file='$NONNULLSIMULATION_DE/3_truth/truth_human_non_null.txt' astalavista.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf' gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' flattened.gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff' missing.annot.file=NULL" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_human_nonnull_de.Rout

## Simulate reads for 10 samples (40 million pairs/sample),

for n in 1 2 3 4 5 6 7 8 9 10
do
$RSEM/rsem-simulate-reads \
$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71 \
$REFERENCEDIR/rsem_model/SRR493366.stat/SRR493366.highQ_adapted.model \
/Volumes/HDKoen2/data/dtu/diff_splice_paper_Kvdb/hsapiens/diffexpression/non_null_simulation_dge/non_null_simulation/1_reads/rsem_files/sample${n}.txt \
0.05 40000000 /Volumes/HDKoen2/data/dtu/diff_splice_paper_Kvdb/hsapiens/diffexpression/non_null_simulation_dge/non_null_simulation/1_reads/reads/sample${n}/sample_${n} \
--seed 123
done

