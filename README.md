# stageWiseTestingPaper

This repository contains all required code to reproduce the analyses for the stage-wise testing in RNA-Seq paper. The repository consists of two major folders: DGE (differential gene expression) and DTU_DTE (differential transcript usage, differential transcript expression), with the code to reproduce analyses from the respective applications. Both repositories contain a subfolder simulation and caseStudy.

## Differential gene expression

The simulations for the DGE analysis were based on the framework provided by the edgeR robust paper:
X Zhou, H Lindsay and MD Robinson. Robustly detecting differential expression in RNA sequencing data using observation weights. [Nucleic Acids Research 42 (11): e91](http://nar.oxfordjournals.org/content/42/11/e91)

The original code for the differential gene expression simulation study can be found on their project’s website:
http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/


In the DGE case study, we analysed the Hammer dataset. We downloaded the raw count table from the [ReCount](http://biorxiv.org/content/early/2016/08/08/068478) project website:
http://bowtie-bio.sourceforge.net/recount/


## Differential transcript usage

In the introduction, we analyse the sim2_human simulated data from
C Soneson, MI Love, MD Robinson: Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences [version 2; referees: 2 approved]. [F1000Research 2016, 4:1521](https://f1000research.com/articles/4-1521/v2).
Details and required files for this simulation study can be found at [the supplementary data of that paper](https://f1000researchdata.s3.amazonaws.com/datasets/7563/315e2602-541f-4781-ab6e-76635dab0360_Sim_2_Quantification.html).

All other simulations for DTU or DTE analysis are based upon the framework provided by
C Soneson\*, KL Matthes\*, M Nowicka, CW Law & MD Robinson: Differential transcript usage from RNA-seq data: isoform pre-filtering improves performance of count-based methods. [Genome Biology 17:12 (2016)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0862-3).

The simulation study requires some files to be downloaded a priori and we therefore suggest to consult their GitHub repository (https://github.com/markrobinsonuzh/diff_splice_paper).
Our repository is similar to their repository, with some adaptations to the bash scripts and R scripts. Details on the changes are reported in the paper.

The DTU case study considers the dataset provided by

S Ren\*, Z Peng\*, J Mao\*, Y Yu, C Yin, X Gao, Z Cui, J Zhang, K Yi, W Xu, C Chen, F Wang, X Guo, J Lu, J Yang, M Wei, Z Tian, Y Guan, L Tang, C Xu, L Wang, X Gao, W Tian, J Wang, H Yang, J Wang and Y Sun. RNA-seq analysis of prostate cancer in the Chinese population identifies recurrent gene fusions, cancer-associated long noncoding RNAs and aberrant alternative splicings. [Cell Research (2012) 22:806–821.](http://www.nature.com/cr/journal/v22/n5/full/cr201230a.html)

We downloaded the raw, unnormalized kallisto processed data from the [Bear's lair](http://biorxiv.org/content/early/2016/05/31/056200) project website (http://pachterlab.github.io/lair/).

Please contact koen.vandenberge@ugent.be or raise an issue on GitHub for questions.


