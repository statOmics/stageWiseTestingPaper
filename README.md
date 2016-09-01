# stageWiseTestingPaper

This repository contains all required code to reproduce the analyses for the stage-wise testing in RNA-Seq paper. The repository consists of two major folders: DGE (differential gene expression) and DTU (differential transcript usage), with the code to reproduce analyses from the respective applications. Both repositories contain a subfolder simulation and caseStudy.

## Differential gene expression

The simulations for the DGE analysis were based on the framework provided by the edgeR robust paper:
Xiaobei Zhou, Helen Lindsay and Mark D. Robinson. Robustly detecting differential expression in RNA sequencing data using observation weights. [Nucleic Acids Research 42 (11): e91](http://nar.oxfordjournals.org/content/42/11/e91)

The original code for the differential gene expression simulation study can be found on their project’s website:
http://imlspenticton.uzh.ch/robinson_lab/edgeR_robust/


In the DGE case study, we analysed the Hammer dataset. We downloaded the raw count table from the [ReCount](http://biorxiv.org/content/early/2016/08/08/068478) project website:
http://bowtie-bio.sourceforge.net/recount/


## Differential transcript usage