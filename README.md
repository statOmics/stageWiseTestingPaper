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

The simulations for DTU analysis are based upon the framework provided by
C Soneson\*, KL Matthes\*, M Nowicka, CW Law & MD Robinson: Differential transcript usage from RNA-seq data: isoform pre-filtering improves performance of count-based methods. [Genome Biology 17:12 (2016)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0862-3).

The simulation study requires some files to be downloaded a priori and we therefore suggest to consult their GitHub repository (https://github.com/markrobinsonuzh/diff_splice_paper).
Our repository is similar to their repository, with some adaptations to the bash scripts and R scripts. Details on the changes are reported in the paper.


The DTU case study considers the dataset provided by

Shancheng Ren\*, Zhiyu Peng\*, Jian-Hua Mao\*, Yongwei Yu, Changjun Yin, Xin Gao, Zilian Cui, Jibin Zhang, Kang Yi, Weidong Xu, Chao Chen, Fubo Wang, Xinwu Guo, Ji Lu, Jun Yang, Min Wei, Zhijian Tian, Yinghui Guan, Liang Tang, Chuanliang Xu, Linhui Wang, Xu Gao, Wei Tian, Jian Wang, Huanming Yang, Jun Wang and Yinghao Sun. RNA-seq analysis of prostate cancer in the Chinese population identifies recurrent gene fusions, cancer-associated long noncoding RNAs and aberrant alternative splicings. [Cell Research (2012) 22:806–821.](http://www.nature.com/cr/journal/v22/n5/full/cr201230a.html)

We downloaded the raw, unnormalized kallisto processed data from the [Bear's lair](http://biorxiv.org/content/early/2016/05/31/056200) project website (http://pachterlab.github.io/lair/).

Please contact koen.vandenberge@ugent.be or raise an issue on GitHub for questions.

git commit x
git push stageWiseRep master
