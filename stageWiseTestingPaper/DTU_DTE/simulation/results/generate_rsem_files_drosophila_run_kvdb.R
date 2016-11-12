## ----- generate_rsem_files_drosophila_run
## <<generate_rsem_files_drosophila_run.R>>


#refdir
#/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files

#path_to_generate_rsem_files="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/software/Rcode/generate_rsem_files_function.R"
#seed=123
#isoform_results_file="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files/rsem_model/SRR1501444.isoforms.results"
#nbr_per_group=3
#meandisp.file="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files/Pickrell.Cheung.Mu.Phi.Estimates.rds"
#outdirbase="/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/no_diffexpression"
#librarysize=25000000
#keepchr=NULL
#nbr_diff_spliced=1000
#nbr_diff_expr=0
#fold_changes=NULL


library(biomaRt)

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_generate_rsem_files)
print(seed)
print(isoform_results_file)
print(nbr_per_group)
print(meandisp.file)
print(outdirbase)
print(librarysize)
print(keepchr)
print(nbr_diff_spliced)
print(nbr_diff_expr)

## Generate fold changes if there are differentially expressed genes
if (!is.null(fold_changes) && fold_changes == "expon") {
  set.seed(seed)
  fold_changes <- (2 + rexp(nbr_diff_expr, 
                            rate = 1))^(c(-1, 1)[round(runif(nbr_diff_expr)) + 1])
}

print(head(fold_changes))

## Read file with mean/dispersion relationship
meandisp <- readRDS(meandisp.file)
meanvect <- meandisp$pickrell.cheung.mu
dispvect <- meandisp$pickrell.cheung.phi

## Use the same ENSEMBL version as our gtf/reference to get chromosomes for genes
#ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host = "jan2013.archive.ensembl.org")
# this is the closest version according to http://www.ensembl.org/info/website/archives/index.html, which I got to work. jan2013 gave errors.
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host = "may2012.archive.ensembl.org") 


## Generate RSEM files for the individual samples
source(path_to_generate_rsem_files)
generate_rsem_files(seed = seed, isoform_results_file = isoform_results_file,
                    nbr_per_group = nbr_per_group, meanvect = meanvect, 
                    dispvect = dispvect, outdirbase = outdirbase, 
                    librarysize = librarysize, keepchr = keepchr, 
                    nbr_diff_spliced = nbr_diff_spliced, 
                    nbr_diff_expr = nbr_diff_expr, fold_changes = fold_changes,
                    mart = ensembl)


