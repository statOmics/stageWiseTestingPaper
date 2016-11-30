setwd("/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/drosophila/reference_files/rsem_model/SRR1501444.stat/")
highQModel=readLines("SRR1501444.model")
probMatrix=highQModel[14:114]
probMatList <- sapply(1:length(probMatrix),function(i) as.numeric(strsplit(probMatrix[i],split=" ")[[1]]), simplify=FALSE)
probMatrixAdapted=vector()
for(i in 1:length(probMatList)){ 
     #set initial/trans prob of quality score 2 to zero
    probMatList[[i]][3]=0
    #rescale to sum to 1
    if(sum(probMatList[[i]])==0){probMatList[[i]]=probMatList[[i]]}else{probMatList[[i]] <- probMatList[[i]]/sum(probMatList[[i]])}
    #retransform in original format for file writing
    probMatrixAdapted[i] <- Reduce(paste,as.character(probMatList[[i]]))
}
highQModel[14:114] <- probMatrixAdapted
highQModel[116] <- c("100 5")
highQModel[129:132] <- c("0 0 0 0 0") 
highQModel[720] <- c("0 0 0 0 0") 
writeLines(highQModel,con="SRR1501444.highQ_adapted.model")

