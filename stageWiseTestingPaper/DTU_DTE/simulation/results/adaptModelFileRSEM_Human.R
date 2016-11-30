setwd("/Users/koenvandenberge/PhD_Data/dtu/diff_splice_paper_Kvdb/hsapiens/reference_files/rsem_model/SRR493366.stat/")
highQModel=readLines("SRR493366.model")
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

## additional adaptations
highQModel[14:114] <- probMatrixAdapted
highQModel[116] <- c("100 5")
transMat <- highQModel[c(129:132,720)]
transMatList <- sapply(1:length(transMat),function(i) as.numeric(strsplit(transMat[i],split=" ")[[1]]), simplify=FALSE)
for(i in 1:length(transMatList)) transMatList[[i]][5]=0
transMatListNormalized <- lapply(transMatList, function(x) x/sum(x))
#retransform to original format
transMatListAdapted <- vector()
for(i in 1:length(transMatListNormalized)) transMatListAdapted[i] <- Reduce(paste,as.character(transMatListNormalized[[i]]))
highQModel[c(129:132,720)] <- transMatListAdapted
writeLines(highQModel,con="SRR493366.highQ_adapted.model")

