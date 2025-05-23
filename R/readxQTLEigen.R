readxQTLEigen <- function( eigenLDblockFile, eigenType, cutThresh = 0.995){

  snpInfo = fread(paste0(eigenLDblockFile,".eigen.",eigenType, ".snp.info"))
  ldmInfo = fread(paste0(eigenLDblockFile,".eigen.",eigenType, ".info"))
  if(eigenType == "gene"){
    eigenBlocks = unique(ldmInfo$probeID)
  }else if(eigenType == "ldblock"){
    eigenBlocks = unique(ldmInfo$LDBLOCK)
  }


  UBlocks <- c()
  lambdaBlocks <- c()
  ldm = paste0(eigenLDblockFile,".eigen.",eigenType, ".bin")
  fp = file(ldm, "rb")
  i = 0;
  for(idx in eigenBlocks){
    curBlock = idx;
    cur_m = 0;
    cur_k = 0;
    sumLambda = 0;
    cur_m  = readBin(fp, "integer", size = 4, n = 1)
    cur_k  = readBin(fp, "integer", size = 4, n = 1 )
    sumLambda = readBin(fp, "double", size = 4, n = 1)
    eigenCutoff = readBin(fp, "double", size = 4, n = 1)
    # if(abs(eigenCutoff - 0.9995) > 1e-6 ){
    #   message("Error in block ",i, " in eigenCutoff = ", eigenCutoff)
    #   break;
    # }
    lambda = readBin(fp, what="double", size=4, n=cur_k)
    U = matrix(readBin(fp, what="double", size=4, n=cur_m * cur_k), ncol = cur_k)
    ## cut thresh further
    haveValue = FALSE;
    revIdx = 0;
    if(eigenCutoff != cutThresh & i == 0){
      ## message("Warning: current proportion of variance in LD block is set as " + to_string(cutThresh)+ ". But the proportion of variance is set as "<< to_string(eigenCutoff) + " in "  + svdLDfile, ".\n")
    }

    lambdaOne = lambda
    UOne = U
    if(eigenType == "gene") {
      snplist = ldmInfo[which(ldmInfo$probeID == idx),][["snpInGene"]]
    } else if(eigenType == "ldblock"){
      snplist = ldmInfo[which(ldmInfo$LDBLOCK == idx),][["snpInLdBlock"]]
    } else {

    }
    rownames(UOne) = snplist
    UBlocks[[ as.character(curBlock)]] <- UOne
    lambdaBlocks[[as.character(curBlock )]] <- lambdaOne
    i = i + 1
  }
  close(fp)
  if(eigenType == "gene") {
    message(length(UBlocks)," gene low-rank marices have been read.")
  } else if(eigenType == "ldblock"){
    message(length(UBlocks)," block low-rank marices have been read.")
  } else {

  }

  out <- list(UBlocks = UBlocks,
              lambdaBlocks = lambdaBlocks)
}

