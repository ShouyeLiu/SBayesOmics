#' SBayesOmics - simGWASMainOverlap
#'
#' Copyright (C) 2023–2025 Shouye Liu
#'
#' This file is part of the SBayesOmics R package.
#'
#' SBayesOmics is free software: you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation, either version 3 of the License, or
#' (at your option) any later version.
#'
#' SBayesOmics is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU General Public License for more details.
#'
#' You should have received a copy of the GNU General Public License
#' along with this package. If not, see <https://www.gnu.org/licenses/>.



simGWASMainOverlap <- function(isIndLevelBool = TRUE,
                               seed = 1234,
                               trainBfile,
                               mapSnpGeneAcrossChrFile,
                               mapLd2SnpFile,
                               simModel = "a", ## a: caual
                               geneOverlap = "a",
                               traitType = "a",
                               indNum = 5000,
                               NumCGCau = 5,      ## num of causal genes in causal model
                               NumCVPerGCau = 10,
                               NumCGPle = 5,      ## for causality model
                               NumCVPerGPle = 10,  ##
                               NumCGNull = 5,
                               NumCVPerGNull = 10,
                               NumCVIG = 10, ## for sbayesC
                               realGeno =TRUE,
                               ldwBool = FALSE,
                               calcLDBlool = FALSE,
                               cisOnly = TRUE,
                               h2cis = 0.5,
                               h2snp = 0.5,
                               h2med = 0.3,
                               outPath = "",
                               smrPath = "",
                               isVaryCaus = FALSE,
                               totalOverlapBool = FALSE,
                               smrIndGenePath = smrPath,
                               smrIndFileSuffix = "bo",
                               cauPleWithinGene = FALSE,
                               debugBool = FALSE
){
  ########################################################
  ########################################################
  ### Step 1. Parameter setting
  ########################################################
  ########################################################
  seed = seed
  if(seed==0){print("Note: seed is not set in dateGenerate function")
  } else {
    print(paste0("Note: seed is set in dataGenerate function as ",seed))
    set.seed(seed)
  }

  # if(h2med == 0){
  #   print("Causal model was not used");
  #   simModel <- sub('a','',simModel)
  # }
  # if(h2med = h2snp){
  #   print("Only causal model was used")
  #   simModel <- sub('bc','',simModel)
  # }
  ###############################################
  simModel = simModel
  modelType <- c()
  overType <- c("a","b","c")
  NumCG <- c()
  NumCV <- c()
  if(grepl("a",simModel)){
    modelType <- c(modelType,"cau")
    NumCG[["cau"]] <- NumCGCau
    if(isVaryCaus){
      NumCV[["cau"]] <- rpois(NumCGCau, NumCVPerGCau) + 1
    }else {
      NumCV[["cau"]] <- rep(NumCVPerGCau,NumCGCau)
    }

  }
  if(grepl("b",simModel)){
    modelType <- c(modelType,"ple")
    NumCG[["ple"]] <- NumCGPle
    if(isVaryCaus){
      NumCV[["ple"]] <- rpois(NumCGPle, NumCVPerGPle) + 1
    }else {
      NumCV[["ple"]] <- rep(NumCVPerGPle,NumCGPle)
    }
  }

  if(grepl("c",simModel)){
    modelType <- c(modelType,"noteQTL")
    NumCV[["noteQTL"]] <- NumCVIG
    if(cisOnly){
      print("cisOnly doesn't work when snps out of gene are used")
      print("set cisOnly = FALSE instead")
      cisOnly = FALSE
    }
  }

  if(grepl("d",simModel)){
    modelType <- c(modelType,"nullGene")
    # minNumCVInNullGene <- NumCVPerGNull + 1
    ## number of genes and the number causal variats per gene depends on each simulation
    NumCG[["nullGene"]] = NumCGNull
    if(isVaryCaus){
      NumCV[["nullGene"]] <- rpois(NumCGNull, NumCVPerGNull) + 1
    }else {
      NumCV[["nullGene"]]  <- rep(NumCVPerGNull,NumCGNull)
    }
  }

  #####################################
  if(!grepl("a",simModel)){
    NumCG[["cau"]] = 0
    NumCV[["cau"]] = 0
    h2med = 0
    NumCGCau = 0      ## for causality model
    NumCVPerGCau = 0  ## for causality model
  }


  if(!grepl("b",simModel)){
    NumCG[["ple"]] = 0
    NumCV[["ple"]] = 0
    NumCGPle = 0      ## for causality model
    NumCVPerGPle = 0  ##
  }
  if(!grepl("c",simModel)){
    NumCG[["noteQTL"]] <- 0
    NumCV[["noteQTL"]] <- 0
    NumCVIG  = 0
    print("snps in intergenic region are not selected.")
    print("set cisOnly = TRUE instead")
    cisOnly = TRUE
  }

  ## for null gene
  if(!grepl("d",simModel)){
    NumCG[["nullGene"]] <- 0
    NumCV[["nullGene"]] <- 0
    NumCVPerGNull  = 0
    NumCGNull     = 0
  }


  if(!grepl("a",simModel) & !grepl("b",simModel) & grepl("c",simModel) ){
    print(paste0("Only SNPs in cis-region and causal model are used, so set ",
                 "h2cis = ",0))
    h2cis = 0
  }
  if(grepl("a",simModel) & !grepl("b",simModel) & !grepl("c",simModel)){
    print(paste0("Only SNPs in cis-region and causal model are used, so set ",
                 "h2med = h2snp = ",h2snp))
    h2med = h2snp
  }
  realGeno =  realGeno
  ldwBool = ldwBool
  indNum = indNum
  h2cis = h2cis
  h2snp = h2snp
  h2med = h2med
  traitType = traitType
  geneOverlap = geneOverlap
  cisOnly = cisOnly ## remove snps out of gene region
  mapSnpGeneAcrossChr <- readRDS(mapSnpGeneAcrossChrFile)
  mapLd2Snp <- readRDS(file = mapLd2SnpFile)  ## only used in step 6
  smrPathBool = FALSE
  if(nchar(smrPath) !=0){smrPathBool = TRUE}
  ########################################################
  ########################################################
  ### Step 2. select candidate gene and snp sets
  ## In this step, we need to deal with gene overlap situation
  ## and different model problems
  ########################################################
  ##c: intergenic model:  snp =/= gene; gene =/= trait; snp ==> trait;
  ##d: nullGene  model:   snp ==> gene; gene =/= trait; snp =/= trait;
  ##b: pleiotropic model: snp ==> gene; gene =/= trait; snp ==> trait;
  ##a: causal model:      snp ==> gene; gene ==> trait; snp =/= trait;
  ########################################################
  genelist <- names(mapSnpGeneAcrossChr$mapGene2Snp)
  totalGeneNum <- length(genelist)
  mapGene2Snp <- mapSnpGeneAcrossChr$mapGene2Snp
  mapSnp2Gene <- mapSnpGeneAcrossChr$mapSnp2Gene
  ########################################################
  ### Step 2.1 causal model and (or) pleiotropic model.
  ########################################################

  geneRemoved <- c()
  snpRemoved  <- c()
  snplistAllInSimGene <- c()
  genelistForModel <- c()
  snplistSetForModel <- c()
  snplistMapForModel <- c()

  for(mod in modelType){
    if(mod == "noteQTL" ) {next;}
    print("###########################################")
    print(paste0("[",mod,"] model is used for simulation"))
    print("###########################################")
    genelistSlted  <- c()
    snplistSetSlted <- c()
    snplistMapslted <- c()

    NumCauGenes <- NumCG[[mod]]
    NumCauVariants <- NumCV[[mod]]
    ge <- NumCauGenes
    ###########################################################
    ### a: not overlap situation.
    ## we need to select gene first and then find all possible
    ## genes overlapped with this gene;
    ## finally, gene with at the start will be chosen and
    ## other related genes are muted.
    ###########################################################
    if(geneOverlap == overType[1]){
      if(isVaryCaus){stop("non-overlap situation not works for vary caus")}
      while(ge){
        ## gene <- sample(genelist[!genelist %in% unique(c(geneRemoved, genelistSlted)) ],1)
        gene <- sample(setdiff(genelist, unique(c(geneRemoved, genelistSlted)) ),1)
        snpInGene <- names(mapGene2Snp[[gene]])
        geneOverlapped <- c()
        for(singleSnp in snpInGene ){
          geneOverlapped <- unique(c(geneOverlapped,mapSnp2Gene[[singleSnp]]))
        }
        ## not overlap
        geneCandidate <- geneOverlapped[1]
        ## geneRemoved <- c(geneRemoved,geneOverlapped[!geneOverlapped %in% geneCandidate])
        geneRemoved <- unique(c(geneRemoved, setdiff( geneOverlapped,geneCandidate) ) )
        ## remove gene from end
        if(length(mapGene2Snp[[geneCandidate]]) < NumCauVariants[ge] )
        {
          geneRemoved <- c(geneRemoved,geneCandidate)
        }else {
          ## store gene info
          selectedCount <- NumCauVariants[ge]
          ge <- ge - 1
          genelistSlted <- c(genelistSlted,geneCandidate)
          snps <- names(sample(mapGene2Snp[[geneCandidate]],selectedCount))
          snplistSetSlted <- unique( c(snplistSetSlted,  snps) )
          snplistMapslted[[geneCandidate]] <- snps
          ## store all snps in gene
          snplistAllInSimGene <- c(snplistAllInSimGene, names(mapGene2Snp[[geneCandidate]]) ) ## this store all snps in genes
          print(paste0(ge,"/",totalGeneNum ,"-th GENE [",geneCandidate,
                       "] was assigned [", NumCauVariants[ge] ,"] eqtls effect!"))
        }
      }
    }
    ########################################################
    ### b: overlap situation.
    ########################################################
    if(geneOverlap == overType[2]){
      if(isVaryCaus){stop("overlap situation not works for vary caus")}
      while(ge){
        ## overlap situation
        gene <- sample(setdiff(genelist, unique(c(geneRemoved, genelistSlted)) ),1)
        snpInGene <- mapGene2Snp[[gene]]
        geneOverlapped <- mapSnp2Gene[[names(which(snpInGene == max(snpInGene ))[1])]]

        ## select snps with effect for at least two genes
        tmpGenes <- c()
        for(geneI in geneOverlapped){
          snpsTotal <- names(mapGene2Snp[[geneI]][mapGene2Snp[[geneI]] > 1] ) ## this makes sure that genes are overlapped.
          if(length(snpsTotal) < NumCauVariants[ge] ){
            geneRemoved <- c(geneRemoved,geneI)
          } else {
            tmpGenes <- c(tmpGenes,geneI)
          }
        }
        ## select overlap genes
        if(length(tmpGenes) >= ge)
        {
          ## if current gene numbers larger than or equal to given threshold
          ## geneCandidate <- sample(tmpGenes,ge)
          geneCandidate <- tmpGenes[1:ge]
          geneRemoved <- c(geneRemoved,tmpGenes)
          ## take intersection only for debugging under two genes
          if(totalOverlapBool){
            print("debugging mode for two genes with totally overlap has been used")
            snplistInterSect <- names(mapGene2Snp[[geneCandidate[1]]])
            for(geneI in geneCandidate[c(-1)]   ){
              snplistInterSect <- intersect(snplistInterSect, names(mapGene2Snp[[geneI]] ) )
            }
            for(geneI in geneCandidate){
              snpsTotal <- names(mapGene2Snp[[geneI]])
              snpsTotal <- intersect(snpsTotal,snplistInterSect)
              snps <- (sample(snpsTotal,NumCauVariants[ge]))
              snplistSetSlted <- unique( c(snplistSetSlted,  snps) )
              snplistMapslted[[geneI]] <- snps
              snplistAllInSimGene <- unique(c(snplistAllInSimGene, snplistInterSect))
              print(paste0(ge,"/",totalGeneNum ,"-th GENE [",geneI,
                           "] was assigned [", NumCauVariants[ge] ,"] eqtls effect!"))
            }
          } else {
            for(geneI in geneCandidate){
              snpsTotal <- names(mapGene2Snp[[geneI]][mapGene2Snp[[geneI]] > 1] )
              snps <- (sample(snpsTotal,NumCauVariants[ge]))
              snplistSetSlted <- unique( c(snplistSetSlted,  snps) )
              snplistMapslted[[geneI]] <- snps
              ## store all snps in gene
              snplistAllInSimGene <- c(snplistAllInSimGene, names(mapGene2Snp[[geneI]]) )
              print(paste0(ge,"/",totalGeneNum ,"-th GENE [",geneI,
                           "] was assigned [", NumCauVariants[ge] ,"] eqtls effect!"))
            }
          }

          ## update ge
          ge = 0
        }else if (length(tmpGenes) < ge )
        {
          ## if current gene numbers less than given threshold,
          ## loop will be used
          geneCandidate <- tmpGenes
          geneRemoved <- c(geneRemoved,tmpGenes)
          for(geneI in geneCandidate){
            snpsTotal <- names(mapGene2Snp[[geneI]][mapGene2Snp[[geneI]] > 1] )
            snps <- (sample(snpsTotal,NumCauVariants[ge]))
            snplistSetSlted <- unique( c(snplistSetSlted,  snps) )
            snplistMapslted[[geneI]] <- snps
            ## store all snps in gene
            snplistAllInSimGene <- c(snplistAllInSimGene, names(mapGene2Snp[[geneI]]) )
            print(paste0(ge,"/",totalGeneNum ,"-th GENE [",geneI,
                         "] was assigned [", NumCauVariants[ge] ,"] eqtls effect!"))
          }
          ge <- ge - length(tmpGenes)
          if(ge < 0) {ge = 0}
        }
        genelistSlted <- c(genelistSlted,geneCandidate)
      }
    }
    ########################################################
    ### c: select genes and snps regardless of overlap
    ########################################################
    if(geneOverlap == overType[3]){
      while (ge){
        gene <- sample(setdiff(genelist, unique(c(geneRemoved, genelistSlted)) ),1)
        ## gene <- sample(setdiff(genelist, unique(c( genelistSlted)) ),1)
        snpNumInGene <- length(mapGene2Snp[[gene]])
        if(snpNumInGene > NumCauVariants[ge]){
          ## select causal genes and cv
          genelistSlted <- c(genelistSlted,gene)
          geneRemoved <- c(geneRemoved,gene)
          snps <- names(sample(mapGene2Snp[[gene]],NumCauVariants[ge]))
          snplistSetSlted <- unique( c(snplistSetSlted,  snps) )
          snplistMapslted[[gene]] <- snps
          ## store all snps in gene
          snplistAllInSimGene <- c(snplistAllInSimGene, names(mapGene2Snp[[gene]]) )
          print(paste0(ge,"/",totalGeneNum ,"-th GENE [",gene,
                       "] was assigned [", NumCauVariants[ge] ,"] eqtls effect!"))
          ge <- ge - 1
        }
      }
    }

    genelistForModel[[mod]] <- genelistSlted
    snplistSetForModel[[mod]] <- snplistSetSlted
    snplistMapForModel[[mod]] <- snplistMapslted
    ## genelist <- genelist[!genelist %in% c(genelistForModel[[mod]])]
    names(NumCV[[mod]]) <- rev(genelistSlted)
    genelist <- setdiff(genelist, genelistForModel[[mod]])
  }
  ########################################################
  ### Step 2.3 intgenic model
  ########################################################
  mapSnpOutGene <- mapSnpGeneAcrossChr$mapSnpOutGene[order(mapSnpGeneAcrossChr$mapSnpOutGene$POS),]
  snpFrommapSnpOutGene <- mapSnpOutGene[["SNP"]]
  totalMkNumOutGene <- 0
  if(grepl("c",simModel) ){
    mod = "noteQTL"
    print("###########################################")
    print(paste0(simModel, ": bayesC model is used for simulation") )
    print("###########################################")
    totalMkNumOutGene <-  length(snpFrommapSnpOutGene)
    snplistSetForModel[[mod]] <- c()
    sn <- NumCV[[mod]]
    while(sn){
      snps <- setdiff(snpFrommapSnpOutGene,snplistSetForModel[[mod]])
      snp <- sample(snps,1)
      snplistSetForModel[[mod]] <- unique(c(snplistSetForModel[[mod]],snp) )
      sn <- sn -1
      print(paste0("SNP [",snp," ] was selected at ",
                   sn,"/",length(snpFrommapSnpOutGene),"-th !!"))
    }
  }

  ########################################################
  ########################################################
  ### Step 3. Assign effect size to gene and snps
  ########################################################
  ########################################################
  ########################################################
  ### Step 3.1 Initialize alphaTrue and betaTrue
  ########################################################
  mapGeneOrdered <- data.frame(GENE = names(mapSnpGeneAcrossChr$mapGene2Snp))


  if(!grepl("a",simModel) & !grepl("b",simModel) & !grepl("d",simModel) & grepl("c",simModel) ){
    snpIdxPos <- unique( mapSnpOutGene)
  }else {
    mapSnpInGene <- mapSnpGeneAcrossChr$mapSnpInGene[order(mapSnpGeneAcrossChr$mapSnpInGene$POS),]
    snpIdxPos <- unique(rbind(mapSnpInGene, mapSnpOutGene))
  }
  snpIdxPos  <- snpIdxPos[order(snpIdxPos$POS),]
  snpIdx  <- snpIdxPos[["SNP"]]


  if(grepl("a",simModel) | grepl("b",simModel) | grepl("d",simModel)){
    ## order snplistAllInSimGene: snps from all cis-region
    snplistAllInSimGene <- unique(snplistAllInSimGene)
    snplistAllInSimGene <- merge(mapSnpInGene,data.frame(SNP = snplistAllInSimGene),
                                 by = "SNP",sort=FALSE)[["SNP"]]  ## sort = FALSE means will keep order

    if(cisOnly){ snpIdx = snplistAllInSimGene}
    ## snpInGeneIdx: snps of caus
    # snpInGeneIdx <- unique(c(snplistSetForModel[["cau"]],snplistSetForModel[["ple"]]) )
    # snpInGeneIdx <- merge(mapSnpInGene,data.frame(SNP = snpInGeneIdx),by = "SNP",sort=FALSE)[["SNP"]]

    geneIdx  <- unique(c(genelistForModel[["cau"]],genelistForModel[["ple"]],genelistForModel[["nullGene"]] ))
    ## order geneIdx
    geneIdx <- merge(mapGeneOrdered,data.frame(GENE = geneIdx),
                     by = "GENE",sort=FALSE)[["GENE"]]
    ## for causal and pleiotropic gene and snps
    thetaTrue <- data.frame(theta = rep(0,length(geneIdx)),row.names = geneIdx)
    ## note here, we need to use all snpIdx
    alphaTrue <- matrix(0, nrow= length(snplistAllInSimGene),
                        ncol= length(geneIdx))
    rownames(alphaTrue) <- snplistAllInSimGene
    colnames(alphaTrue) <- geneIdx
  }

  betaTrue <- data.frame(beta = rep(0,length(snpIdx)),row.names = snpIdx)
  ##########################################################
  ### Step 3.2 simulate snp effects for snps from intergenic
  ##########################################################
  cauSnpNumOutGene <- 0
  if(grepl("b",simModel)| grepl("c",simModel)){
    ## sampling snp effect to all selected snps
    snpNotMedIdx <- unique( c(snplistSetForModel[["ple"]],snplistSetForModel[["noteQTL"]]) )
    snpNotMedIdx <-  merge(snpIdxPos,data.frame(SNP = snpNotMedIdx),by = "SNP",sort=FALSE) [["SNP"]]
    NumSnpNotMed <- length(snpNotMedIdx)
    sigmaForNotMed <- (h2snp - h2med)/NumSnpNotMed
    betaTrue[match(snpNotMedIdx,rownames(betaTrue)),] <- rnorm(NumSnpNotMed,0,sqrt( sigmaForNotMed) )
    cauSnpNumOutGene <- length(snplistSetForModel[["noteQTL"]])
  }
  ########################################################
  ### Step 3.3 use h2cis and h2med simulate gene effects
  ### and alpha effects
  ########################################################
  if(grepl("a",simModel) | grepl("b",simModel) | grepl("d",simModel) ){
    ## sigmaForCis <- 1
    ## sampling eqtl effect in genes
    if(NumCG[["cau"]] > 0){
      ## sampling gene effect for causal model,
      ## but for pleiotropic model, gene effect is zero
      sigmaForGene <- h2med/(NumCG[["cau"]] * h2cis)
      thetaTrue[match(genelistForModel[["cau"]],rownames(thetaTrue)),] <- rnorm(NumCG[["cau"]],0,sqrt( sigmaForGene) )
      ## thetaTrue[match(genelistForModel[["cau"]],rownames(thetaTrue)),]  = 1
      ##
      for(i in 1: NumCG[["cau"]])
      {
        gene = genelistForModel[["cau"]][i]
        if(NumCV[["cau"]][[gene]] == 0) {next;}
        sigmaForCis <- h2cis/ (NumCV[["cau"]][[gene]] )
        alphaTrue[match(snplistMapForModel[["cau"]][[gene]],rownames(alphaTrue)),gene] <- rnorm(NumCV[["cau"]][[gene]],0,sqrt(sigmaForCis) )
      }
    }

    if(NumCG[["ple"]] > 0){
      for(i in 1: NumCG[["ple"]])
      {
        gene = genelistForModel[["ple"]][i]
        if(NumCV[["ple"]][[gene]] == 0) {next;}
        sigmaForCis <- h2cis/ ( NumCV[["ple"]][[gene]] )

        ## Step 1. sampling alphaTrue first
        alphaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(alphaTrue)),gene] <-  rnorm(NumCV[["ple"]][[gene]],0,sqrt(sigmaForCis) )
        ## Step 2. calculate correlation between alpha and beta from pleiotropic snps
        ## resampling alpha if the correlation between beta and alpha > 0.1 or < -0.1.
        corBetaAlpah <- cor(betaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(betaTrue)),],
                            alphaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(alphaTrue)),gene] )
        # while(abs(corBetaAlpah) > 1e-3){
          while((corBetaAlpah)^2 > 1e-5){
          alphaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(alphaTrue)),gene] <-  rnorm(NumCV[["ple"]][[gene]],0,sqrt(sigmaForCis) )
          corBetaAlpah <- cor(betaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(betaTrue)),],
                              alphaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(alphaTrue)),gene] )
        }
        # print(paste0("cor: ", corBetaAlpah))
      } ## end of gene
    } ## end of numCG

    if(NumCG[["nullGene"]] > 0){
      for(i in 1: NumCG[["nullGene"]])
      {
        gene = genelistForModel[["nullGene"]][i]
        if(NumCV[["nullGene"]][[gene]] == 0) {next;}
        sigmaForCis <- h2cis/ ( NumCV[["nullGene"]][[gene]] )
        alphaTrue[match(snplistMapForModel[["nullGene"]][[gene]],rownames(alphaTrue)),gene] <-  rnorm(NumCV[["nullGene"]][[gene]],0,sqrt(sigmaForCis) )
      }
    }
  }

  #########################################################
  ### Step 3.5 construct A matrix based on gene-snp set ###
  #########################################################
  snpIdx <- snpIdx
  cauGeneNum <- 0
  totalGeneNum <- 0
  cauSnpNumInGene <- 0
  if(grepl("a",simModel)| grepl("b",simModel) | grepl("d",simModel) ){
    geneIdx <- colnames(alphaTrue)
    cauGeneNum <- cauGeneNum + length(genelistForModel[["cau"]])
    totalGeneNum <- totalGeneNum + length(geneIdx)
    cauSnpNumInGene <- length(snplistSetForModel[["cau"]]) + length(snplistSetForModel[["ple"]]) + length(snplistSetForModel[["nullGene"]])
    ## cauSnpNumInGene <- length(unique( c(snplistSetForModel[["cau"]], snplistSetForModel[["ple"]], snplistSetForModel[["nullGene"]]) ) )
  }
  cauSnpNum <- length(unique(as.vector(unlist(snplistSetForModel)))  )

  ############################################################
  ############################################################
  ###   Step 4 construct cleaned genotype matrix           ###
  ### This step is mainly to scale and impute raw genotype ###
  ############################################################
  ############################################################
  if(realGeno) {
    if(FALSE){
      bedData <-BEDMatrix::BEDMatrix(trainBfile)
      geno <- as.matrix(bedData)
      W <- BGData::preprocess(geno, center = FALSE, scale = FALSE,impute = TRUE)
      ## replace colname of W "rs6600755_A" as "rs6600755"
      colnames(W)  <- gsub("(.+)_(.+)","\\1",colnames(W))
      if(indNum == 0){
        indNum = dim(geno)[1]
        selectedIndName = indNum
      } else{
        selectedIndName <- sample(rownames(W),indNum)
      }
      W <- W[match(selectedIndName,rownames(W)),]
      W <- apply(W, 2, function(x){ (x-mean(x))/sd(x) })
      ## W <- apply(W, 2, function(x){ x-mean(x) })
      print(paste0("REAL genotype is used ", "(N = ",indNum,")"))
    }
    if(TRUE){
      print(paste0("Train file: ", trainBfile))
      bedData <- BGData::as.BGData(BEDMatrix::BEDMatrix(trainBfile))
      bedIds <- gsub("(.+)_(.+)","\\1",colnames(bedData@geno))
      selectedColumns <- match(snpIdx, bedIds)
      if (anyNA(selectedColumns)) stop("Requested SNP missing from BED")
      geno <- as.matrix(bedData@geno[,selectedColumns,drop=FALSE])
      colnames(geno) <- snpIdx
      ## select phenotype
      # selectedIndName <- sample(rownames(geno),indNum)
      if(indNum == 0){
        selectedIndName = rownames(geno)
      } else{
        selectedIndName <- sample(rownames(geno),indNum)
      }
      geno <- geno[match(selectedIndName,rownames(geno)),]
      ## extract genotype information
      genoMap <- bedData@map[ bedData@map$snp_id %in% snpIdx ,]
      rownames(genoMap) = genoMap$snp_id
      pseudoYForGeno = data.frame(FID = gsub("(.+)_(.+)","\\1", rownames(geno) ),
                                  IID = gsub("(.+)_(.+)","\\2", rownames(geno) ),
                                  PAT = 0,MAT = 0,SEX = 0,PHENOTYPE = -9)
      rownames(pseudoYForGeno) = rownames(geno)
      genoInfo <- BGData::BGData(geno = geno,pheno = as.data.frame(pseudoYForGeno),map = genoMap)
      freq  <- BGData::summarize(X = BGData::geno(genoInfo) )
      ## add additional attributes
      genoInfo@map$freq_na      <- freq$freq_na
      genoInfo@map$sd           <- freq$sd
      genoInfo@map$allele_freq  <- freq$allele_freq
      genoInfo@map$na           <- colSums(is.na(geno))
      snp2pq                    <- 2 * freq$allele_freq * ( 1- freq$allele_freq)
      names(snp2pq)             <- snpIdx
      ## finally, center and standardize the data
      W <- BGData::preprocess(geno, center = TRUE, scale = TRUE,impute = TRUE)
      ## W <- apply(W, 2, function(x){ (x-mean(x))/sd(x) })
      ## W <- apply(W, 2, function(x){ x-mean(x) })
      print(paste0("REAL ddd genotype is used ", "(N = ",length(selectedIndName),")"))
    }

  }else {
    mkNum <- length(snpIdx)
    geno <- matrix(sample(0:2, indNum * mkNum, replace = TRUE), nrow = indNum )
    rownames(geno) <- paste0("ind-",1:indNum)
    selectedIndName = rownames(geno)
    colnames(geno) <- snpIdx
    pseudoX <- apply(geno, 2, function(x){ (x-mean(x))/sd(x) })
    ## center only
    # pseudoX <- apply(geno, 2, function(x){ (x-mean(x)) })
    print(paste0("PSEUDO genotype is used ", "(N = ",indNum,")"))
  }
  #######################################################
  #####   Step 4.1 summary genoInfo and geneInfor   #####
  #######################################################
  if(realGeno){
    X <-  W[match(selectedIndName,rownames(geno)),match(snpIdx,colnames(W))]
    remove(W)
  } else {
    X <- pseudoX
    remove(pseudoX)
    snp2pq = apply(X,2,var)
  }
  ########################################################
  #######################################################
  ########## Step 5. Simulation  Main function  #########
  #######################################################
  ########################################################

  #######################################################
  ### Step 5.1 Parameter setting
  #######################################################
  ## heritability
  h2cis <- h2cis
  h2snp <- h2snp
  ## genotypes
  indNum <- dim(X)[1]; mkNum <-  dim(X)[2]
  ## snp
  betaTrue <- as.matrix(betaTrue)
  samlist <- rownames(X)
  ## gene
  if(grepl("a",simModel)| grepl("b",simModel) | grepl("d",simModel)){
    geneIdx <- colnames(alphaTrue);
    thetaTrue <- as.matrix(thetaTrue)
  }
  #######################################################
  #######################################################
  ### Step 5.2 Generating eqtl expression
  ### Simulate gene expression level for each gene:
  ### g = X*alpha + e_g where X is the cis-eQTL genotype
  ### matrix and e_g ~ N(0, var(X*alpha)*(1-h2_g)/h2_g).
  #######################################################
  if(grepl("a",simModel) | grepl("b",simModel)| grepl("d",simModel)){
    eMatrix <- matrix(0,nrow = indNum ,ncol = length(geneIdx))
    colnames(eMatrix) <- geneIdx
    rownames(eMatrix) <- selectedIndName
    for (i in 1:length(geneIdx)) {
      eqtlG <- X[,snplistAllInSimGene] %*% alphaTrue[,i]
      eqtlVg <- var(eqtlG)
      ## A null molecular phenotype has no cis genetic component.  It is a
      ## valid EIEO/AIAO simulation state; use unit environmental variance
      ## instead of propagating a zero variance into rnorm (NA).
      if (!is.finite(eqtlVg) || eqtlVg <= 0) eqtlVg <- 0
      ## print(eqtlVg)
      ## eqtlVe <- 1 - h2cis
      eqtlVe <- if (eqtlVg == 0) 1 else (1-h2cis)/h2cis * eqtlVg
      eqtlExpression <- eqtlG + rnorm(indNum,0,sqrt(eqtlVe))
      # eqtlExpression <- (eqtlExpression - mean(eqtlExpression))/sd(eqtlExpression) # used for downsample strategy
      eqtlExpression <- eqtlExpression - mean(eqtlExpression)
      eMatrix[,i] <- eqtlExpression
      print(paste0("Gene expression for gene ",i, " generated"))
    }
    # eMatrixOri = eMatrix
    varGene <- apply(eMatrix, 2, var)
    # eMatrix = scale(eMatrix)

    print(paste0("Gene expression matrix generated"))
  }

  ########################################################
  ### Step 5.3 Generating trait phenotype
  ### Simulate trait phenotype based on all gene expression:
  ### y = G*theta + e_y where G is the gene expression matrix
  ### and e_y ~ N(0, var(G*theta)*(1-h2_y)/h2_y).
  ########################################################
  gwasG <- 0
  if(grepl("b",simModel) | grepl("c",simModel) ){
    ## plieiotropic model
    ## gwasG <- gwasG +  X %*% betaTrue
  }
  if(grepl("a",simModel)){
    ## causal model
    ## gwasG <- gwasG + eMatrix %*% as.matrix(thetaTrue)
    ## gwasG <- gwasG + X %*% alphaTrue %*% as.matrix(thetaTrue)
    betaTrue[snplistSetForModel[["cau"]],] <- alphaTrue[snplistSetForModel[["cau"]],] %*% thetaTrue
  }
  gwasG <- gwasG +  X %*% betaTrue

  ## generate trait phenotype
  if (traitType == "a") {
    gwasVg <- var(gwasG)
    # here we need to simulate simulation
    ## gwasVe <- 1 - h2snp
    gwasVe <- (1- h2snp)/(h2snp) * (gwasVg)
    y <- gwasG + rnorm(indNum,0,sqrt(gwasVe))
    y <- (y - mean(y))
    # y <- (yOri - mean(yOri))/sd(yOri)
    vary <- var(y)

    print(paste0("Continuous trait vector generated"))
  } else if(traitType == "b") {
    pr = 1/(1+exp(- gwasG))
    y <- rbinom(indNum,1,pr)
    # vary = indNum * sum(y)/indNum *(1- sum(y)/indNum)   # np(1-p)
    print(paste0("binary trait vector generated"))
  } else {

  }
  rownames(y) <- selectedIndName
  # rownames(yOri) <- selectedIndName
  #########################################################
  ### Step 5.4 Generate GWAS summary statistics
  #########################################################
  nGWAS = indNum
  print(paste0("calculate bhat"))
  if (traitType == "a") {
    ## since y is not scaled, here we need to scale it to R
    bhat <- apply(X, 2, function(x){lm(y ~ x)$coef[2]})
    bhatSE <- apply(X, 2, function(x){coef(summary(lm(y ~ x)))[2,2]})
    bhatSMR = bhat /sqrt(snp2pq)
    bhatSESMR = bhatSE /sqrt(snp2pq)
    ## now we scale it
    bhatSqrtScaleFactor = sqrt(1/(indNum * bhatSE*bhatSE + bhat*bhat))
    # bhat <- bhat * bhatSqrtScaleFactor
    # bhatSE = bhatSE * bhatSqrtScaleFactor
    # bhat <- bhat * sqrt(snp2pq)
    # bhatSE = bhatSE * sqrt(snp2pq)

    bhatPvalue <- apply(X, 2, function(x){coef(summary(lm(y ~ x)))[2,4]})
    for(i in 1:length(bhatPvalue)) {if(bhatPvalue[i] <= 1e-8) bhatPvalue[i] = 0}

    print(paste0("calculate bhat from continuous trait"))
  } else if(traitType == "b"){
    bhat <- apply(X, 2, function(x){glm(y ~ x,family = "binomial")$coef[2]})
    bhatSE <- apply(X, 2, function(x){coef(summary(glm(y ~ x,family = "binomial")))[2,2]})
  } else {

  }
  ## calculate vary
  # vary = median(snp2pq * indNum * (indNum * bhatSESMR^2 + bhatSMR^2)/indNum)
  vary = var(y)
  #########################################################
  ### Step 5.5 Generate eqtl summary statistics,
  ### including marginal and joint eqtl effect
  ### only estimate the effects of cis-eQTL for each gene
  #########################################################
  totalMkNumInGeneOverlap <- 0
  neQTLs <- c()
  varGene <- c()
  if(grepl("a",simModel) | grepl("b",simModel) | grepl("d",simModel)){
    AMargin <- matrix(0,nrow=length(snplistAllInSimGene), ncol=length(geneIdx), dimnames=list(snplistAllInSimGene, geneIdx))
    AMarginSE <- matrix(0,nrow=length(snplistAllInSimGene), ncol=length(geneIdx), dimnames=list(snplistAllInSimGene, geneIdx))
    AMarginSqrt2pq <- matrix(0,nrow=length(snplistAllInSimGene), ncol=length(geneIdx), dimnames=list(snplistAllInSimGene, geneIdx))
    AMarginSqrtScaleFactor <- matrix(0,nrow=length(snplistAllInSimGene), ncol=length(geneIdx), dimnames=list(snplistAllInSimGene, geneIdx))
    AMarginPvalue <- matrix(0,nrow=length(snplistAllInSimGene), ncol=length(geneIdx), dimnames=list(snplistAllInSimGene, geneIdx))
    AMarginSMR <- matrix(0,nrow=length(snplistAllInSimGene), ncol=length(geneIdx), dimnames=list(snplistAllInSimGene, geneIdx))
    AMarginSESMR <- matrix(0,nrow=length(snplistAllInSimGene), ncol=length(geneIdx), dimnames=list(snplistAllInSimGene, geneIdx))
    for (i in 1:length(geneIdx))
    {
      if(i == 1){print(paste0("eQTL margianl and joint effect is generated from dataset same as GWAS!")) }
      ## cis-SNP for the gene
      snpID = intersect(names(mapGene2Snp[[geneIdx[i]]]), colnames(X) )
      ######################################################################
      XInGene <- X[,snpID]
      ## beta
      bInGene <-  apply(XInGene, 2, function(x){lm(eMatrix[,i] ~ x)$coef[2]})
      names(bInGene) <- snpID
      ## se
      bInGeneSE <- apply(XInGene, 2, function(x){coef(summary(lm(eMatrix[,i] ~ x)))[2,2]})
      names(bInGeneSE) <- snpID

      AMargin[snpID,geneIdx[i]] <- bInGene
      AMarginSE[snpID,geneIdx[i]] <- bInGeneSE
      AMarginSMR[snpID,geneIdx[i]] <- bInGene /sqrt(snp2pq[snpID])
      AMarginSESMR[snpID,geneIdx[i]] <- bInGeneSE /sqrt(snp2pq[snpID])
      ## we scale it now
      AMarginSqrtScaleFactor[snpID,geneIdx[i]] <- sqrt(1/(indNum * bInGeneSE*bInGeneSE + bInGene*bInGene))
      # AMargin[snpID,geneIdx[i]] <- bInGene * sqrt(snp2pq[snpID])
      # AMarginSE[snpID,geneIdx[i]] <- bInGeneSE *  sqrt(snp2pq[snpID])
      # AMarginSqrt2pq[snpID,geneIdx[i]]  = sqrt(snp2pq[snpID])

      # AMargin[snpID,geneIdx[i]] <- bInGene * AMarginSqrtScaleFactor[snpID,geneIdx[i]]
      # AMarginSE[snpID,geneIdx[i]] <- bInGeneSE *  AMarginSqrtScaleFactor[snpID,geneIdx[i]]
      AMarginSqrt2pq[snpID,geneIdx[i]]  = sqrt(snp2pq[snpID])
      ## pvalue
      bInGenePvalue <- apply(XInGene, 2, function(x){coef(summary(lm(eMatrix[,i] ~ x)))[2,4]})
      # for(i in 1:length(bInGenePvalue)) {if(bInGenePvalue[i] <= 1e-8) bInGenePvalue[i] = 0}
      names(bInGenePvalue) <- snpID
      AMarginPvalue[snpID,geneIdx[i]] <- bInGenePvalue
      ###### calculate varGene
      svarGene = median(snp2pq[snpID] * indNum * (indNum * AMarginSESMR[snpID,geneIdx[i]]^2 + AMarginSMR[snpID,geneIdx[i]]^2)/indNum)
      svarGene = var(eMatrix[,i])
      names(svarGene) = geneIdx[i]
      varGene = c(varGene,svarGene)
    }
    totalMkNumInGeneOverlap <- sum(colSums(AMargin != 0) )
    neQTLs = apply(eMatrix, 2, length)
  }

  #########################################################
  #########################################################
  #######    Step 6 Save LD matrix information    #########
  #########################################################
  #########################################################

  ###########################################################
  ## Step 6.1 Construct maps among gene, snp and LD blocks ##
  ###########################################################
  geneSnpPosInfo <- mapSnpGeneAcrossChr$geneSnpPosInfo
  if(grepl("a",simModel) | grepl("b",simModel)| grepl("d",simModel) ){
    geneSnpPosInfo <- unique(geneSnpPosInfo[geneSnpPosInfo$GENE %in% geneIdx,])
    mapGene2SnpAcrossBlocks <- lapply(geneIdx, function(gene) {
      ids <- names(mapGene2Snp[[gene]])
      ids <- ids[ids %in% snplistAllInSimGene]
      setNames(vapply(ids, function(snp) sum(geneIdx %in% mapSnp2Gene[[snp]]), integer(1)), ids)
    })
    names(mapGene2SnpAcrossBlocks) <- geneIdx
    mapSnp2GeneAcrossBlocks <- setNames(lapply(snplistAllInSimGene, function(snp) {
      geneIdx[vapply(mapGene2SnpAcrossBlocks, function(x) snp %in% names(x), logical(1))]
    }), snplistAllInSimGene)
    print("Gene overlap situation:")
    print( table(unlist(mapGene2SnpAcrossBlocks)))
  }

  #########################################################
  ########  Step 6.2 Generate per gene LD matrix    #######
  #########################################################
  # AhatBlocks <- c()
  Rgene <- c()
  if(grepl("a",simModel) | grepl("b",simModel)| grepl("d",simModel)){
    genes <- unique(geneSnpPosInfo$GENE)
    if(length(genes) != 0){
      oneGeneR <-c()
      for(oneGene in genes){
        ## Do svd for each gene
        snpInGene <- unique(geneSnpPosInfo[which( geneSnpPosInfo$GENE == oneGene ),][,c("SNP","SNPPOS")] )
        snpInGene <- snpInGene[order(snpInGene$SNPPOS),][["SNP"]]
        XBlock <- X[,match(snpInGene,colnames(X))]
        oneGeneR[[oneGene]] <- cor(X[,snpInGene])
      }
      Rgene <-  oneGeneR

    }
  }
  #########################################################
  ######  Step 6.3 divided chromosome by LD blocks  #######
  #########################################################
  RBlocks <- c()
  DBlocks <- c()
  ldblocks = c()
  if(calcLDBlool){
  if(ldwBool){
    print("LD block is used.")
    ldblocks <- unique(mapLd2Snp$LDBLOCK)
    for (lbs in ldblocks){
      if(FALSE){
        snpInLdBlock <- mapLd2Snp[which(mapLd2Snp$LDBLOCK == lbs),]
        snpInGene <- Reduce(intersect,list(snplistAllInSimGene,snpInLdBlock))
        XBlock <- X[,match(snpInGene,colnames(X))]
        RBlocks[[lbs]] <- cor(XBlock) # WGCNA::cor(XBlock)
        DBlocks[[lbs]] <- apply(XBlock, 2, var) * indNum
      }else {
        ld  <- mapLd2Snp[which(mapLd2Snp$LDBLOCK == lbs),]
        snpInLd <- unique( ld[,c("SNP","SNPPOS")])
        snpInLd <- snpInLd[order(snpInLd$SNPPOS),][["SNP"]]
        snpCommon <- Reduce(intersect, list(colnames(X),snpInLd))
        if(length(snpCommon) == 0) next;
        XBlock <- X[,match(snpCommon,colnames(X))]
        RBlocks[[lbs]] <- cor(XBlock) # WGCNA::cor(XBlock)
        DBlocks[[lbs]] <- apply(XBlock, 2, var) * indNum
      }
    } ## end of ld blocks loop
  }else{
    print("LD block is not used.")
    ldblocks = "one"
    lbs = ldblocks
    RBlocks[[lbs]] <- cor(X) # WGCNA::cor(X)
    DBlocks[[lbs]] <- apply(X, 2, var) * indNum
  }
  }

  ####################################################################
  nGene <- length(Rgene)
  nBlock <- length(RBlocks)
  ## svd for per gene
  Ugene <- c()
  lambdaGene <- c()
  if(nGene != 0){
    oneU <- c()
    oneLambda <- c()
    for (oneGene in names(Rgene) ){
      oneR <- Rgene[[oneGene]]
      eig = eigen(oneR, symmetric=TRUE)
      lambda = eig$values[eig$values > 0]
      selected = which(cumsum(lambda)/sum(lambda) < 0.9999)
      selected <- c(selected, max(selected) + 1)
      oneLambda[[oneGene]] = eig$values[selected]
      oneU[[oneGene]] = eig$vectors[,selected]
      rownames(oneU[[oneGene]]) <- colnames(oneR)
    }
    lambdaGene = oneLambda
    Ugene      = oneU
  }

  ## svd for chromo ld blocks
  UBlocks <- c()
  lambdaBlocks <- c()
  if(calcLDBlool){
  for(lbs in 1:nBlock){
    if(length( rownames (RBlocks[[lbs]] )) != 0){
      R <- RBlocks[[lbs]]
      eig = eigen(R, symmetric=TRUE)
      lambda = eig$values[eig$values > 0]
      selected = which(cumsum(lambda)/sum(lambda) < 0.9999) # 0.995
      selected <- c(selected, max(selected) + 1)
      lambdaBlocks[[lbs]] = eig$values[selected]
      UBlocks[[lbs]] = eig$vectors[,selected]
      rownames(UBlocks[[lbs]]) <- colnames(R)
    }else {
      lambdaBlocks[[lbs]] <- NULL
      UBlocks[[lbs]] <- NULL
    }
  }
  }
  ####################################################################

  varGeneSMR = c()
  varySMR = c()
  gwasDat <- c()
  if(realGeno){
    genoInfoAcrossGene <- genoInfo@map[match(snpIdx,rownames(genoInfo@map)),]
    gwasDat <- data.table::data.table(SNP = snpIdx,
                                      A1 =  genoInfoAcrossGene$allele_1,
                                      A2 = genoInfoAcrossGene$allele_2,
                                      freq = genoInfoAcrossGene$allele_freq,
                                      b = bhat,
                                      bhatSMR = bhatSMR[snpIdx],
                                      bhatSESMR = bhatSESMR[snpIdx],
                                      se = bhatSE[snpIdx],
                                      p = bhatPvalue[snpIdx],
                                      N = indNum )

    gwasDat[, snp2pq:=2 * freq * (1 - freq) ]
    gwasDat[, D:= N* snp2pq]
    gwasDat[,bScale:= b ]
    gwasDat[,seScale:= se  ]
    gwasDat[, varps:=D*(N * seScale^2 + bScale^2)/N]
    gwasDat[, varpsSMR:=D*(N * bhatSESMR^2 + bhatSMR^2)/N]
    varySMR = median(gwasDat$varps)

    geneSumDat <- c()
    if(grepl("a",simModel) | grepl("b",simModel) | grepl("d",simModel)){
      for (i in 1:length(geneIdx)){
        gene = unique(geneSnpPosInfo[which(geneSnpPosInfo$GENE == geneIdx[i]),][,c("GCHR","GENE","midPos","GSTARTORI","GENDORI")])
        ## cis-SNP for each gene
        snpID = intersect(names(mapGene2Snp[[geneIdx[i]]]), colnames(X) )
        ## save results into smr format for use
        genoInfoInGene <- genoInfo@map[match(snpID,rownames(genoInfo@map)),]
        esdDat <- data.frame(
          Chr  =  genoInfoInGene$chromosome,
          SNP  =  genoInfoInGene$snp_id,
          Bp   =  genoInfoInGene$base_pair_position,
          A1   =  genoInfoInGene$allele_1,
          A2   =  genoInfoInGene$allele_2,
          Freq = genoInfoInGene$allele_freq,
          Beta = AMargin[snpID,geneIdx[i]], ##  * sqrt(2 * genoInfoInGene$allele_freq * (1- genoInfoInGene$allele_freq) ),
          se   = AMarginSE[snpID,geneIdx[i]],
          p    = AMarginPvalue[snpID,geneIdx[i]]
        )
        snp2pqeQTL <- 2.0 * esdDat$Freq *( 1- esdDat$Freq)
        DeQTL <-  indNum
        ypySrt <- DeQTL * (indNum * esdDat$se * esdDat$se + (esdDat$Beta)^2)
        varpSrt <- ypySrt/indNum
        ypySrt <- sort(ypySrt)
        varPheeQTL <- median(varpSrt)
        varGeneSMR <- c(varGeneSMR, varPheeQTL)
        geneSumDat[[geneIdx[i]]] <- esdDat
      }
      ## save gene information
      names(varGeneSMR) = geneIdx

    }
  }


  #########################################################
  #########################################################
  piSNP <- cauSnpNum/ mkNum
  if(totalMkNumOutGene == 0) {
    piNonEqtl = 0
  } else {
    piNonEqtl = cauSnpNumOutGene / totalMkNumOutGene
    ## piNoneEqtl
  }
  if(totalMkNumInGeneOverlap ==0){
    piEqtl = 0
    ## piEqtl
  } else {
    piEqtl = cauSnpNumInGene / length(snplistAllInSimGene)
  }
  #########################################################
  #########################################################
  ####### Step 7 save simulation to rds objects   #########
  #########################################################
  #########################################################
  type = paste0("sd-",seed,
                "-sim-",simModel,
                "-trait-",traitType,
                "-ol-",geneOverlap,
                "-gen-",realGeno,
                "-cis-",cisOnly,
                "-ldw-",ldwBool,
                "-sam-",indNum,
                "-snp-",mkNum,
                "-cgc-",NumCGCau,
                "-cvpc-",NumCVPerGCau,
                "-cgp-",NumCGPle,
                "-cvpp-",NumCVPerGPle,
                "-cgn-",NumCGNull,
                "-cvpn-",NumCVPerGNull,
                "-cvo-",NumCVIG,
                "-vcau-",isVaryCaus,
                "-h2c-", h2cis,
                "-h2s-", h2snp,
                "-h2m-", h2med
                ## "-",Sys.Date()
  )

  simRes <- list(NumCGC  = NumCG[["cau"]],
                 NumCVPC = NumCV[["cau"]],
                 NumCGP  = NumCG[["ple"]],
                 NumCVPP = NumCV[["ple"]],
                 NumCGN  = NumCG[["nullGene"]],
                 NumCVPN = NumCV[["nullGene"]],
                 NumCVO  = NumCV[["noteQTL"]],
                 genelistForModel = genelistForModel,
                 snplistSetForModel = snplistSetForModel,
                 snplistMapForModel = snplistMapForModel,
                 snplistAllInSimGene = snplistAllInSimGene,
                 isVaryCaus = isVaryCaus,
                 overlap = geneOverlap,
                 realGeno =realGeno,
                 cisOnly = cisOnly,
                 ldwBool = ldwBool,
                 ldblocks = ldblocks,
                 geneSnpPosInfo = geneSnpPosInfo,
                 mapLd2Snp = mapLd2Snp,
                 ## X = X,
                 h2cis = h2cis,
                 h2snp = h2snp,
                 h2med = h2med,
                 traitType = traitType,
                 betaTrue = betaTrue,
                 bhat = bhat,
                 ## gwas = gwasDat,
                 snp2pq = snp2pq,
                 bhatSqrtScaleFactor = bhatSqrtScaleFactor,
                 bhatSE = bhatSE,
                 bhatSMR = bhatSMR,
                 bhatSESMR = bhatSESMR,
                 RBlocks = RBlocks, ##
                 UBlocks = UBlocks,
                 lambdaBlocks = lambdaBlocks,
                 piSNP = piSNP,
                 piEqtl = piEqtl,
                 piNonEqtl = piNonEqtl,
                 vary = vary,
                 varySMR = varySMR,
                 # genoInfo = genoInfo,
                 gwasDat = gwasDat,
                 # geneSumDat = geneSumDat,
                 nGWAS = nGWAS,
                 D = DBlocks,
                 simModel= simModel,
                 seed = seed,
                 isRealGeno = realGeno,
                 mkNum = mkNum,
                 gwasIndList = selectedIndName,
                 cauSnpNum = cauSnpNum,
                 cauSnpNumOutGene = cauSnpNumOutGene,
                 totalMkNumOutGene = totalMkNumOutGene,
                 cauSnpNumInGene = cauSnpNumInGene,
                 totalMkNumInGeneOverlap = totalMkNumInGeneOverlap,
                 # smrPath = smrPath,
                 smrPathBool = smrPathBool,
                 smrIndFileSuffix = smrIndFileSuffix,
                 type = type
  )

  if(isIndLevelBool){
    ## simRes <- c(simRes,list(X = X,y = y ) )
    simRes[["X"]] = X
    # if(realGeno){ simRes[["snp2pq"]] = 2 * genoInfoAcrossGene$allele_freq * (1 - genoInfoAcrossGene$allele_freq)}
    simRes[["y"]] = y
    # simRes[["yOri"]] = yOri
  } else {
    simRes[["R"]] = RBlocks
  }

  if(grepl("a",simModel) | grepl("b",simModel) |grepl("d",simModel) ){
    simResGene <- list(thetaTrue = thetaTrue,
                       mapGene2Snp = mapGene2Snp,
                       neQTLs = neQTLs,
                       cauGeneNum = cauGeneNum,
                       geneNum    = totalGeneNum,
                       alphaTrue = alphaTrue,
                       AMargin = AMargin,
                       AMarginSMR = AMarginSMR,
                       AMarginSE = AMarginSE,
                       AMarginSqrt2pq = AMarginSqrt2pq,
                       AMarginSqrtScaleFactor = AMarginSqrtScaleFactor,
                       AMarginSESMR = AMarginSESMR,
                       Rgene   = Rgene,
                       Ugene   = Ugene,
                       lambdaGene = lambdaGene,
                       varGene = varGene,
                       varGeneSMR = varGeneSMR,
                       mapSnp2GeneAcrossBlocks = mapSnp2GeneAcrossBlocks,
                       mapGene2SnpAcrossBlocks = mapGene2SnpAcrossBlocks
    )

    if(!isIndLevelBool) {
      simResGene[["Rgene"]] = Rgene
    } else {
      simRes[["Z"]] = X[,snplistAllInSimGene]
      simRes[["eMatrix"]] = eMatrix
      # simRes[["eMatrixOri"]] = eMatrixOri
    }

  }else {
    simResGene <- list(thetaTrue = NA,
                       neQTL = NA,
                       A = NA,
                       AMargin = NA,
                       varGene = NA
    )
  }


  simRes <- c(simRes, simResGene)

  if(outPath != "") {
    ## save all simulate parameters
    simulationOutFile <- paste0(outPath,"/",type,".rds")
    saveRDS(simRes, file=  simulationOutFile )
  }

  ## save simulation results for gctb
  if(smrPathBool){
    ############################################################
    ## Step . Summary-level dataset
    ############################################################
    if(grepl("1000g",smrIndFileSuffix)){
      famDat <- data.frame(fid = gsub("(.+)_(.+)_(.+)","\\1", selectedIndName),
                           iid = gsub("(.+)_(.+)_(.+)","\\2_\\3", selectedIndName)
      )
    }else {
      famDat <- data.frame(fid = gsub("(.+)_(.+)","\\1", selectedIndName),
                           iid = gsub("(.+)_(.+)","\\2", selectedIndName))
    }
    if(TRUE){
      sumLevelPath = paste0(smrPath,"/","summary/")
      commandline <- paste0("mkdir -p ", sumLevelPath)
      system(commandline)
      commandline <- paste0("mkdir -p ", sumLevelPath,"/gene/")
      system(commandline)
      ############################################################
      ## 1. gwas
      ############################################################
      print(paste0("save gwas summary information into ma format",
                   "in path: ",sumLevelPath) )
      # genoInfoAcrossGene <- genoInfo@map[match(snpIdx,rownames(genoInfo@map)),]
      gwasSummaryDat <- data.frame(SNP = snpIdx,
                                   A1 =  genoInfoAcrossGene$allele_1,
                                   A2 = genoInfoAcrossGene$allele_2,
                                   freq = genoInfoAcrossGene$allele_freq,
                                   b = bhatSMR[snpIdx],
                                   se = bhatSESMR[snpIdx],
                                   p = bhatPvalue[snpIdx],
                                   N = indNum
      )
      ## save gwas information
      data.table::fwrite(gwasSummaryDat,file=paste0(sumLevelPath,"/","gwas-",smrIndFileSuffix,"-",type,".ma",spe=""), sep = "\t")
      data.table::fwrite(data.frame(snpIdx = gwasSummaryDat$SNP),file=paste0(sumLevelPath,"/","gwas-",smrIndFileSuffix,"-",type,".ma.snplist",sep=""),
                         col.names =FALSE,
                         sep = "\t")
      ### beta true file
      betaTrueDat = data.frame(QTL = rownames(betaTrue),
                               RefAllele = genoInfoAcrossGene$allele_1,
                               Frequency = genoInfoAcrossGene$allele_freq,
                               Effect = as.vector( betaTrue) )
      data.table::fwrite(betaTrueDat,file=paste0(sumLevelPath,"/","gwas-",smrIndFileSuffix,"-",type,".betaTrue",spe=""), sep = "\t")
      ################################################
      #### ldsc summary statis ofrmat
      gwasLDSCDat <- data.frame(SNP = snpIdx,
                                   N = indNum,
                                   Z = bhatSMR[snpIdx]/bhatSESMR[snpIdx],
                                   A1 =  genoInfoAcrossGene$allele_1,
                                   A2 = genoInfoAcrossGene$allele_2 )
      data.table::fwrite(gwasLDSCDat,file=paste0(sumLevelPath,"/","gwas-",smrIndFileSuffix,"-",type,".sumstats",spe=""), sep = "\t")
      ############################################################
      ## 2. store eQTL data
      ############################################################
      if(grepl("a",simModel) | grepl("b",simModel)){
        eQTLType <- "int"
        print(paste0("save internal eQTL summary information into flist and esd format"," in path: ",sumLevelPath) )
        flistColNames <- c("Chr","ProbeID","GeneticDistance","ProbeBp","Gene","Orientation","PathOfEsd")
        flistDat <- data.frame(matrix(ncol=length(flistColNames),nrow=0,
                                      dimnames=list(NULL, flistColNames)))
        geneSamSize <- data.frame(ProbeID = NULL, indNum = NULL)
        for (i in 1:length(geneIdx)){
          gene = unique(geneSnpPosInfo[which(geneSnpPosInfo$GENE == geneIdx[i]),][,c("GCHR","GENE","midPos","GSTARTORI","GENDORI")])
          ## for each gene
          PathOfEsd = paste0("gene/",eQTLType,"-",smrIndFileSuffix,"-",gene$GENE,"-",type,".esd")
          geneFlist <- data.frame(Chr =  gene$GCHR,
                                  ProbeID = gene$GENE,
                                  GeneticDistance = 0,
                                  ProbeBp = gene$midPos,
                                  Gene = gene$GENE,
                                  N    = indNum,
                                  Orientation = "NA",
                                  PathOfEsd = PathOfEsd
          )
          flistDat <- rbind(flistDat,geneFlist)
          geneSamSize <- rbind(geneSamSize,data.frame(ProbeID = gene$GENE, indNum = indNum))
          ## cis-SNP for each gene
          snpID = intersect(names(mapGene2Snp[[geneIdx[i]]]), colnames(X) )
          ## save results into smr format for use
          genoInfoInGene <- genoInfo@map[match(snpID,rownames(genoInfo@map)),]
          esdDat <- data.frame(
            Chr  =  genoInfoInGene$chromosome,
            SNP  =  genoInfoInGene$snp_id,
            Bp   =  genoInfoInGene$base_pair_position,
            A1   =  genoInfoInGene$allele_1,
            A2   =  genoInfoInGene$allele_2,
            Freq = genoInfoInGene$allele_freq,
            Beta = AMarginSMR[snpID,geneIdx[i]],
            se   = AMarginSESMR[snpID,geneIdx[i]],
            p    = format(AMarginPvalue[snpID,geneIdx[i]], scientific = TRUE)
          )
          snp2pqeQTL <- 2.0 * esdDat$Freq *( 1- esdDat$Freq)
          DeQTL <- snp2pqeQTL * indNum
          # esdDat$Beta = esdDat$Beta/ sqrt(snp2pqeQTL)
          # esdDat$se = esdDat$se/ sqrt(snp2pqeQTL)
          ypySrt <- DeQTL * (indNum * esdDat$se * esdDat$se + (esdDat$Beta)^2)
          varpSrt <- ypySrt/indNum
          ypySrt <- sort(ypySrt)
          varPheeQTL <- median(varpSrt)
          # print(paste0("gene: ", geneIdx[i]," varPhe: ", varPheeQTL) )
          ## save snp marginal effect information
          data.table::fwrite(esdDat,file= paste0(sumLevelPath,"/",PathOfEsd) , sep = "\t")
        }
        ## save gene information
        data.table::fwrite(flistDat,file=paste0(sumLevelPath,"/",eQTLType,"-",smrIndFileSuffix,"-",type,".flist",spe=""), sep = "\t")
        ## save sample list file (fid, iid)
        data.table::fwrite(geneSamSize,file=paste0(sumLevelPath,"/",eQTLType,"-",smrIndFileSuffix,"-",type,".gsam",spe=""), sep = "\t",col.names = FALSE)

      }
    }
    ############################################################
    ## Step . Individual-level dataset
    ############################################################
    indLevelPath = paste0(smrPath,"/","ind/")
    indGenePath = paste0(smrPath,"/","ind/gene/")
    system(paste0("mkdir -p ", indLevelPath))
    system(paste0("mkdir -p ", indGenePath))
    ############################################################
    ########## GWAS###############
    ############################################################
    # fidIdx  <- gsub("(.+)_(.+)","\\1",selectedIndName)
    yDat <- data.frame(FID = famDat$fid,IID = famDat$iid, trait = y)
    data.table::fwrite(yDat,file=paste0(indLevelPath,"/",smrIndFileSuffix,"-",type,".pheno",spe="" ), sep = "\t",col.names = FALSE)
    #######
    samlistDat = data.frame(FID = famDat$fid, IID = famDat$fid)
    data.table::fwrite(famDat,file=paste0(indLevelPath,"/",smrIndFileSuffix,"-",type,".indlist",spe="" ), sep = "\t",col.names = FALSE)

    data.table::fwrite(data.frame(snpIdx = snpIdx),file=paste0(indLevelPath,"/",smrIndFileSuffix,"-",type,".ma.snplist",sep=""),
                       col.names =FALSE,
                       sep = "\t")

    ## betaTrue
    # genoInfoAcrossGene <- genoInfo@map[match(snpIdx,rownames(genoInfo@map)),]
    betaTrueDat = data.frame(QTL = snpIdx,
                             RefAllele = genoInfoAcrossGene$allele_1,
                             Frequency = genoInfoAcrossGene$allele_freq,
                             Effect = as.vector( betaTrue) )
    data.table::fwrite(betaTrueDat,file=paste0(indLevelPath,"/",smrIndFileSuffix,"-",type,".betaTrue",spe=""), sep = "\t")

    ############################################################
    ########### gene ######################
    ############################################################
    if(grepl("a",simModel) | grepl("b",simModel) | grepl("d",simModel)){
      ## thetaTrue
      thetaTrueDat = data.frame(gene = rownames(thetaTrue),
                                Effect = as.vector( thetaTrue) )
      data.table::fwrite(thetaTrueDat,file=paste0(indLevelPath,"/",smrIndFileSuffix,"-",type,".thetaTrue",spe=""), sep = "\t")

      geneColNames <- c("chr","start","end","ensgid","genePath")
      geneDat <- data.frame(matrix(ncol=length(geneColNames),nrow=0,
                                   dimnames=list(NULL, geneColNames)))
      for (i in 1:length(geneIdx)){
        gene = unique(geneSnpPosInfo[which(geneSnpPosInfo$GENE == geneIdx[i]),][,c("GCHR","GENE","midPos","GSTARTORI","GENDORI")])
        ## for each gene
        genePath = paste0(indGenePath,"/",smrIndFileSuffix,"-",gene$GENE,"-",type,".pheno")
        if(grepl("1000g",smrIndFileSuffix)){
          geneFinalPath = paste0("gene","/",smrIndFileSuffix,"-",gene$GENE,"-",type,".pheno")
        } else {
          geneFinalPath = paste0(smrIndGenePath,"/",smrIndFileSuffix,"-",gene$GENE,"-",type,".pheno")
        }
        geneFlist <- data.frame(chr =  gene$GCHR,
                                start = gene$GSTARTORI,
                                end  = gene$GENDORI,
                                ensgid = gene$GENE,
                                genePath = geneFinalPath
        )
        geneDat <- rbind(geneDat,geneFlist)

        perGeneDat <- data.frame(FID = famDat$fid,IID = famDat$iid, gene = eMatrix[,gene$GENE])
        ## save snp marginal effect information
        data.table::fwrite(perGeneDat,file= genePath , sep = "\t")
      }
      ## save gene information
      data.table::fwrite(geneDat,file=paste0(indLevelPath,"/",smrIndFileSuffix,"-","",type,".plist",spe=""), sep = "\t")
      ## save sample list file (fid, iid)
    }

    ######## Here we need to test GCTB C++ in local mac computer
    if(smrIndFileSuffix == "mac"){
      data.table::fwrite(yDat,file=paste0(indLevelPath,"/",smrIndFileSuffix,".pheno",spe="" ), sep = "\t",col.names = FALSE)
      data.table::fwrite(samlistDat,file=paste0(indLevelPath,"/",smrIndFileSuffix,".indlist",spe="" ), sep = "\t",col.names = FALSE)
      data.table::fwrite(geneDat,file=paste0(indLevelPath,"/",smrIndFileSuffix,".plist",spe=""), sep = "\t")
      data.table::fwrite(data.frame(snpIdx = snpIdx),file=paste0(indLevelPath,"/",smrIndFileSuffix,".ma.snplist",sep=""),
                         col.names =FALSE,
                         sep = "\t")
    }

    ############################################################
    ########### genotype ######################
    ############################################################

  }

  print(paste0("Date for simulation generated here !!"))
  return (simRes)

}
