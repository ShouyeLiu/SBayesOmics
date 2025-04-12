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
                               simModel = "a",
                               geneOverlap = "a",
                               traitType = "a",
                               indNum = 5000,
                               NumCGCau = 5,
                               NumCVPerGCau = 10,
                               NumCGPle = 5,
                               NumCVPerGPle = 10,
                               NumCGNull = 5,
                               NumCVPerGNull = 10,
                               NumCVIG = 10,
                               realGeno =TRUE,
                               ldwBool = FALSE,
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
  seed = seed
  if(seed==0){print("Note: seed is not set in dateGenerate function")
  } else {
    print(paste0("Note: seed is set in dataGenerate function as ",seed))
    set.seed(seed)
  }
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
    NumCG[["nullGene"]] = NumCGNull
    if(isVaryCaus){
      NumCV[["nullGene"]] <- rpois(NumCGNull, NumCVPerGNull) + 1
    }else {
      NumCV[["nullGene"]]  <- rep(NumCVPerGNull,NumCGNull)
    }
  }

  if(!grepl("a",simModel)){
    NumCG[["cau"]] = 0
    NumCV[["cau"]] = 0
    h2med = 0
    NumCGCau = 0
    NumCVPerGCau = 0
  }


  if(!grepl("b",simModel)){
    NumCG[["ple"]] = 0
    NumCV[["ple"]] = 0
    NumCGPle = 0
    NumCVPerGPle = 0
  }
  if(!grepl("c",simModel)){
    NumCG[["noteQTL"]] <- 0
    NumCV[["noteQTL"]] <- 0
    NumCVIG  = 0
    print("snps in intergenic region are not selected.")
    print("set cisOnly = TRUE instead")
    cisOnly = TRUE
  }

  if(!grepl("d",simModel)){
    NumCG[["nullGene"]] <- 0
    NumCV[["nullGene"]] <- 0
    NumCVPerGNull  = 0
    NumCGNull     = 0
  }


  if(!grepl("a",simModel) & !grepl("b",simModel) & grepl("c",simModel) ){
    print(paste0("Only SNPs in cis-region and causal model are used, so set ",
                 "h2cis = ",0))
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
  cisOnly = cisOnly
  mapSnpGeneAcrossChr <- readRDS(mapSnpGeneAcrossChrFile)
  mapLd2Snp <- readRDS(file = mapLd2SnpFile)
  smrPathBool = FALSE
  if(nchar(smrPath) !=0){smrPathBool = TRUE}
  genelist <- names(mapSnpGeneAcrossChr$mapGene2Snp)
  totalGeneNum <- length(genelist)
  mapGene2Snp <- mapSnpGeneAcrossChr$mapGene2Snp
  mapSnp2Gene <- mapSnpGeneAcrossChr$mapSnp2Gene
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
    if(geneOverlap == overType[1]){
      if(isVaryCaus){stop("not works for vary caus")}
      while(ge){
        gene <- sample(setdiff(genelist, unique(c(geneRemoved, genelistSlted)) ),1)
        snpInGene <- names(mapGene2Snp[[gene]])
        geneOverlapped <- c()
        for(singleSnp in snpInGene ){
          geneOverlapped <- unique(c(geneOverlapped,mapSnp2Gene[[singleSnp]]))
        }
        geneCandidate <- geneOverlapped[1]
        geneRemoved <- unique(c(geneRemoved, setdiff( geneOverlapped,geneCandidate) ) )
        if(length(mapGene2Snp[[geneCandidate]]) < NumCauVariants[ge] )
        {
          geneRemoved <- c(geneRemoved,geneCandidate)
        }else {
          ge <- ge - 1
          genelistSlted <- c(genelistSlted,geneCandidate)
          snps <- names(sample(mapGene2Snp[[geneCandidate]],NumCauVariants[ge]))
          snplistSetSlted <- unique( c(snplistSetSlted,  snps) )
          snplistMapslted[[geneCandidate]] <- snps
          snplistAllInSimGene <- c(snplistAllInSimGene, names(mapGene2Snp[[geneCandidate]]) )
          print(paste0(ge,"/",totalGeneNum ,"-th GENE [",geneCandidate,
                       "] was assigned [", NumCauVariants[ge] ,"] eqtls effect!"))
        }
      }
    }

    if(geneOverlap == overType[2]){
      if(isVaryCaus){stop("overlap situation not works for vary caus")}
      while(ge){
        gene <- sample(setdiff(genelist, unique(c(geneRemoved, genelistSlted)) ),1)
        snpInGene <- mapGene2Snp[[gene]]
        geneOverlapped <- mapSnp2Gene[[names(which(snpInGene == max(snpInGene ))[1])]]
        tmpGenes <- c()
        for(geneI in geneOverlapped){
          snpsTotal <- names(mapGene2Snp[[geneI]][mapGene2Snp[[geneI]] > 1] )
          if(length(snpsTotal) < NumCauVariants[ge] ){
            geneRemoved <- c(geneRemoved,geneI)
          } else {
            tmpGenes <- c(tmpGenes,geneI)
          }
        }
        if(length(tmpGenes) >= ge)
        {
          geneCandidate <- tmpGenes[1:ge]
          geneRemoved <- c(geneRemoved,tmpGenes)
          if(totalOverlapBool){
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
              snplistAllInSimGene <- c(snplistAllInSimGene, names(mapGene2Snp[[geneI]]) )
              print(paste0(ge,"/",totalGeneNum ,"-th GENE [",geneI,
                           "] was assigned [", NumCauVariants[ge] ,"] eqtls effect!"))
            }
          }
          ge = 0
        }else if (length(tmpGenes) < ge )
        {
          geneCandidate <- tmpGenes
          geneRemoved <- c(geneRemoved,tmpGenes)
          for(geneI in geneCandidate){
            snpsTotal <- names(mapGene2Snp[[geneI]][mapGene2Snp[[geneI]] > 1] )
            snps <- (sample(snpsTotal,NumCauVariants[ge]))
            snplistSetSlted <- unique( c(snplistSetSlted,  snps) )
            snplistMapslted[[geneI]] <- snps
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

    if(geneOverlap == overType[3]){
      while (ge){
        gene <- sample(setdiff(genelist, unique(c(geneRemoved, genelistSlted)) ),1)
        snpNumInGene <- length(mapGene2Snp[[gene]])
        if(snpNumInGene > NumCauVariants[ge]){
          genelistSlted <- c(genelistSlted,gene)
          geneRemoved <- c(geneRemoved,gene)
          snps <- names(sample(mapGene2Snp[[gene]],NumCauVariants[ge]))
          snplistSetSlted <- unique( c(snplistSetSlted,  snps) )
          snplistMapslted[[gene]] <- snps
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
    names(NumCV[[mod]]) <- rev(genelistSlted)
    genelist <- setdiff(genelist, genelistForModel[[mod]])
  }

  mapSnpOutGene <- mapSnpGeneAcrossChr$mapSnpOutGene[order(mapSnpGeneAcrossChr$mapSnpOutGene$POS),]
  snpFrommapSnpOutGene <- mapSnpOutGene[["SNP"]]
  totalMkNumOutGene <- 0
  if(grepl("c",simModel) ){
    mod = "noteQTL"
    print("###########################################")
    print(paste0(simModel, ": model is used for simulation") )
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

  mapSnpInGene <- mapSnpGeneAcrossChr$mapSnpInGene[order(mapSnpGeneAcrossChr$mapSnpInGene$POS),]
  mapGeneOrdered <- data.frame(GENE = names(mapSnpGeneAcrossChr$mapGene2Snp))
  if(!grepl("a",simModel) & !grepl("b",simModel) & !grepl("d",simModel) & grepl("c",simModel) ){
    snpIdxPos <- unique( mapSnpOutGene)
  }else {
    snpIdxPos <- unique(rbind(mapSnpInGene, mapSnpOutGene))
  }
  snpIdxPos  <- snpIdxPos[order(snpIdxPos$POS),]
  snpIdx  <- snpIdxPos[["SNP"]]


  if(grepl("a",simModel) | grepl("b",simModel) | grepl("d",simModel)){
    snplistAllInSimGene <- unique(snplistAllInSimGene)
    snplistAllInSimGene <- merge(mapSnpInGene,data.frame(SNP = snplistAllInSimGene),
                                 by = "SNP",sort=FALSE)[["SNP"]]
    if(cisOnly){ snpIdx = snplistAllInSimGene}
    geneIdx  <- unique(c(genelistForModel[["cau"]],genelistForModel[["ple"]],genelistForModel[["nullGene"]] ))
    geneIdx <- merge(mapGeneOrdered,data.frame(GENE = geneIdx),
                     by = "GENE",sort=FALSE)[["GENE"]]
    thetaTrue <- data.frame(theta = rep(0,length(geneIdx)),row.names = geneIdx)
    alphaTrue <- matrix(0, nrow= length(snplistAllInSimGene),
                        ncol= length(geneIdx))
    rownames(alphaTrue) <- snplistAllInSimGene
    colnames(alphaTrue) <- geneIdx
  }

  betaTrue <- data.frame(beta = rep(0,length(snpIdx)),row.names = snpIdx)
  cauSnpNumOutGene <- 0
  if(grepl("b",simModel)| grepl("c",simModel)){
    snpNotMedIdx <- unique( c(snplistSetForModel[["ple"]],snplistSetForModel[["noteQTL"]]) )
    snpNotMedIdx <-  merge(snpIdxPos,data.frame(SNP = snpNotMedIdx),by = "SNP",sort=FALSE) [["SNP"]]
    NumSnpNotMed <- length(snpNotMedIdx)
    sigmaForNotMed <- (h2snp - h2med)/NumSnpNotMed
    betaTrue[match(snpNotMedIdx,rownames(betaTrue)),] <- rnorm(NumSnpNotMed,0,sqrt( sigmaForNotMed) )
    cauSnpNumOutGene <- length(snplistSetForModel[["noteQTL"]])
  }
  if(grepl("a",simModel) | grepl("b",simModel) | grepl("d",simModel) ){
    if(NumCG[["cau"]] > 0){
      sigmaForGene <- h2med/(NumCG[["cau"]] * h2cis)
      thetaTrue[match(genelistForModel[["cau"]],rownames(thetaTrue)),] <- rnorm(NumCG[["cau"]],0,sqrt( sigmaForGene) )
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
        alphaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(alphaTrue)),gene] <-  rnorm(NumCV[["ple"]][[gene]],0,sqrt(sigmaForCis) )
        corBetaAlpah <- cor(betaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(betaTrue)),],
                            alphaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(alphaTrue)),gene] )
          while((corBetaAlpah)^2 > 1e-5){
          alphaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(alphaTrue)),gene] <-  rnorm(NumCV[["ple"]][[gene]],0,sqrt(sigmaForCis) )
          corBetaAlpah <- cor(betaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(betaTrue)),],
                              alphaTrue[match(snplistMapForModel[["ple"]][[gene]],rownames(alphaTrue)),gene] )
        }

      }
    }

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
  snpIdx <- snpIdx
  cauGeneNum <- 0
  totalGeneNum <- 0
  cauSnpNumInGene <- 0
  if(grepl("a",simModel)| grepl("b",simModel) | grepl("d",simModel) ){
    geneIdx <- colnames(alphaTrue)
    cauGeneNum <- cauGeneNum + length(genelistForModel[["cau"]])
    totalGeneNum <- totalGeneNum + length(geneIdx)
    cauSnpNumInGene <- length(snplistSetForModel[["cau"]]) + length(snplistSetForModel[["ple"]]) + length(snplistSetForModel[["nullGene"]])
  }
  cauSnpNum <- length(unique(as.vector(unlist(snplistSetForModel)))  )

  if(realGeno) {
    if(FALSE){
      ukb <-BEDMatrix::BEDMatrix(trainBfile)
      geno <- as.matrix(ukb)
      W <- BGData::preprocess(geno, center = FALSE, scale = FALSE,impute = TRUE)
      colnames(W)  <- gsub("(.+)_(.+)","\\1",colnames(W))
      selectedIndName <- sample(rownames(W),indNum)
      W <- W[match(selectedIndName,rownames(W)),]
      W <- apply(W, 2, function(x){ (x-mean(x))/sd(x) })
      print(paste0("REAL genotype is used ", "(N = ",indNum,")"))
    }
    if(TRUE){
      print(paste0("Train file: ", trainBfile))
      ukb <- BGData::as.BGData(BEDMatrix::BEDMatrix(trainBfile))
      geno <- as.matrix(ukb@geno)
      colnames(geno)  <- gsub("(.+)_(.+)","\\1",colnames(geno))
      geno <- geno[,match(snpIdx,colnames(geno))]
      selectedIndName <- sample(rownames(geno),indNum)
      geno <- geno[match(selectedIndName,rownames(geno)),]
      genoMap <- ukb@map[ ukb@map$snp_id %in% snpIdx ,]
      rownames(genoMap) = genoMap$snp_id
      pseudoYForGeno = data.frame(FID = gsub("(.+)_(.+)","\\1", rownames(geno) ),
                                  IID = gsub("(.+)_(.+)","\\2", rownames(geno) ),
                                  PAT = 0,MAT = 0,SEX = 0,PHENOTYPE = -9)
      rownames(pseudoYForGeno) = rownames(geno)
      genoInfo <- BGData::BGData(geno = geno,pheno = as.data.frame(pseudoYForGeno),map = genoMap)
      freq  <- BGData::summarize(X = BGData::geno(genoInfo) )
      genoInfo@map$freq_na      <- freq$freq_na
      genoInfo@map$sd           <- freq$sd
      genoInfo@map$allele_freq  <- freq$allele_freq
      genoInfo@map$na           <- colSums(is.na(geno))
      snp2pq                    <- 2 * freq$allele_freq * ( 1- freq$allele_freq)
      names(snp2pq)             <- snpIdx
      W <- BGData::preprocess(geno, center = TRUE, scale = TRUE,impute = TRUE)
      print(paste0("REAL ddd genotype is used ", "(N = ",indNum,")"))
    }

  }else {
    mkNum <- length(snpIdx)
    geno <- matrix(sample(0:2, indNum * mkNum, replace = TRUE), nrow = indNum )
    rownames(geno) <- paste0("ind-",1:indNum)
    selectedIndName = rownames(geno)
    colnames(geno) <- snpIdx
    pseudoX <- apply(geno, 2, function(x){ (x-mean(x))/sd(x) })
    print(paste0("PSEUDO genotype is used ", "(N = ",indNum,")"))
  }
  if(realGeno){
    X <-  W[match(selectedIndName,rownames(geno)),match(snpIdx,colnames(W))]
    remove(W)
  } else {
    X <- pseudoX
    remove(pseudoX)
    snp2pq = apply(X,2,var)
  }
  h2cis <- h2cis
  h2snp <- h2snp
  indNum <- dim(X)[1]; mkNum <-  dim(X)[2]
  betaTrue <- as.matrix(betaTrue)
  samlist <- rownames(X)
  if(grepl("a",simModel)| grepl("b",simModel) | grepl("d",simModel)){
    geneIdx <- colnames(alphaTrue);
    thetaTrue <- as.matrix(thetaTrue)
  }
  if(grepl("a",simModel) | grepl("b",simModel)| grepl("d",simModel)){
    eMatrix <- matrix(0,nrow = indNum ,ncol = length(geneIdx))
    colnames(eMatrix) <- geneIdx
    rownames(eMatrix) <- selectedIndName
    for (i in 1:length(geneIdx)) {
      eqtlG <- X[,snplistAllInSimGene] %*% alphaTrue[,i]
      eqtlVg <- var(eqtlG)
      eqtlVe <- (1-h2cis)/h2cis * eqtlVg
      eqtlExpression <- eqtlG + rnorm(indNum,0,sqrt(eqtlVe))
      eqtlExpression <- eqtlExpression - mean(eqtlExpression)
      eMatrix[,i] <- eqtlExpression
      print(paste0("Gene expression for gene ",i, " generated"))
    }
    varGene <- apply(eMatrix, 2, var)
    print(paste0("Gene expression matrix generated"))
  }
  gwasG <- 0
  if(grepl("b",simModel) | grepl("c",simModel) ){
  }
  if(grepl("a",simModel)){
    betaTrue[snplistSetForModel[["cau"]],] <- alphaTrue[snplistSetForModel[["cau"]],] %*% thetaTrue
  }
  gwasG <- gwasG +  X %*% betaTrue
  if (traitType == "a") {
    gwasVg <- var(gwasG)
    gwasVe <- (1- h2snp)/(h2snp) * (gwasVg)
    y <- gwasG + rnorm(indNum,0,sqrt(gwasVe))
    y <- (y - mean(y))

    print(paste0("Continuous trait vector generated"))
  } else if(traitType == "b") {
    pr = 1/(1+exp(- gwasG))
    y <- rbinom(indNum,1,pr)
    print(paste0("binary trait vector generated"))
  } else {

  }
  rownames(y) <- selectedIndName
  nGWAS = indNum
  if (traitType == "a") {
    bhat <- apply(X, 2, function(x){lm(y ~ x)$coef[2]})
    bhatSE <- apply(X, 2, function(x){coef(summary(lm(y ~ x)))[2,2]})
    bhatSMR = bhat /sqrt(snp2pq)
    bhatSESMR = bhatSE /sqrt(snp2pq)
    bhatSqrtScaleFactor = sqrt(1/(indNum * bhatSE*bhatSE + bhat*bhat))

    bhatPvalue <- apply(X, 2, function(x){coef(summary(lm(y ~ x)))[2,4]})
    for(i in 1:length(bhatPvalue)) {if(bhatPvalue[i] <= 1e-8) bhatPvalue[i] = 0}

  } else if(traitType == "b"){
    bhat <- apply(X, 2, function(x){glm(y ~ x,family = "binomial")$coef[2]})
    bhatSE <- apply(X, 2, function(x){coef(summary(glm(y ~ x,family = "binomial")))[2,2]})
  } else {

  }
  vary = 1
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
      snpID = intersect(names(mapGene2Snp[[geneIdx[i]]]), colnames(X) )
      XInGene <- X[,snpID]
      bInGene <-  apply(XInGene, 2, function(x){lm(eMatrix[,i] ~ x)$coef[2]})
      names(bInGene) <- snpID
      bInGeneSE <- apply(XInGene, 2, function(x){coef(summary(lm(eMatrix[,i] ~ x)))[2,2]})
      names(bInGeneSE) <- snpID

      AMarginSMR[snpID,geneIdx[i]] <- bInGene /sqrt(snp2pq[snpID])
      AMarginSESMR[snpID,geneIdx[i]] <- bInGeneSE /sqrt(snp2pq[snpID])
      AMarginSqrtScaleFactor[snpID,geneIdx[i]] <- sqrt(1/(indNum * bInGeneSE*bInGeneSE + bInGene*bInGene))
      AMarginSqrt2pq[snpID,geneIdx[i]]  = sqrt(snp2pq[snpID])
      bInGenePvalue <- apply(XInGene, 2, function(x){coef(summary(lm(eMatrix[,i] ~ x)))[2,4]})
      names(bInGenePvalue) <- snpID
      AMarginPvalue[snpID,geneIdx[i]] <- bInGenePvalue
      svarGene = median(snp2pq[snpID] * indNum * (indNum * AMarginSESMR[snpID,geneIdx[i]]^2 + AMarginSMR[snpID,geneIdx[i]]^2)/indNum)
      svarGene =1
      names(svarGene) = geneIdx[i]
      varGene = c(varGene,svarGene)
    }
    totalMkNumInGeneOverlap <- sum(colSums(AMargin != 0) )
    neQTLs = apply(eMatrix, 2, length)
  }

  geneSnpPosInfo <- mapSnpGeneAcrossChr$geneSnpPosInfo
  if(grepl("a",simModel) | grepl("b",simModel)| grepl("d",simModel) ){
    geneSnpPosInfo <- unique(geneSnpPosInfo[geneSnpPosInfo$GENE %in% geneIdx,])
    mapList <- ConstructMapSnp2Gene(geneSnpPosInfo,snplist = "SNP",genelist = "GENE",value="eQTLIdx",
                                             genePos = "GSTART", snpPos = "SNPPOS")
    mapSnp2GeneAcrossBlocks <- mapList$mapSnp2Gene
    mapGene2SnpAcrossBlocks <- mapList$mapGeneSNP2eQTL
    print("Gene overlap situation:")
    print( table(unlist(mapGene2SnpAcrossBlocks)))
  }

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
                 geneSnpPosInfo = geneSnpPosInfo,
                 mapLd2Snp = mapLd2Snp,
                 h2cis = h2cis,
                 h2snp = h2snp,
                 h2med = h2med,
                 traitType = traitType,
                 betaTrue = betaTrue,
                 bhat = bhat,
                 snp2pq = snp2pq,
                 bhatSqrtScaleFactor = bhatSqrtScaleFactor,
                 bhatSE = bhatSE,
                 bhatSMR = bhatSMR,
                 bhatSESMR = bhatSESMR,
                 vary = vary,
                 genoInfo = genoInfo,
                 nGWAS = nGWAS,
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
                 smrPathBool = smrPathBool,
                 smrIndFileSuffix = smrIndFileSuffix,
                 type = type
  )

  if(isIndLevelBool){
    simRes[["X"]] = X
    simRes[["y"]] = y
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
                       varGene = varGene,
                       mapSnp2GeneAcrossBlocks = mapSnp2GeneAcrossBlocks,
                       mapGene2SnpAcrossBlocks = mapGene2SnpAcrossBlocks
    )

    if(!isIndLevelBool) {
      simResGene[["Rgene"]] = Rgene
    } else {
      simRes[["Z"]] = X[,snplistAllInSimGene]
      simRes[["eMatrix"]] = eMatrix
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
  print(paste0("Date for simulation generated here !!"))
  return (simRes)

}
