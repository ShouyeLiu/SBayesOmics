downsampleN = function(geneSummaryFile,geneLDFile,targetN){
  xqtl = data.table::fread(paste0(geneSummaryFile,".q.query.gz"))
  genelist = unique(xqtl$GeneID)
  rawN = unique(xqtl$N)
  ### read low-rank LD matrix
  lowLDGene =  readxQTLEigen( geneLDFile, eigenType = "gene", cutThresh = 0.995)
  Ugene = lowLDGene$UBlocks
  lambdaGene = lowLDGene$lambdaBlocks
  keptGenes = names(Ugene)
  xqtlPseudo <- c()
  qtlTrue<- c()
  qtlP1K <- c()
  qtlP10K <- c()
  genelist = unique(xqtl$GeneID)
  rawN = unique(xqtl$N)
  for(genei in genelist ){
    xqtlg = xqtl[GeneID == genei]
    xqtlg[,snp2pq:= 2*A1Freq *(1-A1Freq)]
    xqtlg[, D:=2 * A1Freq * (1 - A1Freq) * N]
    xqtlg[, b:= BETA]
    xqtlg[,se:=SE]
    xqtlg[, varps := D*(N * se^2 + b^2)/N]
    xqtlsnp = xqtlg$SNPID
    xqtlEff = xqtlg$BETA * sqrt(xqtlg$snp2pq)
    xqtlSE = xqtlg$SE * sqrt(xqtlg$snp2pq)
    # ## true sample size
    # N = unique(xqtlg$N)
    # sqtlTrue <- data.table(GeneID = genei,SNPID = xqtlsnp,beta = xqtlEff,se = xqtlSE,Type = "True_10K")
    # qtlTrue = rbind(qtlTrue,sqtlTrue)
    ## ld
    ld = Ugene[[genei]] %*% diag(lambdaGene[[genei]]) %*% t(Ugene[[genei]])
    diag(ld) = diag(ld) + 0.2
    ## original sample size
    sqtlP10K = xqtlg %>%
      mutate(GeneID = genei, xqtlEff = BETA * sqrt(snp2pq), xqtlSE = SE * sqrt(snp2pq), n_trn = targetN) %>%
      mutate( betaPse = as.numeric(xqtlEff + sqrt(1/n_trn - 1/N) * Ugene[[genei]] %*% diag(sqrt(lambdaGene[[genei]])) %*% rnorm(length(lambdaGene[[genei]]),0,1) ),
              sePse = sqrt((1 - betaPse* betaPse))/sqrt(n_trn),
              betaR = xqtlEff,
              seR = sqrt(N/n_trn) * xqtlSE,
              NOri = N,
              N = n_trn
      ) %>%
      select(GeneID,SNPID,betaPse,sePse,betaR,seR,N,NOri)

    # pseudoXqtl = xqtlEff + sqrt(1/n_trn - 1/N) * Ugene[[genei]] %*% diag(sqrt(lambdaGene[[genei]])) %*% rnorm(length(lambdaGene[[genei]]),0,1)
    # pseudoXqtl = pseudoXqtl[,1]
    # seR = sqrt(N/n_trn) * xqtlSE
    # sePse = sqrt((1 - pseudoXqtl* pseudoXqtl))/sqrt(n_trn)
    # sqtlP10K <- data.table(GeneID = genei,SNPID = xqtlsnp,
    #                        betaPse = pseudoXqtl, sePse = sePse,
    #                        betaR = xqtlEff,seR = seR,
    #                        N = n_trn)
    qtlP10K = rbind(qtlP10K,sqtlP10K)
  }
  return(qtlP10K)
}
