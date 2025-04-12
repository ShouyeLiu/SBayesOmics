#' SBayesOmics - ConstructMapSnp2Gene
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


ConstructMapSnp2Gene <- function(data,snplist = "SNP",genelist = "GENE",value="eQTL",
                                 genePos = "GSTART", snpPos = "SNPPOS")
{
  ## 1. construc gene and snp ID dataframe
  snpID <- unique(data.frame(snpID = data[[snplist]],snpPosition = data[[snpPos]]) )
  geneID <- unique(data.frame(geneID = data[[genelist]],genePosition = data[[genePos]]) )

  ## 2. order gene by gene start or end point
  geneID <- geneID[order(geneID$genePosition),]
  geneIdx <- geneID$geneID
  ## order snp by snp position
  snpID <- snpID[order(snpID$snpPosition),]
  snpIdx <- snpID$snpID

  mapSnp2Gene <- c()
  ## 3. begin to construct map from snp to gene.
  for (i in 1:length(snpIdx)){
    singleSnpSet <- data[which(data[[snplist]] == snpIdx[i]),]
    ## order gene by gene position
    singleSnpSet <- singleSnpSet [order(singleSnpSet[[genePos]]),][[genelist]]
    # singleSnpSet = singleSnpSet %>%  arrange(GCHR, !!sym(genePos)) %>% filter()
    mapSnp2Gene[[snpIdx[i]]]  <- singleSnpSet
    }

  ## 4. map gene-snp to eqtl to store eQTL values.
  mapGeneSNP2eQTL <- c()
  genelistNotOverlap <- c() ## for gene
  for(i in 1:length(geneIdx)){
    singleGene <- data[which(data[[genelist]] == geneIdx[i]),]
    ## order snp by snpPos
    singleGene <- singleGene[order(singleGene[[snpPos]]),]
    snpNumInGene <- c()
    for(j in 1:dim(singleGene)[1]){
      singleSnp <- singleGene[j,][[snplist]]
      snpNumInGene <- c(snpNumInGene,length( mapSnp2Gene[[singleSnp]] ))
    }
    if(max(snpNumInGene) < 2){genelistNotOverlap <- c(genelistNotOverlap,geneIdx[i])}
    mapGeneSNP2eQTL[[geneIdx[i]]] <- snpNumInGene
    names(mapGeneSNP2eQTL[[geneIdx[i]]]) <- singleGene[[snplist]]
  }
  ## construct snp-gene eQTL matrix
  # ## 3. Initial coefficent matrix A
  # eQTLMatrix <- matrix(0, nrow= length(snpIdx), ncol= length(geneIdx))
  # rownames(eQTLMatrix) <- snpIdx
  # colnames(eQTLMatrix) <- geneIdx
  # ## 4. begin to construct matrix
  # for (i in 1:length(snpIdx))
  # {
  #   snpSetSub <- data[which(data[[snplist]] == snpIdx[i]),]
  #   for (j in 1:dim(snpSetSub)[1])
  #   {
  #     colIdx <- match(snpSetSub[[genelist]][j],geneIdx)
  #     ## eQTLMatrix[i,colIdx] <- snpSetSub$eQTL[j]
  #     eQTLMatrix[i,colIdx] <- snpSetSub[[value]][j]
  #   }
  # }
  geneName = names(mapGeneSNP2eQTL)

  out <- list(mapSnp2Gene = mapSnp2Gene,
              mapGeneSNP2eQTL = mapGeneSNP2eQTL,
              # eQTLMatrix = eQTLMatrix,
              geneName = geneName,
              genelistNotOverlap = genelistNotOverlap)
  return(out)
}
