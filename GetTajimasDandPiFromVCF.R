## GetTajimasDandPiFromVCF
# Authors: Aneil F. Agrawal and Amardeep Singh

# This script will take a formatted VCF file and provide gene level estimates of TajimasD and PiN/piS
# Be sure to check the example file for the correct formatting of the VCF file for the script. Some
# wrangling will be required.

rm(list=ls())

################
## USER INPUTS##
################
# Users should only have to edit the folling three objects, 
# including population name, and file paths input VCF and output file
# Read in VCF file
population = "DPGP3"
inputFileName = paste("/plas1/amardeep.singh/RNA.Seq.Data/piNpiS.Analysis/4fold0fold.parsedVCFs/4Fold0Fold_withFBgn/Chr2L_4fold_FBgn_formatted.vcf")
#inputFileName = paste("/plas1/amardeep.singh/RNA.Seq.Data/DPGP3.Genomic.Data/dpgp3_sequences/test.vcf")
outputFileName = "/plas1/amardeep.singh/RNA.Seq.Data/DGRP.Genomic.Data/DGRP_4Fold_0Fold_Output/Chr2L_4fold_"
######################################
######################################

vcf_file = read.table(inputFileName, sep = "\t", header = FALSE)

if (population == "DPGP3"){
  vcf_file = vcf_file[nchar(as.character(vcf_file$V4)) < 4,]
} else if (population == "DGRP"){
  vcf_file = vcf_file[nchar(as.character(vcf_file$V5)) < 4,]
}
# Specifications of VCF file
first_col = 6   # First column that contains genotypic information
last_col = 210  # Last column that contains genotypic information
geneID_col = 1
minAllowableSampleSize = 75

GetTajimasD<-function(sfs, L, n){  ## following p. 303 of Walsh and Lynch book.
  S<- sum(sfs)/L
  an <- sum(sapply(1:(n-1), function(x) 1/x))  ## p. 301
  bn <-  sum(sapply(1:(n-1), function(x) 1/(x^2))) ## p. 301
  thetaW = S/an
  betaD <- (1/(an^2 + bn))*( 2*(n^2 + n + 3)/(9*n*(n-1)) - (n+2)/(an * n) + bn/(an^2)  ) ## p. 304
  alphaD<- (1/an)*( (n+1)/(3*(n-1)) - (1/an) ) - betaD ## p. 304

  PI = sum(sapply(1:length(sfs), function(x) sfs[x]*2*(x/n)*((n - x)/(n - 1)))) / L

  TajD = (PI - thetaW)/sqrt(alphaD*S + betaD*(S^2))
  return(TajD)
}

# exampleSFS <- c(25, 15, 9, 4, 4, 3, 2, 2, 0, 2, 1, 0, 0, 1, 0)
# GetTajimasD(exampleSFS, 15)



### The function below does the following
### 1) Excludes sites with fewer than some threshold level of samples (used 100)
### and rescales L accordingly, i.e., if 5% of sites are dropped then L_adjusted = 0.95*L.
### 2) Of retained sites, get minimum n across sites (nmin)
### 3) Subsample:For each retained site, sample without replacement nmin times
### 4) Make SFS from this this subsampled data and calculate Tajima's D
### 5) Repeat steps 3-4 100 times and returns the median D and the corresponding S


GetTajimasDFromVCF<-function(vcf, geneID){
      vcf = vcf[vcf[,geneID_col] == geneID, ]
  vcf = vcf[, c(first_col:last_col)]
    L = nrow(vcf)
    sampleSizePerSite<- sapply(1:(dim(vcf)[1]), function(x) sum(!is.na(vcf[x,])))
      minAllowableSampleSize<-minAllowableSampleSize
      numExcludedSites<-sum(sampleSizePerSite < minAllowableSampleSize)

        vcf.Filtered<-vcf[sampleSizePerSite>=minAllowableSampleSize, ]
        sampleSizePerSite.UsableSites<- sapply(1:(dim(vcf.Filtered)[1]), function(x) sum(!is.na(vcf.Filtered[x,])))
          numNonRefSamplesPerSite.UsableSites<- sapply(1:(dim(vcf.Filtered)[1]), function(x) sum(vcf.Filtered[x,], na.rm = TRUE))
          minSampleSize<-min(sampleSizePerSite.UsableSites)

            nSubSamples = 100; listOfTajDValues<-rep(NA, nSubSamples); listOfSValues<-rep(NA, nSubSamples)
          for(i in 1:nSubSamples){
                  nNonRefSNPsPerSite<- sapply(1:(dim(vcf.Filtered)[1]), function(x) sum(sample(c(rep(0, sampleSizePerSite.UsableSites[x]- numNonRefSamplesPerSite.UsableSites[x]),
                                                                                                                                                             rep(1, numNonRefSamplesPerSite.UsableSites[x])), size = minSampleSize, replace = FALSE)))
              sfs<-rep(NA, dim(vcf.Filtered)[2])
                  for(j in 1:length(sfs)) sfs[j] = sum(nNonRefSNPsPerSite == j)
                  # if(i == 1) cat("sfs = ", sfs, "\n")
                  listOfTajDValues[i]<-GetTajimasD(sfs, L*(1- (numExcludedSites/(dim(vcf)[1]))), minSampleSize)
                      listOfSValues[i]<-sum(sfs)

          }

            returnThisTajD<-median(listOfTajDValues)
            indexOfMedian<-which(abs(listOfTajDValues - median(listOfTajDValues)) == min(abs(listOfTajDValues - median(listOfTajDValues))))
              returnThisS<- mean(listOfSValues[indexOfMedian])
              return((c(returnThisTajD, returnThisS)))

}

TajimasD_wrapper = function(vcf, geneID_col){
      UniqueGeneIDs = unique(as.character(vcf[,geneID_col]))

  TajimasD = (t(sapply(1:length(UniqueGeneIDs), function(x) GetTajimasDFromVCF(vcf,UniqueGeneIDs[x]))))
    out = (cbind(UniqueGeneIDs,TajimasD))
    colnames(out) = c("geneID", "D", "S")
      return(out)

}





GetSummarizedInfoFromVCFthatIncludesInvariantSites.V1<-function(vcf){
      vcf = vcf[, c(first_col:last_col)]

  nSites.Starting = (dim(vcf)[1])
    sampleSizePerSite <- sapply(1:(dim(vcf)[1]), function(x) sum(!is.na(vcf[x,])))
    minAllowableSampleSize<-minAllowableSampleSize
      #numExcludedSites<-sum(sampleSizePerSite < minAllowableSampleSize)

      vcf.Filtered<-vcf[sampleSizePerSite>=minAllowableSampleSize, ]


      sampleSizePerSite.UsableSites <- sapply(1:(dim(vcf.Filtered)[1]), function(x) sum(!is.na(vcf.Filtered[x,])))
        numNonRefSamplesPerSite.UsableSites<- sapply(1:(dim(vcf.Filtered)[1]), function(x) sum(vcf.Filtered[x,], na.rm = TRUE))
        minSampleSize<-min(sampleSizePerSite.UsableSites)

          nSites.For.pi.calculation = dim(vcf.Filtered)[1]

          pPerSite = numNonRefSamplesPerSite.UsableSites/sampleSizePerSite.UsableSites
            pi.VecPerSite.AllUseableSitesAndSamples = sapply(1:length(pPerSite), function(x) 2*pPerSite[x]*(1-pPerSite[x])  )
            pi.AllUseableSitesAndSamples = mean(pi.VecPerSite.AllUseableSitesAndSamples)
              VariableSitesPresent = ifelse(sum(pi.AllUseableSitesAndSamples != 0) > 0,  TRUE, FALSE)
              InvariantSitesPresent = ifelse(sum(pi.AllUseableSitesAndSamples == 0) > 0,  TRUE, FALSE)
                minSampleSize.UseableVariantSites = ifelse(VariableSitesPresent, min(sampleSizePerSite.UsableSites[pi.AllUseableSitesAndSamples != 0]), 0)
                minSampleSize.UseableInvariantSites = ifelse(InvariantSitesPresent, min(sampleSizePerSite.UsableSites[pi.AllUseableSitesAndSamples == 0]), 0)



                if(VariableSitesPresent){
                        vcf.FilteredMore<-vcf.Filtered[sampleSizePerSite.UsableSites>=minSampleSize.UseableVariantSites, ]
                    sampleSizePerSite.UsableSites2<- sapply(1:(dim(vcf.FilteredMore)[1]), function(x) sum(!is.na(vcf.FilteredMore[x,])))
                        numNonRefSamplesPerSite.UsableSites2<- sapply(1:(dim(vcf.FilteredMore)[1]), function(x) sum(vcf.FilteredMore[x,], na.rm = TRUE))

                    downsample.for.one.sfs<-function(this.vcf, this.sampleSizePerSite, this.numNonRefSamplesPerSite, this.minSampleSize){
                              nNonRefSNPsPerSite<- sapply(1:(dim(this.vcf)[1]), function(x) sum(sample(c(rep(0, this.sampleSizePerSite[x]- this.numNonRefSamplesPerSite[x]),
                                                                                                                                                                                          rep(1, this.numNonRefSamplesPerSite[x])), size = this.minSampleSize, replace = FALSE)))
                          sampled.sfs<-rep(NA, dim(this.vcf)[2] - 1)
                                for(j in 1:length(sampled.sfs)) sampled.sfs[j] = sum(nNonRefSNPsPerSite == j)
                                nInvariantSites = sum(nNonRefSNPsPerSite == 0 | nNonRefSNPsPerSite == this.minSampleSize)
                                      # nInvariantSites0 = sum(nNonRefSNPsPerSite == 0 )
                                      # nInvariantSites1 = sum(nNonRefSNPsPerSite == this.minSampleSize)
                                      # cat("nNonRefSNPsPerSite = ", nNonRefSNPsPerSite, "\n")
                                      # cat("c(nInvariantSites, nInvariantSites0, nInvariantSites1) = ", c(nInvariantSites, nInvariantSites0, nInvariantSites1), "\n")
                                      return(c(nInvariantSites, sampled.sfs))

                    }

                  one.sample.nInvariantSites.And.sfs<- downsample.for.one.sfs(vcf.FilteredMore, sampleSizePerSite.UsableSites2, numNonRefSamplesPerSite.UsableSites2, minSampleSize.UseableVariantSites)
                      # cat(sampleSizePerSite.UsableSites2, "\n", numNonRefSamplesPerSite.UsableSites2, "\n")
                      nSubSamples = 100; listOfTajDValues<-rep(NA, nSubSamples); listOfSValues<-rep(NA, nSubSamples)
                  for(i in 1:nSubSamples){
                            sfsAndnInvariant<-downsample.for.one.sfs(vcf.FilteredMore, sampleSizePerSite.UsableSites2, numNonRefSamplesPerSite.UsableSites2, minSampleSize.UseableVariantSites)
                        this.L = sum(sfsAndnInvariant)
                              # if(i == 1) cat("sfs = ", sfsAndnInvariant[2:(minSampleSize.UseableVariantSites)], "\n")
                              listOfTajDValues[i]<-GetTajimasD(sfsAndnInvariant[2:(minSampleSize.UseableVariantSites)], this.L, minSampleSize.UseableVariantSites)
                              listOfSValues[i]<-sum(sfsAndnInvariant[2:length(sfsAndnInvariant)])

                  }

                    returnThisTajD<-median(listOfTajDValues)
                      indexOfMedian<-which(abs(listOfTajDValues - median(listOfTajDValues)) == min(abs(listOfTajDValues - median(listOfTajDValues))))
                      returnThisS<- mean(listOfSValues[indexOfMedian])
                        }



                if(!VariableSitesPresent){
                        one.sample.nInvariantSites.And.sfs<-rep(NA, dim(vcf.Filtered)[2])
                    returnThisTajD<-NA
                        returnThisS<-NA

                }
                  return(c(nSites.Starting, nSites.For.pi.calculation, pi.AllUseableSitesAndSamples, minSampleSize.UseableVariantSites,
                                      returnThisTajD, returnThisS, one.sample.nInvariantSites.And.sfs))
}




## This function is just a wrapper so that you can apply the operation of nterest
## to a whole list of genes it and won't crash midway if one gene doesn't work for some reason (other than the vcf isn't two dimensional)
GetSummarizedInfoFromVCFthatIncludesInvariantSites<-function(vcf, geneID){
      #print(c(geneID))
      vcf = vcf[vcf[,geneID_col] == geneID, ]
  x<-rep(NA, 6 + dim(vcf)[2])
    worked = FALSE
    try({x<-GetSummarizedInfoFromVCFthatIncludesInvariantSites.V1(vcf); worked = TRUE})
      if(!worked) cat("Houston, we have a problem! GetSummarizedInfoFromVCFthatIncludesInvariantSites.V1 didn't work.")
      return(x)

}

VCF_summary_wrapper = function(vcf, geneID_column){
      UniqueGeneIDs = unique(as.character(vcf[,geneID_column]))
  results = as.data.frame(matrix(NA, nrow = length(UniqueGeneIDs), ncol = 7 + (last_col - first_col)))
    #results = data.frame(NA, test = length(UniqueGeneIDs), ncol = 7 + (last_col - first_col + 1))
    results$geneID = UniqueGeneIDs
  for (i in 1:length(UniqueGeneIDs)){
          print(UniqueGeneIDs[i])
      x = rep(NA, dim(results)[2] - 1)
          try({x = GetSummarizedInfoFromVCFthatIncludesInvariantSites(vcf,as.character(UniqueGeneIDs[i])) })
          results[i, 1:(dim(results)[2] -1)] = x
              write.table(results, file = paste0(outputFileName, "Summary.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

              #write.table(results, file = paste0(outputFileName, "Summary_April132022.txt", quote = F, row.names = F, col.names = T, sep = "\t"))

  }
    #VCF_summary = (t(sapply(1:length(UniqueGeneIDs), function(x) GetSummarizedInfoFromVCFthatIncludesInvariantSites(vcf,UniqueGeneIDs[x]))))
    #out = (cbind(UniqueGeneIDs,TajimasD))
    #colnames(out) = c("geneID", "D", "S")
    #return(out)
    return(results)

}

output_df = VCF_summary_wrapper(vcf_file, geneID_col)

write.table(output_df, file = paste0(outputFileName, "Summary.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

rm(list=ls())
