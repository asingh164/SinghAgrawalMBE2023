### Here we use "SDIU - Sexual Dimorphism in Isoform Usage" synonymously with "SSS - Sex-specific splicing"


SBGEandSSS.MBE = read.csv(file = "/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/SexDifferences_IsoFormUsage/SBGEandSSSdataForMBE.csv", header = T)

main.df = SBGEandSSS.MBE
main.df$SBGEcat.body.Osada = factor(main.df$SBGEcat.body.Osada, levels = c("extFB", "sFB", "FB", "UB", "MB", "sMB", "extMB"))
main.df$SBGEcat.head.Osada = factor(main.df$SBGEcat.head.Osada, levels = c("extFB", "sFB", "FB", "UB", "MB", "sMB", "extMB"))

# data frame columns
# 1) geneID - FlyBase gene ID
# 2) chrm - chromosomal location
# 3) Is.A - indicator of autosomal gene
# 4) Is.X - indicator of X-linked gene (Is.X = 1 - Is.A)
# 5) Whole.SBGE.Osada - log2 sex effect from DESeq2 analysis of whole body samples of Osada et al. data set
# 6) padj.Whole.SBGE.Osada - adjusted pvalue of sex effect
# 7) SBGEcat.body.Osada - binned by SBGE:  extreme female-bias "extFB", strong female-bias "sFB", modest female-bias "FB", unbiased "UB", modest male-bias "MB", strong male-bias "sMB", extreme male-bias "estMB"
# 8) SDIU.body.geneWisePadj - gene-level adjusted pvalue for sex differences in isoform usage from JunctionSeq of whole body samples of Osada et al. data set
# 9) SDIU.body.sig - indicator whether SDIU.body.geneWisePadj < 0.05 or not
# 10) SDIU.body.p.strength.cat - categorized strength of evidence for SDIU:  "0" -> p >= 0.2; "1" -> 0.2 > p >= 0.01; "1" -> p < 0.01 ; p = SDIU.body.geneWisePadj
# 11-16) As 5-10 but with respect to head
# 17) factor denoting combined categorization with respect to both SDIU and SBGE for body expression
# 18) factor denoting combined categorization with respect to both SDIU and SBGE for head expression
# 19) intersexual genetic correlation in expression (rmf), calculated from Huang et al. 2015 data
# 20) measure of the genetic variation in female expression (variation among DGRP lines) as calculated by Huang et al. 2015
# 21) significance of variance component above
# 22-23) As 20-21 but with respect to male expression
# 24) Indicator of significant genetic variation (p < 0.05) in both sexes 
# 25) Indicator of significant genetic variation (p < 0.05) in either sex
# 26) Number of synonymous polymorphisms taken from Fraisse et al. (2019)
# 27) Number of nonsynonymous polymorphisms taken from Fraisse et al. (2019)
# 28) Number of synonymous interspecific differences taken from Fraisse et al. (2019)
# 29) Number of nonsynonymous interspecific differences taken from Fraisse et al. (2019)
# 30) Direction of Selection: DoS = Dn/(Dn + Ds) - Pn/(Pn + Ps)
# 31) recomb taken from Fraisse et al. (2019)
# 32) gene length calculated as total exon length using featureCounts
# 33) Divergence part of DoS, i.e., Dn/(Dn + Ds)
# 34) Polymorphism part of DoS, i.e., Dn/(Dn + Ds)
# 35) number of 4-fold degenerate sites used for calculation of pi
# 36) pi for 0-fold degenerate sites
# 37) pi for 4-fold degenerate sites
# 38) Tajima's D based on 0-fold degenerate sites (not used)
# 39) Tajima's D based on 4-fold degenerate sites 
# 40) piNfract = piN/(piN + piS)
# 41) log ofexpression across larval and adjult tissues (based on Fly Atlas 2 data; see Methods)
# 42) tissue specificity across larval and adult tissues (based on Fly Atlas 2 data; see Methods)
# 43-45) first 3 principal components of among-tissue expression variation in males (based on Fly Atlas 2 data; see Methods)
# 46-48) first 3 principal components of among-tissue expression variation in females (based on Fly Atlas 2 data; see Methods)



#############
library(lmPerm)
library(car)
library(boot)
#####

##############################################

## does one permutation to give difference
## in mean value between 2 groups
## given a list of values from 2 groups (no "NA" allowed)
## where the first n.group1 entries are values from group 1
get.one.perm.diff<- function(x.vec, n.group1){
  perm.vec = sample(x.vec, size = length(x.vec), replace = FALSE)
  return(mean(perm.vec[1:n.group1]) - mean(perm.vec[(n.group1 +1):length(perm.vec)]))
}


# summary function given vector of numeric values; returns
# 1) sample size
# 2) mean
# 3-4) bootstrap 95% CI bounds
# 5) SE of mean
# returns NAs for 2-5 if sample size is less than 10
summary.func.A<-function(z){
  z = na.omit(z)
  if(length(z) < 10) x = c(length(z), rep(NA, 4))
  if(length(z) >= 10){
    nboot = 10000
    bootvals.vec = rep(NA, nboot)
    for(i in 1:nboot) {
      boot.sample = sample(z, size = length(z), replace = TRUE)
      bootvals.vec[i] = mean(boot.sample)
    }
    x = as.numeric(c(length(z), mean(z), quantile(bootvals.vec, probs = c(0.025, 0.975)), sqrt(var(z)/length(z))))
  }
  return(x)
}


### Returns a df the meaning of the columns depends on the summary function used
### This is used provide a summary of a (continuous) variable (e.g., Tajima's D) with respect
### to another factor (e.g., SBGE)
MakeSummaryOfVarByFactor<-function(varName, factorName, df, summaryfunc = summary.func.A){
  factor.col.num = which(names(df) == factorName)
  if(is.factor(df[,factorName])) myfactor = df[,factorName]
  if(!is.factor(df[,factorName])){
    myfactor = as.factor(df[,factorName])
  } 
  n = nlevels(myfactor)
  if(n > 30) cat("\n\n Are you sure you specified a factor \n or something that can be sensibly converted to a factor? \n\n")
  
  binNames = levels(myfactor)
  for(i in 1:n) {
    z = df[myfactor == binNames[i] ,varName]
    y = summaryfunc(z)
    if(i == 1){
      m = matrix(nrow = n, ncol = length(y))
      m[1,] = y
    }
    if(i > 1) m[i,] = y
  }
  
  m.2 = rbind(summaryfunc(df[,varName]), m)
  binNames.2 = c("All", binNames)
  summary.df = as.data.frame(m.2)
  summary.df$BinName = factor(binNames.2, levels = binNames.2)
  # summary.df = as.data.frame(m)
  # summary.df$BinName = factor(binNames, levels = binNames)
  summary.df = summary.df[,c(dim(summary.df)[2], 1:(dim(summary.df)[2]-1))]
  return(summary.df)
}

## Returns a summary table for a given variable
## (1) Category of gene (wrt to SBGE and SDIU)
## (2) number of genes in category
## (3) mean
## (4-5) bootstrap 95% CIs
MakeSummaryTableForSDIUxSBGE.forSpecificVariable<-function(this.varName, this.SBGEcat.name, this.SDIUxSBGE.name, df){
  if(this.SBGEcat.name %in% c("SBGEcat.body.Osada", "SBGEcat.body.FA2", "SBGEcat.body.Mishra")) SDIUcol.name = "SDIU.body.sig"
  if(this.SBGEcat.name %in% c("SBGEcat.head.Osada", "SBGEcat.head.FA2", "SBGEcat.head.Mishra")) SDIUcol.name = "SDIU.head.sig"
  summary.Z.by.SDIU = MakeSummaryOfVarByFactor(this.varName, SDIUcol.name, df)
  levels(summary.Z.by.SDIU$BinName)<-c("All", "NotSDIU", "SDIU")
  pre.summary.Z.by.SBGE = MakeSummaryOfVarByFactor(this.varName, this.SBGEcat.name, df)
  summary.Z.by.SBGE = pre.summary.Z.by.SBGE[2:nrow(pre.summary.Z.by.SBGE),]
  pre.summary.Z.by.SDIUxSBGE = MakeSummaryOfVarByFactor(this.varName, this.SDIUxSBGE.name, df)
  summary.Z.by.SDIUxSBGE = pre.summary.Z.by.SDIUxSBGE[2:nrow(pre.summary.Z.by.SDIUxSBGE),]
  combined.summary.Z.by.SDIUxSBGE = rbind(summary.Z.by.SDIU[c(1,3,2),], summary.Z.by.SDIUxSBGE, summary.Z.by.SBGE)
  combined.summary.Z.by.SDIUxSBGE$BinName = as.character(combined.summary.Z.by.SDIUxSBGE$BinName)
  ordered.factor.names = c("All", "NotSDIU", "SDIU",
                           "extFB", "NotSDIU_extFB", "SDIU_extFB",
                           "sFB", "NotSDIU_sFB", "SDIU_sFB",
                           "FB", "NotSDIU_FB", "SDIU_FB",
                           "UB", "NotSDIU_UB", "SDIU_UB",
                           "MB", "NotSDIU_MB", "SDIU_MB",
                           "sMB", "NotSDIU_sMB", "SDIU_sMB",
                           "extMB", "NotSDIU_extMB", "SDIU_extMB")
  combined.summary.Z.by.SDIUxSBGE$BinName = factor(combined.summary.Z.by.SDIUxSBGE$BinName, levels = ordered.factor.names)
  sorted.combined.summary.Z.by.SDIUxSBGE = combined.summary.Z.by.SDIUxSBGE[order(combined.summary.Z.by.SDIUxSBGE$BinName),]
  return(sorted.combined.summary.Z.by.SDIUxSBGE)
}




Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body<-function(summaryTable, y.varLabel, this.ylim = c(-1,1)){
  vec.col.for.fig.body = rep(c("black", "grey", "purple"), 8)
  xvals.for.fig.body = as.vector(sapply(1:8, function(x) x + 0.2*(0:2)))
  plot(xvals.for.fig.body,
       summaryTable[,3], col = vec.col.for.fig.body,
       xlim = c(0.9, 8.4),
       ylim = this.ylim,
       xlab = "SBGE category", ylab = y.varLabel, pch = 19, xaxt = 'n')
  abline(v = 1.7, lwd = 3)
  abline(v = 2.7, lwd = 0.5); abline(v = 3.7, lwd = 0.5); abline(v = 4.7, lwd = 0.5)
  abline(v = 5.7, lwd = 0.5); abline(v = 6.7, lwd = 0.5); abline(v = 7.7, lwd = 0.5)
  arrows(x0=xvals.for.fig.body, y0=summaryTable[,4],
         x1 = xvals.for.fig.body, y1=summaryTable[,5], 
         code=3, angle = 90, length = 0.05,
         col = vec.col.for.fig.body, lwd=2)
}


Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Head<-function(summaryTable, y.varLabel, this.ylim = c(-1,1)){
  summaryTable = summaryTable[c(1:3, 10:18),]
  vec.col.for.fig.head = rep(c("black",  "grey", "purple"), 4)
  xvals.for.fig.head = as.vector(sapply(1:4, function(x) x + 0.2*(0:2)))
  plot(xvals.for.fig.head,
       summaryTable[,3], col = vec.col.for.fig.head,
       xlim = c(0.9, 4.4),
       ylim = this.ylim,
       xlab = "SBGE category", ylab = y.varLabel, pch = 19, xaxt = 'n')
  abline(v = 1.7, lwd = 3)
  abline(v = 2.7, lwd = 0.5); abline(v = 3.7, lwd = 0.5); 
  arrows(x0=xvals.for.fig.head, y0=summaryTable[,4],
         x1 = xvals.for.fig.head, y1=summaryTable[,5], 
         code=3, angle = 90, length = 0.05,
         col = vec.col.for.fig.head, lwd=2)
}


##### Counts for Table 1

summary(main.df$SDIU.head.sig | main.df$SDIU.body.sig)
summary(!is.na(main.df$SDIU.head.sig) | !is.na(main.df$SDIU.body.sig))

Tested.SBGEinBoth.vec = (!is.na(main.df$Whole.SBGE.Osada)) & (!is.na(main.df$Head.SBGE.Osada))
Tested.SDIUinBoth.vec = (!is.na(main.df$SDIU.body.sig)) & (!is.na(main.df$SDIU.head.sig))
Tested.SBGEandSDIUinBoth.vec = Tested.SBGEinBoth.vec & Tested.SDIUinBoth.vec
sum(Tested.SBGEinBoth.vec)
sum(Tested.SBGEandSDIUinBoth.vec)

geneID.for.Tested.SBGEandSDIUinBoth = main.df$geneID[Tested.SBGEandSDIUinBoth.vec]

SDIUinBoth.vec = (main.df$SDIU.body.sig) & (main.df$SDIU.head.sig)
Not.SDIUinBoth.vec = (!main.df$SDIU.body.sig) & (!main.df$SDIU.head.sig)
Not.Tested.Both.vec = (is.na(main.df$SDIU.body.sig)) & (is.na(main.df$SDIU.head.sig))

SDIUinBody.Not.SDIUinHead.vec = (main.df$SDIU.body.sig) & (!main.df$SDIU.head.sig)
SDIUinBody.NotTestedHead.vec = (main.df$SDIU.body.sig) & (is.na(main.df$SDIU.head.sig))
Not.SDIUinBody.SDIUinHead.vec = (!main.df$SDIU.body.sig) & (main.df$SDIU.head.sig)
Not.SDIUinBody.NotTestedHead.vec = (!main.df$SDIU.body.sig) & (is.na(main.df$SDIU.head.sig))
SDIUinHead.NotTestedBody.vec = (main.df$SDIU.head.sig) & (is.na(main.df$SDIU.body.sig))
Not.SDIUinHead.NotTestedBody.vec = (!main.df$SDIU.head.sig) & (is.na(main.df$SDIU.body.sig))

vector.sums =c(sum(SDIUinBoth.vec, na.rm = T),
               sum(Not.SDIUinBoth.vec, na.rm = T),
               sum(Not.Tested.Both.vec, na.rm = T),
               sum(SDIUinBody.Not.SDIUinHead.vec, na.rm = T),
               sum(SDIUinBody.NotTestedHead.vec, na.rm = T),
               sum(Not.SDIUinBody.SDIUinHead.vec, na.rm = T),
               sum(Not.SDIUinBody.NotTestedHead.vec, na.rm = T),
               sum(SDIUinHead.NotTestedBody.vec, na.rm = T),
               sum(Not.SDIUinHead.NotTestedBody.vec, na.rm = T))

values.For.Table1 = rbind(vector.sums[c(1,4,5)], vector.sums[c(6,2, 7)], vector.sums[c(8,9,3)])
colnames(values.For.Table1)<-c("Head.SDIU", "Head.NotSDIU", "Head.NA")
rownames(values.For.Table1)<-c("Body.SDIU", "Body.NotSDIU", "Body.NA")
values.For.Table1

### summary of body SSS
rowSums(values.For.Table1)
rowSums(values.For.Table1)[1]/sum(rowSums(values.For.Table1)[1:2])

### summary of head SSS
colSums(values.For.Table1)
colSums(values.For.Table1)[1]/sum(colSums(values.For.Table1)[1:2])

### summary of SSS in body OR head
(rowSums(values.For.Table1)[1] + colSums(values.For.Table1)[1] - values.For.Table1[1,1])
(rowSums(values.For.Table1)[1] + colSums(values.For.Table1)[1] - values.For.Table1[1,1])/(sum(values.For.Table1) - values.For.Table1[3,3])



############  Association of SDIU and SBGE
###### BODY: Proportion SDIU wrt to SBGE
## here I make a plot that shows the proportion of SDIU in body
## Black points are just wrt to SBGE category
## The 3 adjacement points are low, mid, and high baseMean expression
## The expression categories were chosen simply by lowest 1/3, mid 1/3, and top 1/3 for this tissue (i.e., body).


MakeSimpleFig.BodyAndHead<-function(df, xvec = 1:8, this.ylim){
  xvals.for.fig = xvec
  xadj = c(round(xvec[1:(length(xvec)/2)]), nrow(df)/2 + round(xvec[1:(length(xvec)/2)]))
  minSampleSize = 20
  vecToExcludeLowSampleSizePoints = sapply(df[,2], function(x) ifelse(x > minSampleSize, 1, NA))
  xadj = vecToExcludeLowSampleSizePoints*xadj
  vec.col.for.fig = c(rep("orange", length(xvec)/2), rep("blue", length(xvec)/2))
  plot(xvals.for.fig, df[xadj,3], 
       col = vec.col.for.fig,
       ylim = this.ylim, #c(min(df[xvec,4]), max(df[xvec,5])), 
       #xlab = "SBGE category", ylab = "Proportion SDIU", 
       pch = 16, xaxt = 'n')
  arrows(x0=xvals.for.fig, y0= df[xadj,4],
         x1 = xvals.for.fig, y1= df[xadj,5], 
         code=3, angle = 90, length = 0.05, 
         vec.col.for.fig, lwd=2)
  abline(v = 1.5, lwd = 2)
  abline(v = 0.5 + 2:7, lwd = 0.5)
}

summaryPropSDIUbySBGE.body.Osada = MakeSummaryOfVarByFactor("SDIU.body.sig", "SBGEcat.body.Osada", main.df)
summaryPropSDIUbySBGE.head.Osada = MakeSummaryOfVarByFactor("SDIU.head.sig", "SBGEcat.head.Osada", main.df)

summaryPropXbySDIU.body = MakeSummaryOfVarByFactor("Is.X", "SDIU.body.sig", main.df)
summaryPropXbySDIU.head = MakeSummaryOfVarByFactor("Is.X", "SDIU.head.sig", main.df)


MakeSimpleFig.BodyAndHead(rbind(summaryPropSDIUbySBGE.body.Osada, summaryPropSDIUbySBGE.head.Osada),
                    xvec= c(-0.1 + (1:8), 0.1 + (1:8)), this.ylim = c(0, 0.8))

dt.SDIUwrtSBGE.body.Osada = main.df[,c("geneID","SDIU.body.sig", "Whole.SBGE.Osada", "SBGEcat.body.Osada")]
dt.SDIUwrtSBGE.body.Osada = na.omit(dt.SDIUwrtSBGE.body.Osada)
dt.SDIUwrtSBGE.body.Osada.ExlExtSB = dt.SDIUwrtSBGE.body.Osada[dt.SDIUwrtSBGE.body.Osada$SBGEcat.body.Osada %in% c("sFB", "FB", "UB", "MB", "sMB"),]

dt.SDIUwrtSBGE.head.Osada = main.df[,c("geneID", "SDIU.head.sig","Head.SBGE.Osada", "SBGEcat.head.Osada")]
dt.SDIUwrtSBGE.head.Osada = na.omit(dt.SDIUwrtSBGE.head.Osada)
dt.SDIUwrtSBGE.head.Osada.ExlExtSB = dt.SDIUwrtSBGE.head.Osada[dt.SDIUwrtSBGE.head.Osada$SBGEcat.head.Osada %in% c("sFB", "FB", "UB", "MB", "sMB"),]


ProportionSDIUwrtSBGE.mod.Body.Only.Osada.ExclExtSB = glm(SDIU.body.sig ~  Whole.SBGE.Osada + I(Whole.SBGE.Osada^2),
                                                                                data = dt.SDIUwrtSBGE.body.Osada.ExlExtSB, family = binomial)
summary(ProportionSDIUwrtSBGE.mod.Body.Only.Osada.ExclExtSB)

boot_case.ProportionSDIUwrtSBGE.mod.Body.Only.Osada.ExclExtSB <- Boot(ProportionSDIUwrtSBGE.mod.Body.Only.Osada.ExclExtSB, method = "case", R = 9999)
confint(boot_case.ProportionSDIUwrtSBGE.mod.Body.Only.Osada.ExclExtSB, type = "perc")

ProportionSDIUwrtSBGE.mod.Head.Only.Osada.ExclExtSB = glm(SDIU.head.sig ~  Head.SBGE.Osada + I(Head.SBGE.Osada^2),
                                                                                data = dt.SDIUwrtSBGE.head.Osada.ExlExtSB, family = binomial)
summary(ProportionSDIUwrtSBGE.mod.Head.Only.Osada.ExclExtSB)

boot_case.ProportionSDIUwrtSBGE.mod.Head.Only.Osada.ExclExtSB <- Boot(ProportionSDIUwrtSBGE.mod.Head.Only.Osada.ExclExtSB, method = "case", R = 9999)
confint(boot_case.ProportionSDIUwrtSBGE.mod.Head.Only.Osada.ExclExtSB, type = "perc")


#############  examining only genes that are tested for SDIU and SBGE in both tissues
dt.SDIUwrtSBGE.Osada = main.df[,c("geneID","SDIU.body.sig", "SDIU.head.sig",
           "Whole.SBGE.Osada", "Head.SBGE.Osada", "SBGEcat.body.Osada", "SBGEcat.head.Osada")]
dt.SDIUwrtSBGE.Osada.GenesInBothTissues = na.omit(dt.SDIUwrtSBGE.Osada)
dt.SDIUwrtSBGE.Osada.GenesInBothTissues.ExlExtSB = dt.SDIUwrtSBGE.Osada.GenesInBothTissues[dt.SDIUwrtSBGE.Osada.GenesInBothTissues$SBGEcat.body.Osada %in% c("sFB", "FB", "UB", "MB", "sMB") & dt.SDIUwrtSBGE.Osada.GenesInBothTissues$SBGEcat.head.Osada %in% c("sFB", "FB", "UB", "MB", "sMB"),]


summaryPropSDIUbySBGE.body.Osada.GenesInBothTissues = MakeSummaryOfVarByFactor("SDIU.body.sig", "SBGEcat.body.Osada", dt.SDIUwrtSBGE.Osada.GenesInBothTissues)
summaryPropSDIUbySBGE.head.Osada.GenesInBothTissues = MakeSummaryOfVarByFactor("SDIU.head.sig", "SBGEcat.head.Osada", dt.SDIUwrtSBGE.Osada.GenesInBothTissues)


MakeSimpleFig.BodyAndHead(rbind(summaryPropSDIUbySBGE.body.Osada.GenesInBothTissues, summaryPropSDIUbySBGE.head.Osada.GenesInBothTissues),
                    xvec= c(-0.1 + (1:8), 0.1 + (1:8)), this.ylim = c(0, 0.8))

ProportionSDIUwrtSBGE.mod.Body.OnlyGenesTestedBothTissues.Osada.ExclExtSB = glm(SDIU.body.sig ~  Whole.SBGE.Osada + I(Whole.SBGE.Osada^2),
                                                                data = dt.SDIUwrtSBGE.Osada.GenesInBothTissues.ExlExtSB, family = binomial)
summary(ProportionSDIUwrtSBGE.mod.Body.OnlyGenesTestedBothTissues.Osada.ExclExtSB)

boot_case.ProportionSDIUwrtSBGE.mod.Body.Osada.ExclExtSB <- Boot(ProportionSDIUwrtSBGE.mod.Body.OnlyGenesTestedBothTissues.Osada.ExclExtSB, method = "case", R = 9999)
confint(boot_case.ProportionSDIUwrtSBGE.mod.Body.Osada.ExclExtSB, type = "perc")

ProportionSDIUwrtSBGE.mod.Head.OnlyGenesTestedBothTissues.Osada.ExclExtSB = glm(SDIU.head.sig ~  Head.SBGE.Osada + I(Head.SBGE.Osada^2),
                                                                                data = dt.SDIUwrtSBGE.Osada.GenesInBothTissues.ExlExtSB, family = binomial)
summary(ProportionSDIUwrtSBGE.mod.Head.OnlyGenesTestedBothTissues.Osada.ExclExtSB)

boot_case.ProportionSDIUwrtSBGE.mod.Head.Osada.ExclExtSB <- Boot(ProportionSDIUwrtSBGE.mod.Head.OnlyGenesTestedBothTissues.Osada.ExclExtSB, method = "case", R = 9999)
confint(boot_case.ProportionSDIUwrtSBGE.mod.Head.Osada.ExclExtSB, type = "perc")

###############



################## Proportion X-linked

summaryPropXbySDIU.body = MakeSummaryOfVarByFactor("Is.X", "SDIU.body.sig", main.df)
summaryPropXbySDIU.head = MakeSummaryOfVarByFactor("Is.X", "SDIU.head.sig", main.df)


#### Body - proportion X
## Use permutation test to  test if proportion X
## is different between SDIU and Not.SDIU genes
SDIU.Is.X.body.vec = main.df$Is.X[main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$Is.X)]
Not.SDIU.Is.X.body.vec = main.df$Is.X[!main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$Is.X)]

Diff.mean.Is.X.SDIUvsNotSDIU.body = mean(SDIU.Is.X.body.vec) - mean(Not.SDIU.Is.X.body.vec)
Diff.mean.Is.X.SDIUvsNotSDIU.body

SDIU.NotSDIU.Is.X.body.vec = c(SDIU.Is.X.body.vec, Not.SDIU.Is.X.body.vec)

Diff.mean.Is.X.SDIUvsNotSDIU.body.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.Is.X.body.vec, length(SDIU.Is.X.body.vec)))
sum(Diff.mean.Is.X.SDIUvsNotSDIU.body.null.dist <= Diff.mean.Is.X.SDIUvsNotSDIU.body)/length(Diff.mean.Is.X.SDIUvsNotSDIU.body.null.dist)
sum(Diff.mean.Is.X.SDIUvsNotSDIU.body.null.dist >= Diff.mean.Is.X.SDIUvsNotSDIU.body)/length(Diff.mean.Is.X.SDIUvsNotSDIU.body.null.dist)


## FB genes are enriched on X
SBGEandSDIU.Is.X.body.for.simple.fig = main.df[ ,c("geneID","Is.X","SDIU.body.sig","SDIU.body.p.strength.cat", "Whole.SBGE.Osada", "SBGEcat.body.Osada", "SDIUxSBGE.body.Osada")]

summary(main.df$SDIUxSBGE.body.Osada)

## These plots use 95% bootstrap CIs for error bars.
##  (black - all; grey- NotSDIU; purple - SDIU)
sorted.combined.summary.PropX.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("Is.X","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.Is.X.body.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.PropX.by.SDIUxSBGE, "Proportion X-linked", this.ylim = c(0, 0.45))


#### HEAD - proportion X-lined
## Use permutation test to test if proportion X
## is different between SDIU and Not.SDIU genes
SDIU.Is.X.head.vec = main.df$Is.X[main.df$SDIU.head.sig & !is.na(main.df$SDIU.head.sig) & !is.na(main.df$Is.X)]
Not.SDIU.Is.X.head.vec = main.df$Is.X[!main.df$SDIU.head.sig & !is.na(main.df$SDIU.head.sig) & !is.na(main.df$Is.X)]

Diff.mean.Is.X.SDIUvsNotSDIU.head = mean(SDIU.Is.X.head.vec) - mean(Not.SDIU.Is.X.head.vec)
Diff.mean.Is.X.SDIUvsNotSDIU.head

SDIU.NotSDIU.Is.X.head.vec = c(SDIU.Is.X.head.vec, Not.SDIU.Is.X.head.vec)

Diff.mean.Is.X.SDIUvsNotSDIU.head.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.Is.X.head.vec, length(SDIU.Is.X.head.vec)))
sum(Diff.mean.Is.X.SDIUvsNotSDIU.head.null.dist <= Diff.mean.Is.X.SDIUvsNotSDIU.head)/length(Diff.mean.Is.X.SDIUvsNotSDIU.head.null.dist)
sum(Diff.mean.Is.X.SDIUvsNotSDIU.head.null.dist >= Diff.mean.Is.X.SDIUvsNotSDIU.head)/length(Diff.mean.Is.X.SDIUvsNotSDIU.head.null.dist)


## summarize with respect to SDIU and SBGE
## These plots use 95% bootstrap CIs for error bars.
##  (black - all; grey- NotSDIU; purple - SDIU)

SBGEandSDIU.Is.X.head.for.simple.fig = main.df[ ,c("geneID","Is.X","SDIU.head.sig","SDIU.head.p.strength.cat", "Head.SBGE.Osada", "SBGEcat.head.Osada", "SDIUxSBGE.head.Osada")]


sorted.combined.summary.PropX.by.SDIUxSBGE.Head = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("Is.X", "SBGEcat.head.Osada", "SDIUxSBGE.head.Osada", SBGEandSDIU.Is.X.head.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Head(sorted.combined.summary.PropX.by.SDIUxSBGE.Head, "Proportion X-linked", this.ylim = c(0, 0.85))

sorted.combined.summary.PropX.by.SDIUxSBGE.Head.adjForFig = sorted.combined.summary.PropX.by.SDIUxSBGE.Head
sorted.combined.summary.PropX.by.SDIUxSBGE.Head.adjForFig[7, c(2:6)]=NA
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.PropX.by.SDIUxSBGE.Head.adjForFig, "Proportion X-linked", this.ylim = c(0, 0.85))


### body glm
dt.ProportionX.wrt.Body.SBGE = na.omit(main.df[,c("geneID","Is.X", "Whole.SBGE.Osada", "SBGEcat.body.Osada")])
dt.ProportionX.wrt.Body.SBGE.ExlExtSB = dt.ProportionX.wrt.Body.SBGE[dt.ProportionX.wrt.Body.SBGE$SBGEcat.body.Osada %in% c("sFB", "FB", "UB", "MB", "sMB") ,]

ProportionX.wrt.Body.SBGE.ExlExtSB = glm(Is.X ~  Whole.SBGE.Osada + I(Whole.SBGE.Osada^2),
                                                                              data = dt.ProportionX.wrt.Body.SBGE.ExlExtSB, family = binomial)
summary(ProportionX.wrt.Body.SBGE.ExlExtSB)

boot_case.ProportionX.wrt.Body.SBGE.ExlExtSB <- Boot(ProportionX.wrt.Body.SBGE.ExlExtSB, method = "case", R = 9999)
confint(boot_case.ProportionX.wrt.Body.SBGE.ExlExtSB, type = "perc")

### head glm
dt.ProportionX.wrt.Head.SBGE = na.omit(main.df[,c("geneID","Is.X", "Head.SBGE.Osada", "SBGEcat.head.Osada")])
dt.ProportionX.wrt.Head.SBGE.ExlExtSB = dt.ProportionX.wrt.Head.SBGE[dt.ProportionX.wrt.Head.SBGE$SBGEcat.head.Osada %in% c("sFB", "FB", "UB", "MB", "sMB") ,]

ProportionX.wrt.Head.SBGE.ExlExtSB = glm(Is.X ~  Head.SBGE.Osada + I(Head.SBGE.Osada^2),
                                         data = dt.ProportionX.wrt.Head.SBGE.ExlExtSB, family = binomial)
summary(ProportionX.wrt.Head.SBGE.ExlExtSB)

boot_case.ProportionX.wrt.Head.SBGE.ExlExtSB <- Boot(ProportionX.wrt.Head.SBGE.ExlExtSB, method = "case", R = 9999)
confint(boot_case.ProportionX.wrt.Head.SBGE.ExlExtSB, type = "perc")


####################  average expression: log.RPKM.avg
### BODY - log.avgExp.AdultLarva (calculated from Fly Atlas 2 [FA2])

## Use permutation test to test if average expression
## is different for non-SDIU than SDIU genes
SDIU.avgExp.AdultLarva.body.vec = main.df$log.avgExp.AdultLarva[main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$log.avgExp.AdultLarva)]
Not.SDIU.avgExp.AdultLarva.body.vec = main.df$log.avgExp.AdultLarva[!main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$log.avgExp.AdultLarva)]

Diff.mean.log.avgExp.AdultLarva.SDIUvsNotSDIU.body = mean(SDIU.avgExp.AdultLarva.body.vec) - mean(Not.SDIU.avgExp.AdultLarva.body.vec)
Diff.mean.log.avgExp.AdultLarva.SDIUvsNotSDIU.body

SDIU.NotSDIU.avgExp.AdultLarva.body.vec = c(SDIU.avgExp.AdultLarva.body.vec, Not.SDIU.avgExp.AdultLarva.body.vec)


Diff.mean.log.avgExp.AdultLarva.body.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.avgExp.AdultLarva.body.vec, length(SDIU.Is.X.body.vec)))
sum(Diff.mean.log.avgExp.AdultLarva.body.null.dist <= Diff.mean.log.avgExp.AdultLarva.SDIUvsNotSDIU.body)/length(Diff.mean.log.avgExp.AdultLarva.body.null.dist)
sum(Diff.mean.log.avgExp.AdultLarva.body.null.dist >= Diff.mean.log.avgExp.AdultLarva.SDIUvsNotSDIU.body)/length(Diff.mean.log.avgExp.AdultLarva.body.null.dist)

#head
SDIU.avgExp.AdultLarva.head.vec = main.df$log.avgExp.AdultLarva[main.df$SDIU.head.sig & !is.na(main.df$SDIU.head.sig) & !is.na(main.df$log.avgExp.AdultLarva)]
Not.SDIU.avgExp.AdultLarva.head.vec = main.df$log.avgExp.AdultLarva[!main.df$SDIU.head.sig & !is.na(main.df$SDIU.head.sig) & !is.na(main.df$log.avgExp.AdultLarva)]

Diff.mean.log.avgExp.AdultLarva.SDIUvsNotSDIU.head = mean(SDIU.avgExp.AdultLarva.head.vec) - mean(Not.SDIU.avgExp.AdultLarva.head.vec)
Diff.mean.log.avgExp.AdultLarva.SDIUvsNotSDIU.head

SDIU.NotSDIU.avgExp.AdultLarva.body.vec = c(SDIU.avgExp.AdultLarva.head.vec, Not.SDIU.avgExp.AdultLarva.head.vec)


Diff.mean.log.avgExp.AdultLarva.head.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.avgExp.AdultLarva.body.vec, length(SDIU.Is.X.head.vec)))
sum(Diff.mean.log.avgExp.AdultLarva.head.null.dist <= Diff.mean.log.avgExp.AdultLarva.SDIUvsNotSDIU.head)/length(Diff.mean.log.avgExp.AdultLarva.head.null.dist)
sum(Diff.mean.log.avgExp.AdultLarva.head.null.dist >= Diff.mean.log.avgExp.AdultLarva.SDIUvsNotSDIU.head)/length(Diff.mean.log.avgExp.AdultLarva.head.null.dist)


### Body - log.avgExp.AdultLarva from FA2 
SBGEandSDIU.logRPKM.FA2.body.for.simple.fig = main.df[ ,c("geneID","log.avgExp.AdultLarva","SDIU.body.sig","SDIU.body.p.strength.cat", "Whole.SBGE.Osada", "SBGEcat.body.Osada", "SDIUxSBGE.body.Osada")]
sorted.combined.summary.logRPKM.FA2.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("log.avgExp.AdultLarva","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.logRPKM.FA2.body.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.logRPKM.FA2.by.SDIUxSBGE, "log(avgExp.AdultsLarva)", this.ylim = c(0, 3))

### HEAD - log.avgExp.AdultLarva from FA2 
SBGEandSDIU.logRPKM.FA2.for.simple.fig.head = main.df[ ,c("geneID","log.avgExp.AdultLarva","SDIU.head.sig","SDIU.head.p.strength.cat", "Head.SBGE.Osada", "SBGEcat.head.Osada", "SDIUxSBGE.head.Osada")]
sorted.combined.summary.logRPKM.FA2.by.SDIUxSBGE.head = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("log.avgExp.AdultLarva","SBGEcat.head.Osada", "SDIUxSBGE.head.Osada", SBGEandSDIU.logRPKM.FA2.for.simple.fig.head)

sorted.combined.summary.logRPKM.FA2.by.SDIUxSBGE.head.AdjForFig = sorted.combined.summary.logRPKM.FA2.by.SDIUxSBGE.head;
sorted.combined.summary.logRPKM.FA2.by.SDIUxSBGE.head.AdjForFig[7, c(2:6)] = NA
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.logRPKM.FA2.by.SDIUxSBGE.head.AdjForFig, "log(avgExp.AdultsLarva)", this.ylim = c(0, 4.5))


### body -  linear model  for expression
dt.log.avgExp.AdultLarva.wrt.Body.SBGE = na.omit(main.df[,c("geneID","log.avgExp.AdultLarva", "Whole.SBGE.Osada", "SBGEcat.body.Osada")])
dt.log.avgExp.AdultLarva.wrt.Body.SBGE.ExlExtSB = dt.log.avgExp.AdultLarva.wrt.Body.SBGE[dt.log.avgExp.AdultLarva.wrt.Body.SBGE$SBGEcat.body.Osada %in% c("sFB", "FB", "UB", "MB", "sMB") ,]

log.avgExp.AdultLarva.wrt.Body.SBGE.ExlExtSB = lm(log.avgExp.AdultLarva ~  Whole.SBGE.Osada + I(Whole.SBGE.Osada^2),
                                         data = dt.log.avgExp.AdultLarva.wrt.Body.SBGE.ExlExtSB)
summary(log.avgExp.AdultLarva.wrt.Body.SBGE.ExlExtSB)

boot_case.log.avgExp.AdultLarva.wrt.Body.SBGE.ExlExtSB <- Boot(log.avgExp.AdultLarva.wrt.Body.SBGE.ExlExtSB, method = "case", R = 9999)
confint(boot_case.log.avgExp.AdultLarva.wrt.Body.SBGE.ExlExtSB, type = "perc")

### head -  linear model for expression
dt.log.avgExp.AdultLarva.wrt.Head.SBGE = na.omit(main.df[,c("geneID","log.avgExp.AdultLarva", "Head.SBGE.Osada", "SBGEcat.body.Osada")])
dt.log.avgExp.AdultLarva.wrt.Head.SBGE.ExlExtSB = dt.log.avgExp.AdultLarva.wrt.Head.SBGE[dt.log.avgExp.AdultLarva.wrt.Head.SBGE$SBGEcat.body.Osada %in% c("sFB", "FB", "UB", "MB", "sMB") ,]

log.avgExp.AdultLarva.wrt.Head.SBGE.ExlExtSB = lm(log.avgExp.AdultLarva ~  Head.SBGE.Osada + I(Head.SBGE.Osada^2),
                                                  data = dt.log.avgExp.AdultLarva.wrt.Head.SBGE.ExlExtSB)
summary(log.avgExp.AdultLarva.wrt.Head.SBGE.ExlExtSB)

boot_case.log.avgExp.AdultLarva.wrt.Head.SBGE.ExlExtSB <- Boot(log.avgExp.AdultLarva.wrt.Head.SBGE.ExlExtSB, method = "case", R = 9999)
confint(boot_case.log.avgExp.AdultLarva.wrt.Head.SBGE.ExlExtSB, type = "perc")





################## expression pleiotropy: tissue specificity: tao.FML.FA2
###BODY
## Use permutation test to test if mean tao.FML.FA2 
## is  significantly lower for SDIU genes than non-SDIU genes
SDIU.taoFML.vec = na.omit(main.df$tao.FML.FA2[main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$tao.FML.FA2)])
Not.SDIU.taoFML.vec = na.omit(main.df$tao.FML.FA2[!main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$tao.FML.FA2)])

Diff.mean.taoFML.SDIUvsNotSDIU = mean(SDIU.taoFML.vec) - mean(Not.SDIU.taoFML.vec)
Diff.mean.taoFML.SDIUvsNotSDIU

SDIU.NotSDIU.taoFML.vec = c(SDIU.taoFML.vec, Not.SDIU.taoFML.vec)
Diff.mean.taoFML.SDIUvsNotSDIU.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.taoFML.vec, length(SDIU.taoFML.vec)))
sum(Diff.mean.taoFML.SDIUvsNotSDIU.null.dist <= Diff.mean.taoFML.SDIUvsNotSDIU)/length(Diff.mean.taoFML.SDIUvsNotSDIU.null.dist)
sum(Diff.mean.taoFML.SDIUvsNotSDIU.null.dist >= Diff.mean.taoFML.SDIUvsNotSDIU)/length(Diff.mean.taoFML.SDIUvsNotSDIU.null.dist)

SBGEandSDIU.taoFML.body.for.simple.fig = main.df[ ,c("geneID","tao.FML.FA2","SDIU.body.sig","SDIU.body.p.strength.cat", "Whole.SBGE.Osada", "SBGEcat.body.Osada", "SDIUxSBGE.body.Osada")]
sorted.combined.summary.taoFML.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("tao.FML.FA2","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.taoFML.body.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.taoFML.by.SDIUxSBGE, "tao.FML.FA2", this.ylim = c(0, 1))

### HEAD
SDIU.taoFML.head.vec = na.omit(main.df$tao.FML.FA2[main.df$SDIU.head.sig & !is.na(main.df$SDIU.head.sig) & !is.na(main.df$tao.FML.FA2)])
Not.SDIU.taoFML.head.vec = na.omit(main.df$tao.FML.FA2[!main.df$SDIU.head.sig & !is.na(main.df$SDIU.head.sig) & !is.na(main.df$tao.FML.FA2)])

Diff.mean.taoFML.SDIUvsNotSDIU.head = mean(SDIU.taoFML.head.vec) - mean(Not.SDIU.taoFML.head.vec)
Diff.mean.taoFML.SDIUvsNotSDIU.head

SDIU.NotSDIU.taoFML.head.vec = c(SDIU.taoFML.head.vec, Not.SDIU.taoFML.head.vec)
Diff.mean.taoFML.SDIUvsNotSDIU.head.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.taoFML.head.vec, length(SDIU.taoFML.head.vec)))
sum(Diff.mean.taoFML.SDIUvsNotSDIU.head.null.dist <= Diff.mean.taoFML.SDIUvsNotSDIU.head)/length(Diff.mean.taoFML.SDIUvsNotSDIU.head.null.dist)
sum(Diff.mean.taoFML.SDIUvsNotSDIU.head.null.dist >= Diff.mean.taoFML.SDIUvsNotSDIU.head)/length(Diff.mean.taoFML.SDIUvsNotSDIU.head.null.dist)



SBGEandSDIU.taoFML.head.for.simple.fig = main.df[ ,c("geneID","tao.FML.FA2","SDIU.head.sig","SDIU.head.p.strength.cat","baseMean.head.Osada", "Head.SBGE.Osada", "SBGEcat.head.Osada", "SDIUxSBGE.head.Osada")]
sorted.combined.summary.taoFML.by.SDIUxSBGE.head = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("tao.FML.FA2","SBGEcat.head.Osada", "SDIUxSBGE.head.Osada", SBGEandSDIU.taoFML.head.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.taoFML.by.SDIUxSBGE.head, "tao.FML.FA2", this.ylim = c(0, 1))



### body - linear model
dt.tao.FML.FA2.wrt.Body.SBGE = na.omit(main.df[,c("geneID","tao.FML.FA2", "Whole.SBGE.Osada", "SBGEcat.body.Osada")])
dt.tao.FML.FA2.wrt.Body.SBGE.ExlExtSB = dt.tao.FML.FA2.wrt.Body.SBGE[dt.tao.FML.FA2.wrt.Body.SBGE$SBGEcat.body.Osada %in% c("sFB", "FB", "UB", "MB", "sMB") ,]

tao.FML.FA2.wrt.Body.SBGE.ExlExtSB = lm(tao.FML.FA2 ~  Whole.SBGE.Osada + I(Whole.SBGE.Osada^2),
                                                  data = dt.tao.FML.FA2.wrt.Body.SBGE.ExlExtSB)
summary(tao.FML.FA2.wrt.Body.SBGE.ExlExtSB)

boot_case.tao.FML.FA2.wrt.Body.SBGE.ExlExtSB <- Boot(tao.FML.FA2.wrt.Body.SBGE.ExlExtSB, method = "case", R = 9999)
confint(boot_case.tao.FML.FA2.wrt.Body.SBGE.ExlExtSB, type = "perc")

### head - lm tissue specificity wrt SBGE
dt.tao.FML.FA2.wrt.Head.SBGE = na.omit(main.df[,c("geneID","tao.FML.FA2", "Head.SBGE.Osada", "SBGEcat.body.Osada")])
dt.tao.FML.FA2.wrt.Head.SBGE.ExlExtSB = dt.tao.FML.FA2.wrt.Head.SBGE[dt.tao.FML.FA2.wrt.Head.SBGE$SBGEcat.body.Osada %in% c("sFB", "FB", "UB", "MB", "sMB") ,]

tao.FML.FA2.wrt.Head.SBGE.ExlExtSB = lm(tao.FML.FA2 ~  Head.SBGE.Osada + I(Head.SBGE.Osada^2),
                                                  data = dt.tao.FML.FA2.wrt.Head.SBGE.ExlExtSB)
summary(tao.FML.FA2.wrt.Head.SBGE.ExlExtSB)

boot_case.tao.FML.FA2.wrt.Head.SBGE.ExlExtSB <- Boot(tao.FML.FA2.wrt.Head.SBGE.ExlExtSB, method = "case", R = 9999)
confint(boot_case.tao.FML.FA2.wrt.Head.SBGE.ExlExtSB, type = "perc")


###################################################
ExclExtSB.Body.Osada.vec = main.df$SBGEcat.body.Osada %in% c("sFB", "FB", "UB", "MB", "sMB")
summary(main.df$SBGEcat.body.Osada)
summary(main.df$SBGEcat.body.Osada[ExclExtSB.Body.Osada.vec])

#######  rmf patterns

## note we will do some filtering by Vg 
## first let's just inspect it
summary(main.df$Vg.Sig0.05.Either)
summary(main.df$Vg.Sig0.05.Both)
summary(main.df$SBGEcat.body.Osada)
summary(main.df[main.df$Vg.Sig0.05.Either & !is.na(main.df$Vg.Sig0.05.Either) ,"SBGEcat.body.Osada"])
summary(main.df[main.df$Vg.Sig0.05.Both & !is.na(main.df$Vg.Sig0.05.Both) ,"SBGEcat.body.Osada"])


## Use permutation test to test if mean rmf  is significantly lower for SDIU genes than non-SDIU genes
SDIU.rmf.vec = na.omit(main.df$rmf[main.df$Vg.Sig0.05.Either & main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig)])
Not.SDIU.rmf.vec = na.omit(main.df$rmf[main.df$Vg.Sig0.05.Either & !main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig)])

Diff.mean.rmf.SDIUvsNotSDIU = mean(SDIU.rmf.vec) - mean(Not.SDIU.rmf.vec)
Diff.mean.rmf.SDIUvsNotSDIU

SDIU.NotSDIU.rmf.vec = c(SDIU.rmf.vec, Not.SDIU.rmf.vec)

Diff.mean.rmf.SDIUvsNotSDIU.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.rmf.vec, length(SDIU.rmf.vec)))
sum(Diff.mean.rmf.SDIUvsNotSDIU <= Diff.mean.rmf.SDIUvsNotSDIU.null.dist)/length(Diff.mean.rmf.SDIUvsNotSDIU.null.dist)
sum(Diff.mean.rmf.SDIUvsNotSDIU >= Diff.mean.rmf.SDIUvsNotSDIU.null.dist)/length(Diff.mean.rmf.SDIUvsNotSDIU.null.dist)


## Below is a table that shows this difference between SDIU and non-SDIU
### persists across (most SBGE categories)
## In general, rmf is strongest for moderately MB genes 

SBGEandSDIU.rmf.body.for.simple.fig = main.df[main.df$Vg.Sig0.05.Either  & !is.na(main.df$Vg.Sig0.05.Either) ,c("geneID","rmf","log.avgExp.AdultLarva","SDIU.body.sig","SDIU.body.p.strength.cat", "Whole.SBGE.Osada", "SBGEcat.body.Osada", "SDIUxSBGE.body.Osada")]
sorted.combined.summary.rmf.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("rmf","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.rmf.body.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.rmf.by.SDIUxSBGE, "rmf", this.ylim = c(0, 1))


### version requiring Vg Sig for both sexes
SBGEandSDIU.rmf.body.for.simple.fig.VgSigBoth = main.df[main.df$Vg.Sig0.05.Both  & !is.na(main.df$Vg.Sig0.05.Both) ,c("geneID","rmf","log.avgExp.AdultLarva","SDIU.body.sig","SDIU.body.p.strength.cat", "Whole.SBGE.Osada", "SBGEcat.body.Osada", "SDIUxSBGE.body.Osada")]
sorted.combined.summary.rmf.by.SDIUxSBGE.VgSigBoth = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("rmf","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.rmf.body.for.simple.fig.VgSigBoth)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.rmf.by.SDIUxSBGE.VgSigBoth, "rmf", this.ylim = c(0, 1))




####### rmf linear model (excludes extreme SB)

SBGEandSDIU.Body.rmf.mod.data.complete.ExclExtSB.withPCs = main.df[ExclExtSB.Body.Osada.vec & main.df$Vg.Sig0.05.Either ,c("geneID","rmf","SDIU.body.p.strength.cat","log.avgExp.AdultLarva", "Whole.SBGE.Osada", 
                                                                                                                      "Is.X",  "tao.FML.FA2",
                                                                                                                      "m.Comp.1", "m.Comp.2", "m.Comp.3",
                                                                                                                      "f.Comp.1", "f.Comp.2", "f.Comp.3" )]
SBGEandSDIU.Body.rmf.mod.data.complete.ExclExtSB.withPCs = na.omit(SBGEandSDIU.Body.rmf.mod.data.complete.ExclExtSB.withPCs)

options(contrasts = rep("contr.sum", 2)) ## this should make lmp and lm results comparable
mod.DPGP3.SDIU.rmf.withPCs = lm(rmf ~ SDIU.body.p.strength.cat 
                                   + log.avgExp.AdultLarva 
                                   + Whole.SBGE.Osada 
                                   + I(Whole.SBGE.Osada^2) 
                                   + Is.X + tao.FML.FA2 
                                   + m.Comp.1 + m.Comp.2 + m.Comp.3 
                                   + f.Comp.1 + f.Comp.2 + f.Comp.3, 
                                   data = SBGEandSDIU.Body.rmf.mod.data.complete.ExclExtSB.withPCs)
summary(mod.DPGP3.SDIU.rmf.withPCs)

mod.lmp.rmf.ExclExtSB.withPCs = lmp(rmf ~ SDIU.body.p.strength.cat 
                                    + log.avgExp.AdultLarva 
                                    + Whole.SBGE.Osada 
                                    + I(Whole.SBGE.Osada^2) 
                                    + Is.X + tao.FML.FA2  
                                    + m.Comp.1 + m.Comp.2 + m.Comp.3 
                                    + f.Comp.1 + f.Comp.2 + f.Comp.3,
                            center = FALSE,
                            perm = "Prob",
                            Ca = 1e-3,
                            maxIter = 1e6,
                            data = SBGEandSDIU.Body.rmf.mod.data.complete.ExclExtSB.withPCs)
summary(mod.lmp.rmf.ExclExtSB.withPCs)


boot_case.mod.rmf.ExclExtSB.withPCs <- Boot(mod.DPGP3.SDIU.rmf.withPCs, method = "case", R = 9999)
confint(boot_case.mod.rmf.ExclExtSB.withPCs, type = "perc")





##################################  end: rmf   ##############################################################





##################### begin:  tajD.S  ##############################

## Use permutation test to test if mean tajD.S  
## is  significantly different for SDIU genes than non-SDIU genes
Above4FoldMin.vec = main.df$n4fold.usedForPI >= 30; summary(Above4FoldMin.vec)
Above4FoldMin.vec = Above4FoldMin.vec & !is.na(Above4FoldMin.vec); summary(Above4FoldMin.vec)

SDIU.tajD.S.vec = na.omit(main.df$tajD.S[Above4FoldMin.vec & main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$tajD.S)])
Not.SDIU.tajD.S.vec = na.omit(main.df$tajD.S[Above4FoldMin.vec & !main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$tajD.S)])

Diff.mean.tajD.S.SDIUvsNotSDIU = mean(SDIU.tajD.S.vec) - mean(Not.SDIU.tajD.S.vec)
Diff.mean.tajD.S.SDIUvsNotSDIU

SDIU.NotSDIU.tajD.S.vec = c(SDIU.tajD.S.vec, Not.SDIU.tajD.S.vec)

Diff.mean.tajD.S.SDIUvsNotSDIU.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.tajD.S.vec, length(SDIU.tajD.S.vec)))
sum(Diff.mean.tajD.S.SDIUvsNotSDIU.null.dist <= Diff.mean.tajD.S.SDIUvsNotSDIU)/length(Diff.mean.tajD.S.SDIUvsNotSDIU.null.dist)
sum(Diff.mean.tajD.S.SDIUvsNotSDIU.null.dist >= Diff.mean.tajD.S.SDIUvsNotSDIU)/length(Diff.mean.tajD.S.SDIUvsNotSDIU.null.dist)



SBGEandSDIU.Body.tajD.S.for.simple.fig = main.df[Above4FoldMin.vec, ]

sorted.combined.summary.tajD.S.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("tajD.S","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.Body.tajD.S.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.tajD.S.by.SDIUxSBGE, "tajD.S", this.ylim = c(-0.14,-0.08))

## head fig
sorted.combined.summary.tajD.S.by.SDIUxSBGE.head = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("tajD.S","SBGEcat.head.Osada", "SDIUxSBGE.head.Osada", SBGEandSDIU.Body.tajD.S.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.tajD.S.by.SDIUxSBGE.head, "tajD.S", this.ylim = c(-0.18,-0.05))


#### models Tajima's D
SBGEandSDIU.Body.tajD.S.mod.data.complete.ExclExtSB.withPCs = main.df[ExclExtSB.Body.Osada.vec & Above4FoldMin.vec ,c("geneID","tajD.S","SDIU.body.p.strength.cat","log.avgExp.AdultLarva", "Whole.SBGE.Osada", 
                                                                                                                        "Is.X", "recomb", "totalExonLength", "tao.FML.FA2",
                                                                                                                        "m.Comp.1", "m.Comp.2", "m.Comp.3",
                                                                                                                        "f.Comp.1", "f.Comp.2", "f.Comp.3")]
SBGEandSDIU.Body.tajD.S.mod.data.complete.ExclExtSB.withPCs = na.omit(SBGEandSDIU.Body.tajD.S.mod.data.complete.ExclExtSB.withPCs)

options(contrasts = rep("contr.sum", 2)) ## this should make lmp and lm results comparable
mod.DPGP3.SDIU.tajD.S.withPCs = lm(tajD.S ~ SDIU.body.p.strength.cat 
                                     + log.avgExp.AdultLarva 
                                     + Whole.SBGE.Osada 
                                     + I(Whole.SBGE.Osada^2) 
                                     + Is.X + tao.FML.FA2  
                                     + log(recomb) + log(totalExonLength)
                                     + m.Comp.1 + m.Comp.2 + m.Comp.3
                                     + f.Comp.1 + f.Comp.2 + f.Comp.3, 
                                     data = SBGEandSDIU.Body.tajD.S.mod.data.complete.ExclExtSB.withPCs)
summary(mod.DPGP3.SDIU.tajD.S.withPCs)

mod.lmp.DPGP3.SDIU.tajD.S.withPCs = lmp(tajD.S ~ SDIU.body.p.strength.cat 
                                + log.avgExp.AdultLarva 
                                + Whole.SBGE.Osada 
                                + I(Whole.SBGE.Osada^2) 
                                + Is.X + tao.FML.FA2  
                                + log(recomb) + log(totalExonLength)
                                + m.Comp.1 + m.Comp.2 + m.Comp.3
                                + f.Comp.1 + f.Comp.2 + f.Comp.3,
                                center = FALSE,
                                perm = "Prob",
                                Ca = 1e-3,
                                maxIter = 1e6,
                                data = SBGEandSDIU.Body.tajD.S.mod.data.complete.ExclExtSB.withPCs)
summary(mod.lmp.DPGP3.SDIU.tajD.S.withPCs)


# Compute 95 % bootstrap confidence intervals
boot_case.mod.tajD.S.withPCs <- Boot(mod.DPGP3.SDIU.tajD.S.withPCs, method = "case", R = 9999)
confint(boot_case.mod.tajD.S.withPCs, type = "perc")


############################ end:  tajD.S ########################################################


############# piNfract  = piN / (piN + piS)

## Use permutation test to test if mean piNfract  
## is does not differ between SDIU genes vs non-SDIU genes

SDIU.piNfract.vec = main.df$piNfract[Above4FoldMin.vec & main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$piNfract)]
Not.SDIU.piNfract.vec = main.df$piNfract[Above4FoldMin.vec & !main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$piNfract)]

Diff.mean.piNfract.SDIUvsNotSDIU = mean(SDIU.piNfract.vec) - mean(Not.SDIU.piNfract.vec)
Diff.mean.piNfract.SDIUvsNotSDIU

SDIU.NotSDIU.piNfract.vec = c(SDIU.piNfract.vec, Not.SDIU.piNfract.vec)

Diff.mean.piNfract.SDIUvsNotSDIU.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.piNfract.vec, length(SDIU.piNfract.vec)))
sum(Diff.mean.piNfract.SDIUvsNotSDIU.null.dist <= Diff.mean.piNfract.SDIUvsNotSDIU)/length(Diff.mean.piNfract.SDIUvsNotSDIU.null.dist)
sum(Diff.mean.piNfract.SDIUvsNotSDIU.null.dist >= Diff.mean.piNfract.SDIUvsNotSDIU)/length(Diff.mean.piNfract.SDIUvsNotSDIU.null.dist)


SBGEandSDIU.Body.piNfract.for.simple.fig = main.df[Above4FoldMin.vec, ]

sorted.combined.summary.piNfract.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("piNfract","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.Body.piNfract.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.piNfract.by.SDIUxSBGE, "piNfract", this.ylim = c(0.05, 0.25))

### head fig
sorted.combined.summary.piNfract.by.SDIUxSBGE.head = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("piNfract","SBGEcat.head.Osada", "SDIUxSBGE.head.Osada", SBGEandSDIU.Body.piNfract.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.piNfract.by.SDIUxSBGE.head, "piNfract", this.ylim = c(0.05, 0.25))


#### models piNfract

SBGEandSDIU.Body.piNfract.mod.data.complete.ExclExtSB.withPCs = main.df[ExclExtSB.Body.Osada.vec & Above4FoldMin.vec ,c("geneID","piNfract","SDIU.body.p.strength.cat","log.avgExp.AdultLarva", "Whole.SBGE.Osada", 
                                                                                                                   "Is.X", "recomb", "totalExonLength", "tao.FML.FA2",
                                                                                                                   "m.Comp.1", "m.Comp.2", "m.Comp.3",
                                                                                                                   "f.Comp.1", "f.Comp.2", "f.Comp.3")]
SBGEandSDIU.Body.piNfract.mod.data.complete.ExclExtSB.withPCs = na.omit(SBGEandSDIU.Body.piNfract.mod.data.complete.ExclExtSB.withPCs)

options(contrasts = rep("contr.sum", 2)) ## this should make lmp and lm results comparable
mod.DPGP3.SDIU.piNfract.withPCs = lm(piNfract ~ SDIU.body.p.strength.cat 
                                + log.avgExp.AdultLarva 
                                + Whole.SBGE.Osada 
                                + I(Whole.SBGE.Osada^2) 
                                + Is.X + tao.FML.FA2  
                                + log(recomb) + log(totalExonLength)
                                + m.Comp.1 + m.Comp.2 + m.Comp.3
                                + f.Comp.1 + f.Comp.2 + f.Comp.3, 
                                data = SBGEandSDIU.Body.piNfract.mod.data.complete.ExclExtSB.withPCs)
summary(mod.DPGP3.SDIU.piNfract.withPCs)

mod.lmp.DPGP3.SDIU.piNfract.withPCs = lmp(piNfract ~ SDIU.body.p.strength.cat 
                                        + log.avgExp.AdultLarva 
                                        + Whole.SBGE.Osada 
                                        + I(Whole.SBGE.Osada^2) 
                                        + Is.X + tao.FML.FA2  
                                        + log(recomb) + log(totalExonLength)
                                        + m.Comp.1 + m.Comp.2 + m.Comp.3
                                        + f.Comp.1 + f.Comp.2 + f.Comp.3,
                                        center = FALSE,
                                        perm = "Prob",
                                        Ca = 1e-3,
                                        maxIter = 1e6,
                                        data = SBGEandSDIU.Body.piNfract.mod.data.complete.ExclExtSB.withPCs)
summary(mod.lmp.DPGP3.SDIU.piNfract.withPCs)


# Compute 95 % bootstrap confidence intervals:
boot_case.mod.piNfract.withPCs <- Boot(mod.DPGP3.SDIU.piNfract.withPCs, method = "case", R = 9999)
confint(boot_case.mod.piNfract.withPCs, type = "perc")



############################ end:  piNfract #######################################################

############# Direction of Selection DoS  ##########

## Use permutation test to test if mean DoS  
## differs between SDIU genes vs non-SDIU genes
SDIU.DoS.vec = main.df$DoS[Above4FoldMin.vec & main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$DoS)]
Not.SDIU.DoS.vec = main.df$DoS[Above4FoldMin.vec & !main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$DoS)]

Diff.mean.DoS.SDIUvsNotSDIU = mean(SDIU.DoS.vec) - mean(Not.SDIU.DoS.vec)
Diff.mean.DoS.SDIUvsNotSDIU

SDIU.NotSDIU.DoS.vec = c(SDIU.DoS.vec, Not.SDIU.DoS.vec)

Diff.mean.DoS.SDIUvsNotSDIU.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.DoS.vec, length(SDIU.DoS.vec)))
sum(Diff.mean.DoS.SDIUvsNotSDIU.null.dist <= Diff.mean.DoS.SDIUvsNotSDIU)/length(Diff.mean.DoS.SDIUvsNotSDIU.null.dist)
sum(Diff.mean.DoS.SDIUvsNotSDIU.null.dist >= Diff.mean.DoS.SDIUvsNotSDIU)/length(Diff.mean.DoS.SDIUvsNotSDIU.null.dist)



SBGEandSDIU.Body.DoS.for.simple.fig = main.df[Above4FoldMin.vec, ]
sorted.combined.summary.DoS.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("DoS","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.Body.DoS.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.DoS.by.SDIUxSBGE, "DoS", this.ylim = c(0.0, 0.3))



#### models DoS

SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs = main.df[ExclExtSB.Body.Osada.vec & Above4FoldMin.vec ,c("geneID","DoS","SDIU.body.p.strength.cat","log.avgExp.AdultLarva", "Whole.SBGE.Osada", 
                                                                                                                   "Is.X", "recomb", "totalExonLength", "tao.FML.FA2",
                                                                                                                   "m.Comp.1", "m.Comp.2", "m.Comp.3",
                                                                                                                   "f.Comp.1", "f.Comp.2", "f.Comp.3")]
SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs = na.omit(SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs)

options(contrasts = rep("contr.sum", 2)) ## this should make lmp and lm results comparable
mod.DPGP3.SDIU.DoS.withPCs = lm(DoS ~ SDIU.body.p.strength.cat 
                        + log.avgExp.AdultLarva 
                        + Whole.SBGE.Osada 
                        + I(Whole.SBGE.Osada^2) 
                        + Is.X + tao.FML.FA2  
                        + log(recomb) + log(totalExonLength)
                        + m.Comp.1 + m.Comp.2 + m.Comp.3
                        + f.Comp.1 + f.Comp.2 + f.Comp.3, 
                        data = SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs)
summary(mod.DPGP3.SDIU.DoS.withPCs)

mod.lmp.DPGP3.SDIU.DoS.withPCs = lmp(DoS ~ SDIU.body.p.strength.cat 
                                        + log.avgExp.AdultLarva 
                                        + Whole.SBGE.Osada 
                                        + I(Whole.SBGE.Osada^2) 
                                        + Is.X + tao.FML.FA2  
                                        + log(recomb) + log(totalExonLength)
                                        + m.Comp.1 + m.Comp.2 + m.Comp.3
                                        + f.Comp.1 + f.Comp.2 + f.Comp.3,
                                        center = FALSE,
                                        perm = "Prob",
                                        Ca = 1e-3,
                                        maxIter = 1e6,
                                        data = SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs)
summary(mod.lmp.DPGP3.SDIU.DoS.withPCs)


# Compute 95 % bootstrap confidence intervals:
boot_case.mod.DoS.withPCs <- Boot(mod.DPGP3.SDIU.DoS.withPCs, method = "case", R = 9999)
confint(boot_case.mod.DoS.withPCs, type = "perc")


############################ end:  DoS #############################################################

#### Remember that Var(DoS) = Var(DivPartDoS) + Var(PolyPartDoS) + 2*Cov(DivPartDoS, PolyPartDoS)
###  covariance matrix below shows that variance in both Div and Poly part contribute substantially
### to variance in DoS
SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs.DivAndPoly = main.df[ExclExtSB.Body.Osada.vec & Above4FoldMin.vec ,c("geneID","DoS","DivPartDoS","PolyPartDoS","SDIU.body.p.strength.cat","log.avgExp.AdultLarva", "Whole.SBGE.Osada", 
                                                                                                                   "Is.X", "recomb", "totalExonLength", "tao.FML.FA2",
                                                                                                                   "m.Comp.1", "m.Comp.2", "m.Comp.3",
                                                                                                                   "f.Comp.1", "f.Comp.2", "f.Comp.3")]
SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs.DivAndPoly = na.omit(SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs.DivAndPoly)

cov(SBGEandSDIU.Body.DoS.mod.data.complete.ExclExtSB.withPCs.DivAndPoly[,c("DoS","DivPartDoS","PolyPartDoS")])


############# divergence part of DoS

## Use permutation test to test if mean DivPartDoS  
## differs between SDIU genes vs non-SDIU genes
SDIU.DivPartDoS.vec = main.df$DivPartDoS[Above4FoldMin.vec & main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$DivPartDoS)]
Not.SDIU.DivPartDoS.vec = main.df$DivPartDoS[Above4FoldMin.vec & !main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$DivPartDoS)]

Diff.mean.DivPartDoS.SDIUvsNotSDIU = mean(SDIU.DivPartDoS.vec) - mean(Not.SDIU.DivPartDoS.vec)
Diff.mean.DivPartDoS.SDIUvsNotSDIU

SDIU.NotSDIU.DivPartDoS.vec = c(SDIU.DivPartDoS.vec, Not.SDIU.DivPartDoS.vec)

Diff.mean.DivPartDoS.SDIUvsNotSDIU.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.DivPartDoS.vec, length(SDIU.DivPartDoS.vec)))
sum(Diff.mean.DivPartDoS.SDIUvsNotSDIU.null.dist <= Diff.mean.DivPartDoS.SDIUvsNotSDIU)/length(Diff.mean.DivPartDoS.SDIUvsNotSDIU.null.dist)
sum(Diff.mean.DivPartDoS.SDIUvsNotSDIU.null.dist >= Diff.mean.DivPartDoS.SDIUvsNotSDIU)/length(Diff.mean.DivPartDoS.SDIUvsNotSDIU.null.dist)


## DivPartDoS is elevated at extreme SBGE categories 
SBGEandSDIU.Body.DivPartDoS.for.simple.fig = main.df[Above4FoldMin.vec, ]

sorted.combined.summary.DivPartDoS.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("DivPartDoS","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.Body.DivPartDoS.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.DivPartDoS.by.SDIUxSBGE, "DivPartDoS", this.ylim = c(0.2, 0.5))




SBGEandSDIU.Body.DivPartDoS.mod.data.complete.ExclExtSB.withPCs = main.df[ExclExtSB.Body.Osada.vec & Above4FoldMin.vec ,c("geneID","DivPartDoS","SDIU.body.p.strength.cat","log.avgExp.AdultLarva", "Whole.SBGE.Osada", 
                                                                                                                   "Is.X", "recomb", "totalExonLength", "tao.FML.FA2",
                                                                                                                   "m.Comp.1", "m.Comp.2", "m.Comp.3", 
                                                                                                                   "f.Comp.1", "f.Comp.2", "f.Comp.3")]
SBGEandSDIU.Body.DivPartDoS.mod.data.complete.ExclExtSB.withPCs = na.omit(SBGEandSDIU.Body.DivPartDoS.mod.data.complete.ExclExtSB.withPCs)

options(contrasts = rep("contr.sum", 2)) ## this should make lmp and lm results comparable
mod.DPGP3.SDIU.DivPartDoS.withPCs = lm(DivPartDoS ~ SDIU.body.p.strength.cat 
                                + log.avgExp.AdultLarva 
                                + Whole.SBGE.Osada 
                                + I(Whole.SBGE.Osada^2) 
                                + Is.X + tao.FML.FA2  
                                + log(recomb) + log(totalExonLength)
                                + m.Comp.1 + m.Comp.2 + m.Comp.3
                                + f.Comp.1 + f.Comp.2 + f.Comp.3, 
                                data = SBGEandSDIU.Body.DivPartDoS.mod.data.complete.ExclExtSB.withPCs)
summary(mod.DPGP3.SDIU.DivPartDoS.withPCs)

mod.lmp.DPGP3.SDIU.DivPartDoS.withPCs = lmp(DivPartDoS ~ SDIU.body.p.strength.cat 
                                     + log.avgExp.AdultLarva 
                                     + Whole.SBGE.Osada 
                                     + I(Whole.SBGE.Osada^2) 
                                     + Is.X + tao.FML.FA2  
                                     + log(recomb) + log(totalExonLength)
                                     + m.Comp.1 + m.Comp.2 + m.Comp.3
                                     + f.Comp.1 + f.Comp.2 + f.Comp.3,
                                     center = FALSE,
                                     perm = "Prob",
                                     Ca = 1e-3,
                                     maxIter = 1e6,
                                     data = SBGEandSDIU.Body.DivPartDoS.mod.data.complete.ExclExtSB.withPCs)
summary(mod.lmp.DPGP3.SDIU.DivPartDoS.withPCs)


# Compute 95 % bootstrap confidence intervals:
boot_case.mod.DivPartDoS.withPCs <- Boot(mod.DPGP3.SDIU.DivPartDoS.withPCs, method = "case", R = 9999)
confint(boot_case.mod.DivPartDoS.withPCs, type = "perc")

############################ end:  DivPartDoS #######################################################

############# polymorphism part of DoS: Pn/(Pn + Ps)

## Use permutation test to test if mean PolyPartDoS  
## differs between SDIU genes vs non-SDIU genes
SDIU.PolyPartDoS.vec = main.df$PolyPartDoS[Above4FoldMin.vec & main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$PolyPartDoS)]
Not.SDIU.PolyPartDoS.vec = main.df$PolyPartDoS[Above4FoldMin.vec & !main.df$SDIU.body.sig & !is.na(main.df$SDIU.body.sig) & !is.na(main.df$PolyPartDoS)]

Diff.mean.PolyPartDoS.SDIUvsNotSDIU = mean(SDIU.PolyPartDoS.vec) - mean(Not.SDIU.PolyPartDoS.vec)
Diff.mean.PolyPartDoS.SDIUvsNotSDIU

SDIU.NotSDIU.PolyPartDoS.vec = c(SDIU.PolyPartDoS.vec, Not.SDIU.PolyPartDoS.vec)

Diff.mean.PolyPartDoS.SDIUvsNotSDIU.null.dist = sapply(1:10000, function(x) get.one.perm.diff(SDIU.NotSDIU.PolyPartDoS.vec, length(SDIU.PolyPartDoS.vec)))
sum(Diff.mean.PolyPartDoS.SDIUvsNotSDIU.null.dist <= Diff.mean.PolyPartDoS.SDIUvsNotSDIU)/length(Diff.mean.PolyPartDoS.SDIUvsNotSDIU.null.dist)
sum(Diff.mean.PolyPartDoS.SDIUvsNotSDIU.null.dist >= Diff.mean.PolyPartDoS.SDIUvsNotSDIU)/length(Diff.mean.PolyPartDoS.SDIUvsNotSDIU.null.dist)


## PolyPartDoS is elevated at extreme SBGE categories 
SBGEandSDIU.Body.PolyPartDoS.for.simple.fig = main.df[Above4FoldMin.vec, ]

sorted.combined.summary.PolyPartDoS.by.SDIUxSBGE = MakeSummaryTableForSDIUxSBGE.forSpecificVariable("PolyPartDoS","SBGEcat.body.Osada", "SDIUxSBGE.body.Osada", SBGEandSDIU.Body.PolyPartDoS.for.simple.fig)
Make.Figure.SDIUxSBGE.ForSpecifiedVariable.Body(sorted.combined.summary.PolyPartDoS.by.SDIUxSBGE, "PolyPartDoS", this.ylim = c(0.05, 0.25))


#### models PolyPartDoS

SBGEandSDIU.Body.PolyPartDoS.mod.data.complete.ExclExtSB.withPCs = main.df[ExclExtSB.Body.Osada.vec & Above4FoldMin.vec ,c("geneID","PolyPartDoS","SDIU.body.p.strength.cat","log.avgExp.AdultLarva", "Whole.SBGE.Osada", 
                                                                                                                          "Is.X", "recomb", "totalExonLength", "tao.FML.FA2",
                                                                                                                          "m.Comp.1", "m.Comp.2", "m.Comp.3",
                                                                                                                          "f.Comp.1", "f.Comp.2", "f.Comp.3")]
SBGEandSDIU.Body.PolyPartDoS.mod.data.complete.ExclExtSB.withPCs = na.omit(SBGEandSDIU.Body.PolyPartDoS.mod.data.complete.ExclExtSB.withPCs)

options(contrasts = rep("contr.sum", 2)) ## this should make lmp and lm results comparable
mod.DPGP3.SDIU.PolyPartDoS.withPCs = lm(PolyPartDoS ~ SDIU.body.p.strength.cat 
                                       + log.avgExp.AdultLarva 
                                       + Whole.SBGE.Osada 
                                       + I(Whole.SBGE.Osada^2) 
                                       + Is.X + tao.FML.FA2  
                                       + log(recomb) + log(totalExonLength)
                                       + m.Comp.1 + m.Comp.2 + m.Comp.3
                                       + f.Comp.1 + f.Comp.2 + f.Comp.3, 
                                       data = SBGEandSDIU.Body.PolyPartDoS.mod.data.complete.ExclExtSB.withPCs)
summary(mod.DPGP3.SDIU.PolyPartDoS.withPCs)

mod.lmp.DPGP3.SDIU.PolyPartDoS.withPCs = lmp(PolyPartDoS ~ SDIU.body.p.strength.cat 
                                            + log.avgExp.AdultLarva 
                                            + Whole.SBGE.Osada 
                                            + I(Whole.SBGE.Osada^2) 
                                            + Is.X + tao.FML.FA2  
                                            + log(recomb) + log(totalExonLength)
                                            + m.Comp.1 + m.Comp.2 + m.Comp.3
                                            + f.Comp.1 + f.Comp.2 + f.Comp.3,
                                            center = FALSE,
                                            perm = "Prob",
                                            Ca = 1e-3,
                                            maxIter = 1e6,
                                            data = SBGEandSDIU.Body.PolyPartDoS.mod.data.complete.ExclExtSB.withPCs)
summary(mod.lmp.DPGP3.SDIU.PolyPartDoS.withPCs)


# Compute 95 % bootstrap confidence intervals:
boot_case.mod.PolyPartDoS.withPCs <- Boot(mod.DPGP3.SDIU.PolyPartDoS.withPCs, method = "case", R = 9999)
confint(boot_case.mod.PolyPartDoS.withPCs, type = "perc")



############################ end:  PolyPartDoS  #######################################################




