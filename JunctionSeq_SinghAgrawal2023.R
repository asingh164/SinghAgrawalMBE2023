nSeq Script for Differential Exon Usage Analysis  #####
#####################################################################

## Author: Amardeep Singh -- amardeep.singh[at]alum.utoronto.ca
## This program makes use of public releases of the Drosophila genome from Ensembl (BDGP6.28)
## To retrieve the genome assembly file yourself:
# Fasta file
# wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/drosophila_melanogaster/dna//Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz
# GTF file
# wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/drosophila_melanogaster/

# This script requires the 'JunctionSeq' R packages which can be downloaded through Bioconductor and an older version of
### NOTE: Junctonseq (as of July 2020) will require that you install an older version of DESeq2 (version 1.10.1)
##                https://www.bioconductor.org/packages/3.2/bioc/src/contrib/DESeq2_1.10.1.tar.gz
## NOTE: This same script was used to independently examine differential exon usage in both whole body and brain tissue but the code only refers to a single tissue.

### Script Details ###
# This script is written primarily in R, Bash commands are denoted by #---- Linux ---- and end with # ---- ----

# ---- LINUX ----
# Export PATHS
export PATH_TO_GFF=/PATH/TO/GFF/Drosophila_melanogaster.BDGP6.28.99.gtf.gz
export OUT_PATH_FLAT_GFF=/OUTPUT/FILE/PATH/HERE

# Make a flat gff file from the genome annotation provided by Ensembl
QoRTs makeFlatGff --stranded \
  ${PATH_TO_GFF} \
  ${OUT_PATH_FLAT_GFF}/Dmel.JunctionSeq.flat.gff.gz

## Processing bam files
# First, run QoRTs to perform QC checks on aligned bam files. This needs to be done so that the size factor for normalization can be done
# NOTE: QoRTs seems to be pretty memory hungry. I've allocated java 10gb of memory here, you may need to change this depending on available resources
ls *.bam | parallel -j 50 \
"mkdir {.} && java -Xmx10G -jar /plas1/amardeep.singh/apps/QoRTs-STABLE.jar \
                    --stranded --maxReadLength 101 \
                    --singleEnded {} ${PATH_TO_GFF} {.}/"

# Make a 'decoder' file that will use the sample names as the directories to find outputs
# I had these in the following folder, but this might change depending on where and what is in the directory.
ls /PATH/TO/SAMPLES/BAM/QoRTS/OUTPUT/FILES > QoRTs.decoder.file.txt
# ----/LINUX----

# ---- R Code ---
require(QoRTs)

# Read QoRTs outputs into R to generate a sizeFactor by which to normalize read counts
# Tell QoRTs which directory to look into for directory that contain outputs. Each directory should be named to correspond to a sample (which should be decoded by the decoder file) and give it the location of your decoder file
qorts.results <- read.qc.results.data("PATH/TO/qorts.output/",
                                      decoder.files = "/PATH/TO/DECODER/FILE/",
                                      calc.DESeq2 = TRUE)

# Save size factors for each sample in a text file
get.size.factors(qorts.results, outfile = "/PATH/TO/SAVE/OUTPUTS/SIZEFACTOR_FILE_NAME.txt");

# ----/R Code ---

# ---- Linux ----
# Create novel junction splice sites
java -Xmx10G -jar /plas1/amardeep.singh/apps/QoRTs-STABLE.jar mergeNovelSplices --minCount 20 --stranded \
/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/RAL.mapped.files.to.ensembl/bam.files/samtools.sort.out/qorts.output \
/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/RAL.size.Factors.GEO.txt \
/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/gtf/Drosophila_melanogaster.BDGP6.28.99.gtf.gz \
/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files
# ----/LINUX----

# Running JunctionSeq in R

# ---- R Code ----

### JunctionSeq Script for Body Tissue

rm(list=ls())
require(DESeq2)
require(JunctionSeq)
require(BiocParallel) # This is a package used for paralellization of jobs

# Set global variables
numCores = 50
FDRThreshold = 0.01
mappedReadsThreshold = 50


# Load in decoder files and add fields for conditions

# Run this only once to edit the decoder file
#decoder.file =  read.delim("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.txt", header = FALSE, stringsAsFactors=FALSE)
#colnames(decoder.file) = "unique.ID"
#sex.tmp = gsub("SAM.sorted.trimmed.RAL.....","", decoder.file$unique.ID)
#decoder.file$sex =  gsub("\\..*","", sex.tmp)
#write.table(decoder.file, "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", quote = F, row.names = FALSE, col.names = TRUE, sep = "\t")

# Loading in decoder file
decoder <- read.table("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/QoRTs.decoder.file.for.JunctionSeq.txt", header = TRUE, stringsAsFactors = FALSE)
# Here I am subsetting only the samples that I want to use (for this analysis I started just with body tissue and only a single female replicate)
decoder.for.junctionseq <- decoder[!(grepl("head|replicate.2", decoder$unique.ID)),]

#Providing the directory for the count files:
countFiles <- paste0("/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files/body.only.replicate.1/", decoder.for.junctionseq$unique.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

###  Run the differential exon usage (DEU) analysis:
# Creast a design dataframe which has a column for the condition you are testng (This will need to change when adding factors not sure how yet)
design.df <- data.frame(condition = factor(decoder.for.junctionseq$sex))

# Building the count set object that JunctionSeq will analyze and add to it all of the parameters of the analysis
count.set.object <- readJunctionSeqCounts(countfiles = countFiles,
                                          samplenames = decoder.for.junctionseq$unique.ID,
                                          design = design.df,
                                          flat.gff.file = "/plas1/amardeep.singh/Ensembl.Dmel.Genome.Release/JunctionSeq.files/count.files/withNovel.forJunctionSeq.gff.gz",
                                          nCores = numCores,
                                          verbose = TRUE)

# Generate size factors for normalization and load them into the count.set.object
count.set.object <- estimateJunctionSeqSizeFactors(count.set.object)

# Generate test specific dispersion estimates and load into count.set.object
count.set.object <- estimateJunctionSeqDispersions(count.set.object, nCores = numCores)

# Fit the observed dispersions to a regression to create a fitted dispersion
count.set.object <- fitJunctionSeqDispersionFunction(count.set.object)

# Perform the hypothesis tests to test for differential splice junction/exon usage (DEU)
count.set.object <- testForDiffUsage(count.set.object, nCores = numCores)

# Calculate effect sizes and parameter estimates
count.set.object <- estimateEffectSizes(count.set.object)

# Save output to file
writeCompleteResults(count.set.object, "/PATH/TO/OUTPUT/file.txt",
                    gzip.output = TRUE,
                    FDR.threshold = FDR_Threshold,
                    save.allGenes = TRUE, save.sigGenes = TRUE,
                    save.fit = FALSE, save.VST = FALSE,
                    save.bedTracks = TRUE,
                    save.jscs = TRUE,
                    bedtrack.format = c("BED", "GTF", "GFF3"),
                    verbose = TRUE)
