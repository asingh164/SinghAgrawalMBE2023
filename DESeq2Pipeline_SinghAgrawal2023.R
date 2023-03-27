#############################################################################
######              Differential Gene Expression Pipeline               #####
#############################################################################

## Author: Amardeep Singh -- amardeep.singh[at]alum.utoronto.ca
## This script makes use of SAM files that I created separately see: STARAlignmentPipeline_SinghAgrawal2023 script

### Script Details ###
# This script is written primarily in R with bash commands.
# R commands are (sometimes) denoted by #---- R ---- and end with # ----/----
# Bash commands are denoted by #---- BASH ---- and end with # ----/----
# The script only refers to samples taken from bodies, but head samples were also run with this pipelinne

##
# Path to sam files: /PATH/TO/SAM/FILES/
# Path to gtf /PATH/TO/GTF/

#####################################################################
###     Creating count matrices of each aligned transcriptome     ###
#####################################################################
require("GenomicFeatures")
require("Rsamtools")
require("GenomicAlignments")
require("DESeq2")

#############################
###   GENOMIC FEATURES    ###
#############################
# Read in gene model from GFF file
gtf.TxDb = makeTxDbFromGFF("/PATH/TO/GENE/MODEL/FILE", format="gtf") #A gff or gtf file is fine

# Creates a list of exons grouped by gene
exonsByGene = exonsBy(gtf.TxDb, by="gene")

# Specify BAM files to be used for generating read count matrix  by pointing to the correct directory
# I specified the pattern of the file names, but you might need to change this or not include it
bam_files.body = list.files("/PATH/TO/DIRECTORY/WITH/BAM/FILES", pattern="SAM.sorted.trimmed.RAL", full=TRUE)

# Indicate that bamLst are BAM files
bamLst.body = BamFileList(bam_files.body, yieldSize=100000)

###############################
###   GENOMIC ALIGNMENTS    ###
###############################
# Call sumarizeOvarlaps to count reads and produce sum_exp which gets inputted into DESeq2
sum_exp.body <- summarizeOverlaps(exonsByGene, bamLst.body, mode="Union", singleEnd=TRUE, ignore.strand=TRUE)

############################################
###      Prepping files for DESEQ2       ###
############################################
# Differential gene expression analysis pipeline
# Provide metadata about the individuals found the count matrix
# (each column of matrix = individual, each row of RNA_Data = info about individual)
# Making sure order of individuals in RNA Data and in summarized experiment match

# Make a condition table
# Body
condition.table.body = as.data.frame(factor(colnames(sum_exp.body)))
colnames(condition.table.body) = "SampleName"
condition.table.body$condition = factor(rep(c("female", "male"), 18))
condition.table.body$condition = as.factor(condition.table.body$condition)
#condition.table = as.data.frame(condition.table)
colData(sum_exp.body) = DataFrame(condition.table.body)

###########################################
###     Running the DESeq2 pipeline     ###
###########################################

# build the DESeqDataSet from a SummarizedExperiment

# Whole Body
dds_interaction.body = DESeqDataSet(sum_exp.body, design = ~ condition)
# Filter out genes with fewer than 50 reads mapping across all samples (more stringent filtering is applied downstream)
keep.body = rowMeans(counts(dds_interaction.body)) >= 50
dds_interaction.body = dds_interaction.body[keep.body,]
#run DESeq2
dds_interaction.body = DESeq(dds_interaction.body)
#name results
all_genes_interaction.body = results(dds_interaction.body, contrast=c("condition", "male", "female")) # This specifies the output such that the log2FC is male/female

# FlyBaseID as first column
all_genes_interaction.body$FlyBaseID= rownames(all_genes_interaction.body)
all_genes_interaction.body = all_genes_interaction.body[,c(7,1:6)]

# Save output files for both whole body and head
write.table(all_genes_interaction.body, file = "/PATH/TO/OUTPUT/FILE/OUTPUTFILE.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Save DESeq2 RDS files
saveRDS(dds_interaction.body, file = "/PATH/TO/OUTPUT/FILE/OUTPUTFILE.rds")



##
