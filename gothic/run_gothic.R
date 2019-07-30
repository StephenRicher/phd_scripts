
library(BSgenome.HSapienscustom.ensembl.hg38)
library('GOTHiC')

#forgeBSgenomeDataPkg("bsgenome_seed_file.txt")

custom_genome = BSgenome.HSapienscustom.ensembl.hg38
#bam_R1 = "/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment/custom_bam/HB2_CL4_1.R1.GNG12_AS1_DIRAS3.bam"
#bam_R2 = "/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment/custom_bam/HB2_CL4_1.R2.GNG12_AS1_DIRAS3.bam"
bam_R1 = "/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment/custom_bam/HB2_CL4_1.R1.captured.bam"
bam_R2 = "/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment/custom_bam/HB2_CL4_1.R2.captured.bam"
sample = "HB2_CL4_1"

results = GOTHiC(fileName1 = bam_R1, 
       fileName2 = bam_R2, 
       sampleName = sample, 
       res = 5000,
       BSgenomeName='BSgenome.HSapienscustom.ensembl.hg38',
       genome = BSgenome.HSapienscustom.ensembl.hg38, 
       restrictionSite='^GATC',
       enzyme='MBo1',
       cistrans='cis',
       filterdist = 0,
       DUPLICATETHRESHOLD = 1, 
       fileType='BAM', 
       parallel = FALSE, 
       cores = NULL)
