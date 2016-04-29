#Package ChIPseeker
library('ChIPseeker')
library('GenomicRanges')

browseVignettes("ChIPseeker")
browseVignettes("GenomicRanges")


files_example <- getSampleFiles()
peak_example <- readPeakFile(files_example[[4]])
covplot(peak_example, weightCol="X8.44")

peak_example=GenomicRanges::GRangesList(CBX6=readPeakFile(files_example[[4]]), CBX7=readPeakFile(files_example[[5]]))
p <- covplot(peak_example)
print(p)


#Add all files from GSE68349 directory:
        GSE68349_path <- "~/Desktop/PhD/Paper analysis/Dominguez-Sola et al., 2015/GSE68349_RAW"
   GSE68349_all_files <- list.files(path = GSE68349_path)

#Remove files with '.gz' extension:
 GSE68349_files_no.gz <- GSE68349_all_files[grep("gz",GSE68349_all_files, invert = T)]
       GSE68349_files <- as.list(paste(GSE68349_path, "/", GSE68349_files_no.gz, sep = ""))
names(GSE68349_files) <- c("RK059","RK051","RK040","RK050")

peak_RK059 <- readPeakFile(GSE68349_files[["RK059"]], header = F)
peak_RK051 <- readPeakFile(GSE68349_files[["RK051"]], header = F)
peak_RK040 <- readPeakFile(GSE68349_files[["RK040"]], header = F)
peak_RK050 <- readPeakFile(GSE68349_files[["RK050"]], header = F)

peak_files <- as.list(peak_RK059,peak_RK051,peak_RK040,peak_RK050)
lapply(GSE68349_files, covplot, weightCol = "V4")

peak<-c(RK059=readPeakFile(files[[1]]), RK051=readPeakFile(files[[2]]))
covplot(peak, weightCol = "V4", chrs = c("chr14", "chr16"))
#or specify the full path:
#peak_RK059 <- readPeakFile("~/Desktop/PhD/Paper analysis/Dominguez-Sola et al., 2015/GSE68349_RAW/GSM1668935_Peaks.CB4_FOXO1_RK059_vs_Input_RK063_p5.bedgraph.gz")


#Plot peaks:
covplot(peak_RK059, weightCol = "V4")

covplot(peak_RK051, weightCol = "V4", chrs = c("chr14", "chr22"), #xlim=c(4.5e7, 5e7)
        )
#chr14 = human Igh locus
#chr22 = human Igl locus

# 4.2 Profile of ChIP peaks binding to TSS regions:
#First we need to create the vector txdb. Call TxDb... package
library('TxDb.Hsapiens.UCSC.hg19.knownGene')

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

promoter <- getPromoters(TxDb = txdb, upstream = 5000, downstream = 4000)
tagMatrix <- getTagMatrix(peak_RK059, windows = promoter)

tagHeatmap(tagMatrix, xlim = c(-5000, 4000), color = "red")

#ChIPseeker provide a one step function to generate this figure from bed file.
#The following function will generate the same figure as above:
peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000, color="red")

plotAvgProf(tagMatrix, xlim=c(-5000, 4000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

#Add confident interval to previous function:
#plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 500)

#Peak annotation:
peakAnno <- annotatePeak(peak_RK059, tssRegion = c(-5000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
head(as.GRanges(peakAnno))

peakAnno_dataframe <- as.data.frame(peakAnno)
    peakAnno_chr14 <- peakAnno_dataframe[grep('chr14', peakAnno_dataframe$seqnames),]


# Profile of several ChIP peak data binding to TSS region:
     promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrixList <- lapply(GSE68349_files, getTagMatrix, windows = promoter)

 
plotAvgProf(tagMatrixList, xlim = c(-3000, 3000))

plotAvgProf(tagMatrixList, xlim = c(-3000, 3000), conf=0.95,resample=100, facet="row")



