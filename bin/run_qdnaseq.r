#!/usr/bin/env Rscript

# Title     : Run QDNAseq analusis
# Objective : Perform CN analysis of input bam file
# Created by: rjf
# Created on: 05/08/21
# INCOMPLETE, see documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf


library(argparser)

p <- arg_parser("Run a CN anlysis via QDNAseq for hg38")
p <- add_argument(p, "--bam", help="Input BAM file (ref hg38)")
p <- add_argument(p, "--out_prefix", help="Prefix for output file names")
p <- add_argument(p, "--qdnaseq_hg38", help="HG38 reference set R package (zip)")
p <- add_argument(p, "--method", help="cutoff or CGHcall", default="cutoff")
p <- add_argument(p, "--cutoff", help="CN cutoff modifer, 0.0-1.0", default='0.5')
p <- add_argument(p, "--cutoffDEL", help="CN cutoff threshold for deletion (double loss)", default=0.5)
p <- add_argument(p, "--cutoffLOSS", help="CN cutoff threshold for loss", default=1.5)
p <- add_argument(p, "--cutoffGAIN", help="CN cutoff threshold for gain", default=2.5)
p <- add_argument(p, "--cellularity", help="CGHcall cellularity", default=1.0)
p <- add_argument(p, "--reference", help="QDNAseq GC/mappability bins reference. Defaults to qdnaseq_hg38", default="hg38")
p <- add_argument(p, "--binsize", help="bin size in MBp", default=500)

argv <- parse_args(p)
if (argv$cutoff %in% c('none','None','NONE')) {
    argv$cutoff <- 'none'
} else {
    argv$cutoff <- as.numeric(as.character(argv$cutoff))
}

options(future.globals.maxSize = 1048576000)

#Set-up fig outputs
pdf_file <- paste(argv$out_prefix, 'plots.pdf', sep="_")
pdf(pdf_file)

#Import QDNAseq lib
library(QDNAseq)


if (argv$reference == "hg38") {
  library(QDNAseq.hg38)
  bins <- getBinAnnotations(binSize=argv$binsize, genome="hg38")
} else if (argv$reference == "hg19") {
    library(QDNAseq.hg19)
  bins <- getBinAnnotations(binSize=argv$binsize, genome="hg19")
}


## Primary analysis

#Import BAM files
readCounts <- binReadCounts(bins, bamfiles=argv$bam)

#Plot a raw copy number profile (read counts across the genome)
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))

#Highlight bins that will be removed with default filtering
highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)

#Apply filters and plot median read counts as a function of GC content and mappability
autosomalReadCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

bedout <- paste(argv$out_prefix, "raw_bins.bed", sep="_")
exportBins(readCounts, file=bedout, format="bed", type="copynumber")

#Plot GC / mappability / readcount
isobarPlot(autosomalReadCountsFiltered)

#Estimate the GC / mappability correction
autosomalReadCountsFiltered <- estimateCorrection(autosomalReadCountsFiltered)

#Noise plot
noisePlot(autosomalReadCountsFiltered)

#Create copy numbers object
readCountsFiltered <- applyFilters(autosomalReadCountsFiltered, chromosomes=NA)

#Apply correction for GC content and mappability
copyNumbers <- correctBins(readCountsFiltered)

#Normalize to median
copyNumbersNormalized <- normalizeBins(copyNumbers)

#Smooth outliers
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

#Plot smoothed copy number profile
plot(copyNumbersSmooth)

#Export bins
bedout <- paste(argv$out_prefix, "bins.bed", sep="_")
exportBins(copyNumbersSmooth, bedout, format="bed") #ADD VARIABLE NAME


##Secondary analysis

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)

if (argv$method=='cutoff') {
    if (argv$cutoff == 'none') {
        copyNumbersCalled <- callBins(copyNumbersSegmented, method = 'cutoff')
    }
    else {
        cutoffDEL <- argv$cutoffDEL + 0.5 - argv$cutoff
        cutoffLOSS <- argv$cutoffLOSS + 0.5 - argv$cutoff
        cutoffGAIN <- argv$cutoffGAIN + argv$cutoff - 0.5
        copyNumbersCalled <- callBins(copyNumbersSegmented, method = 'cutoff', cutoffs=log2(c(deletion = cutoffDEL, loss = cutoffLOSS, gain = cutoffGAIN, amplification = 10)/2))
    }
}
if (argv$method=='CGHcall') {
    copyNumbersCalled <- callBins(copyNumbersSegmented, method = 'CGHcall', cellularity=argv$cellularity )
}
plot(copyNumbersCalled)

#Create PNG output
png_file <- paste(argv$out_prefix, 'cov.png', sep="_")
png(png_file)
plot(copyNumbersCalled)

noise_png_file <- paste(argv$out_prefix, 'noise_plot.png', sep="_")
png(noise_png_file)
noisePlot(autosomalReadCountsFiltered)

isobar_png_file <- paste(argv$out_prefix, 'isobar_plot.png', sep="_")
png(isobar_png_file)
isobarPlot(autosomalReadCountsFiltered)

#write outputs
bedout <- paste(argv$out_prefix, "calls.bed", sep="_")
exportBins(copyNumbersCalled, file=bedout, format="bed", type="calls")
vcfout <- paste(argv$out_prefix, "calls.vcf", sep="_")
exportBins(copyNumbersCalled, file=vcfout, format="vcf", type="calls")

segout <- paste(argv$out_prefix, "segs.seg", sep="_")
exportBins(copyNumbersCalled, file=segout, format="seg", type="segments")
segout <- paste(argv$out_prefix, "segs.bed", sep="_")
exportBins(copyNumbersCalled, file=segout, format="bed", type="segments")
segout <- paste(argv$out_prefix, "segs.vcf", sep="_")
exportBins(copyNumbersCalled, file=segout, format="vcf", type="segments")
