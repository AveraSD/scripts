###############################################################################
## PANC-CNV
###############################################################################
# Copyright (c) 2016 Tobias Meissner

# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
# THE SOFTWARE.

#!/usr/bin/env Rscript

###############################################################################
# command line options
###############################################################################
packages <- function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x, character.only=TRUE, quietly=TRUE)){
    install.packages(pkgs=x, repos="http://cran.r-project.org", quiet=TRUE)
    require(x, character.only=TRUE, quietly=TRUE)
  }
}
packages(optparse)
option_list = list(
  make_option(c("-t", "--tumor"), type="character", default=NULL, 
              help="sambamba depth tumor output", metavar="character"),
  make_option(c("-g", "--germline"), type="character", default=NULL, 
              help="sambamba depth germline output(s)", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="bed file", metavar="character"),
  make_option(c("-a", "--amplification_cut"), type="numeric", default=1.5, 
              help="Cutoff to define amlification [default= %default]", 
              metavar="numeric"),
  make_option(c("-d", "--deletion_cut"), type="numeric", default=0.5, 
              help="Cutoff to define deletion [default= %default]", 
              metavar="numeric"),
  make_option(c("-c", "--min_cov"), type="numeric", default=100, 
              help="min. coverage for CNV call [default= %default]", 
              metavar="numeric"),
  make_option(c("-p", "--tumor_purity"), type="numeric", default=1, 
              help="purity of tumor sample [default= %default]", 
              metavar="numeric"),
  make_option(c("-o", "--out"), type="character", default="cnv.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$tumor)){
  print_help(opt_parser)
  stop("A tumor sample must be supplied", call.=FALSE)
}
if (is.null(opt$germline)){
  print_help(opt_parser)
  stop("(A) germline sample(s) must be supplied", call.=FALSE)
}
if (is.null(opt$bed)){
  print_help(opt_parser)
  stop("A bed file must be supplied", call.=FALSE)
}

# if there are multiple germline input files, split by seperator (comma)
opt$germline <- unlist(strsplit(opt$germline, split=','))

###############################################################################
# load and install required packages if needed
###############################################################################
packages(plyr)
packages(GenomicRanges)
packages(SomaticSignatures)
packages(BSgenome.Hsapiens.UCSC.hg19)
packages(diffloop)

###############################################################################
# functions
###############################################################################

MyMerge <- function(x, y){
  df <- merge(x, y, by= "ID", all=T)
  return(df)
}

# adjust GC content
adjustGC <- function(granges, id, plot=TRUE) {
  gc <- gcContent(granges, BSgenome.Hsapiens.UCSC.hg19)
  count <- values(granges)[,id]
  rough <- loess(count ~ gc, span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final <- loess(predict(rough, i) ~ i, span = 0.3)
  normv <- predict(final, gc)
  
  if(plot) {
    plot(count ~ gc, 
         ylim = quantile(count, c(1e-04, 0.999),na.rm=T), 
         xlim = c(0, 1), 
         pch = ".")
    points(count ~ gc, col = rgb(1, 0, 0, 0.3), pch = ".")
    lines(i, predict(rough, i), col = "green")
    points(gc, normv, col = "red", pch = ".")
  }
  
  countgcloess <- count/(normv/median(normv, na.rm = TRUE))
  countgcloess
}

dpToM <- function(dp, files) {
  dpM <- Reduce(MyMerge, dp)
  colnames(dpM)[-1] <- basename(files)
  # remove duplicated IDs if they exist
  if(any(duplicated(dpM$ID))) {
    dpM <- dpM[-which(duplicated(dpM$ID)),]
  }
  rownames(dpM) <- dpM[,1]
  dpTumorM <- dpM[,-1]
  dpM[is.na(dpM)] <- 0
  dpM
}

dpMToGR <- function(dpM) {
  dd <- do.call(rbind, strsplit(rownames(dpM), '_'))
  dd <- cbind(dd, dd[,2])
  dpM <- cbind(dd, dpM)
  colnames(dpM)[1:3] <- c('chrom', 'start', 'end')
  dpM <- as.data.frame(dpM, stringsAsFactors=F)
  dpM$start <- as.numeric(as.vector(dpM$start))
  dpM$end <- as.numeric(as.vector(dpM$end))
  grdpM <- as(as.data.frame(dpM), 'GRanges')
  grdpM
}

myBin <- function(grBed, grTG) {
  overlaps <- findOverlaps(grBed, grTG)
  samples <- colnames(mcols(grTG))
  sampleList <- list()
  for (i in 1:length(samples)) {
    sampleList[[i]] <- sapply(
      splitAsList(mcols(grTG)[[samples[i]]][subjectHits(overlaps)], 
                  factor(queryHits(overlaps), seq_len(length(grBed)))), 
      as.numeric)
    #sampleList[[i]] <- round(sapply(sampleList[[i]], mean),0)
    sampleList[[i]] <- round(sapply(sampleList[[i]], median),0)
  }
  tt <- DataFrame(sampleList)
  colnames(tt) <- samples
  gr <- grBed
  values(gr) <- cbind(values(gr), tt)
  gr
}

cnvPat <- function(grCNVp) {
  pair <- grCNVp[ , c(1, which(colnames(mcols(grCNVp))=='TUMOR'), 
                      which(colnames(mcols(grCNVp))=='NORMAL'))]
  if(any(apply(as.matrix(values(pair)[, 2:3]), 1, function(x) any(x==0 | is.na(x))))) {
    pair <- pair[-which(apply(as.matrix(values(pair)[, 2:3]), 1, 
                              function(x) any(x==0 | is.na(x)))),]
  }
  if(any(apply(as.matrix(values(pair)[, 2:3]), 1, function(x) any(is.nan(x))))) {
    pair <- pair[-which(apply(as.matrix(values(pair)[, 2:3]), 1, 
                              function(x) any(is.nan(x)))),]  
  }
  if(any(apply(as.matrix(values(pair)[, 2:3]), 1, function(x) any(is.na(x))))) {
    pair <- pair[-which(apply(as.matrix(values(pair)[, 2:3]), 1, 
                              function(x) any(is.na(x)))),]
  }
  ratio <- round(2^(log2(mcols(pair)[,2]) - log2(mcols(pair)[,3])),1) # mcols(pair)[,2] / mcols(pair)[,3])
  rea <- ifelse(ratio>=1.5, 'amplification', 
                ifelse(ratio<=0.5, 'deletion', 'neutral'))
  values(pair) <- cbind(values(pair), DataFrame(ratio, rea))
  values(pair)[,c(2,3)] <- DataFrame(round(as.matrix(values(pair)[,c(2,3)]),0))
  
  pair
}

getGene <- function(x) {
  x <- unlist(strsplit(x, split=';'))
  unique(gsub('gene_name=', '', x[grep('gene_name', x)]))
}

###############################################################################
# read in and process data
###############################################################################

# read bed file and convert to granges object
cat('Reading in bed file ...\n')
bed <- read.csv2(opt$bed, sep='\t', header=F, stringsAsFactors = F)[,c(1:4)]
bed$id <- paste(bed$V1, bed$V2, sep = '_')
colnames(bed) <- c('chrom', 'start', 'end', 'gene', 'id')
#remove duplicted regions if they exist
if(any(duplicated(bed$id))) {
  bed <- bed[-which(duplicated(bed$id)),]
}
grBed <- as(bed, 'GRanges')

# read in tumor / germline coverage data
cat('Reading in data ...\n')
files <- c(opt$tumor, opt$germline)
dp <- lapply(files, function(x) {
  cat(paste0('   ', x, '\n'))
  #x <- read.csv2(x, sep='\t', stringsAsFactors = F)
  x <- read.csv2(x, sep='\t', stringsAsFactors = F, skip=1, header=F)[,c(1:9)]
  colnames(x) <- c('chrom',	'chromStart'	,'chromEnd'	,'F3'	,'F4',	'F5',	'readCount'	,'meanCoverage'	,'sampleName')
  x <- cbind(ID=paste(x$chrom, x$chromStart, sep='_'), x)
  x
})
dp2 <- lapply(dp, function(x) x[,c(1,8)])

# convert to matrix
cat('Converting data ...\n')
dpM <- dpToM(dp2, files)

# convert to GRanges object and merge Tumor Germline into one object
grTG <- dpMToGR(dpM)

# reduce to bed file regions with bin size == exons
cat('Binning data ...\n')
grBin <- myBin(grBed, grTG)
grBin <- diffloop::addchr(grBin) ##add chr

# nromalize data based on library size
cat('Normalizing data ...\n')
librarySize <- colSums(dpM[,2:dim(dpM)[2]]) / 1000000
normFactor <- max(librarySize) / librarySize

grNorm <- grBin
values(grNorm)[, 4:dim(values(grNorm))[2]] <- DataFrame(
  sweep(as.matrix(values(grNorm)[, 4:dim(values(grNorm))[2]]),
        MARGIN=2,normFactor,`*`))

# adjust GC content
cat('Adjusting GC content ...\n')
grGC <- grNorm
gcadjust <- list()
for (i in 1:(dim(values(grNorm))[2]-3)) {
  gcadjust[[i]] <- adjustGC(grGC, i+3, plot=F)
}
tt <- DataFrame(gcadjust)
colnames(tt) <- colnames(mcols(grTG))[-1]
values(grGC)[4:dim(values(grNorm))[2]] <- tt

# get CNV regions only
#cat('Extracting CNV regions ...\n')
#grCNV <- grGC[grep('CNV',grGC$gene),]

# get population mean for germline samples
cat('Computing germline population average ...\n')
grCNVp <- grGC
values(grCNVp) <- cbind(values(grCNVp)[,1:4], 
                        DataFrame(rowMeans(
                          as.matrix(values(grGC)[5:dim(values(grGC))[2]]))))
names(mcols(grCNVp))[4:5] <- c('TUMOR', 'NORMAL')

# remove low coverage regions
rem <- which(as.numeric(as.matrix(values(grCNVp)['TUMOR'])) < opt$min_cov & 
               as.numeric(as.matrix(values(grCNVp)['NORMAL'])) < opt$min_cov)
if(length(rem)!=0) {
  grCNVp <- grCNVp[-rem,]  
}

# adjust tumor purity
cat('Adjusting tumor purity ...\n')
grP <- grCNVp
grP$TUMOR <- grP$TUMOR / opt$tumor_purity

# compute CNV ratios
cat('Computing CNV ratios ...\n')
res <- cnvPat(grP)

#values(res) <- DataFrame(cbind(as.data.frame(values(res)), Gene=unlist(lapply(res$gene, getGene))))

# writting results
cat(paste0('Writing results to '), opt$out, '\n')
write.csv(as(res[which(res$ratio >= opt$amplification_cut | 
                         res$ratio <= opt$deletion_cut)][,c(1:5)], 
             'data.frame'), 
          opt$out, 
          row.names = F, 
          quote = F)