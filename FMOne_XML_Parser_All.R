require(XML)
require(WriteXLS)
require(stringr)
#path <- '~/averalajolla/CCB/fm/'
#path <- '/data/fm_xml/'
path <- 'c:/data/fm_xml/'
files <- dir(path, pattern = '.xml')
#out_path = '~/Downloads'
out_path = 'c:/data/'

var <- NULL
cnv <- NULL
rea <- NULL
bio <- NULL
samples <- NULL
nhc <- NULL
for (i in 1:length(files)) {
  dat <- xmlParse(paste0(path,files[i]))
  xmlData <- xmlToList(dat)
  ## variants
  vars <- lapply(xmlData$`variant-report`$`short-variants`, function(x) {
    x <- unlist(x)
    names(x) <- gsub('.attrs.', '', names(x))
    if(any(names(x)=='dna-evidence.sample')) {
      x <- x[-which(names(x)=='dna-evidence.sample')]
    }
    names(x) <- gsub('.attrs.', '', names(x))
    temp <- vector(length=12, mode='character')
    names(temp) <- c('allele-fraction', 'cds-effect', 'depth', 'functional-effect', 'gene', 'percent-reads',
                     'position', 'protein-effect', 'status', 'strand', 'subclonal', 'transcript')
    temp['allele-fraction'] <- x['allele-fraction']
    temp['cds-effect'] <- x['cds-effect']
    temp['depth'] <- x['depth']
    temp['functional-effect'] <- x['functional-effect']
    temp['gene'] <- x['gene']
    temp['percent-reads'] <- x['percent-reads']
    temp['position'] <- x['position']
    temp['protein-effect'] <- x['protein-effect']
    temp['status'] <- x['status']
    temp['strand'] <- x['strand']
    temp['subclonal'] <- x['subclonal']
    temp['transcript'] <- x['transcript']
    temp
  })
  if(length(vars)==0) {
    var <- var
  } else {
    varDf <- data.frame(do.call(rbind, vars),stringsAsFactors = F)
    varDf <- cbind(varDf,
                   Sample=gsub('.xml', '', files[i]),
                   disease=toupper(as.vector(xmlData$`variant-report`$.attrs['disease'])),
                   disease_ontology=as.vector(xmlData$`variant-report`$.attrs['disease-ontology']),
                   gender=as.vector(xmlData$`variant-report`$.attrs['gender']),
                   pathology_diagnosis=as.vector(xmlData$`variant-report`$.attrs['pathology-diagnosis']),
                   percent_tumor_nuclei=as.vector(xmlData$`variant-report`$.attrs['percent-tumor-nuclei']),
                   purity_assessment=as.vector(xmlData$`variant-report`$.attrs['purity-assessment']),
                   tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin']),
                   study=as.vector(xmlData$`variant-report`$.attrs['study']),
                   test_type=as.vector(xmlData$`variant-report`$.attrs['test-type'])
    )
    var <- rbind(var, varDf)
  }
  # cnvs
  cnvs <- lapply(xmlData$`variant-report`$`copy-number-alterations`, function(x) {
    x <- unlist(x)
    names(x) <- gsub('.attrs.', '', names(x))
    if(any(names(x)=='dna-evidence.sample')) {
      x <- x[-which(names(x)=='dna-evidence.sample')]
    }
    x
  })
  if(length(cnvs)==0) {
    cnv <- cnv
  } else {
    cnvDf <- data.frame(do.call(rbind, cnvs), stringsAsFactors = F)
    cnvDf <- cbind(cnvDf,
                   Sample=gsub('.xml', '', files[i]),
                   disease=toupper(as.vector(xmlData$`variant-report`$.attrs['disease'])),
                   disease_ontology=as.vector(xmlData$`variant-report`$.attrs['disease-ontology']),
                   gender=as.vector(xmlData$`variant-report`$.attrs['gender']),
                   pathology_diagnosis=as.vector(xmlData$`variant-report`$.attrs['pathology-diagnosis']),
                   percent_tumor_nuclei=as.vector(xmlData$`variant-report`$.attrs['percent-tumor-nuclei']),
                   purity_assessment=as.vector(xmlData$`variant-report`$.attrs['purity-assessment']),
                   tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin']),
                   study=as.vector(xmlData$`variant-report`$.attrs['study']),
                   test_type=as.vector(xmlData$`variant-report`$.attrs['test-type'])
    )
    if(!any(colnames(cnvDf)=='copy.number') | any(colnames(cnvDf)=='copy-number')) {
      cnvDf <- cbind('copy.number'=NA, cnvDf)
    }
    cnv <- rbind(cnv, cnvDf)
  }
  ## rearangements
  reas <- lapply(xmlData$`variant-report`$rearrangements, function(x) {
    x <- unlist(x)
    names(x) <- gsub('.attrs.', '', names(x))
    if(any(names(x)=='dna-evidence.sample')) {
      x <- x[-which(names(x)=='dna-evidence.sample')]
    }
    if(any(names(x)=='rna-evidence.sample')) {
      x <- x[-which(names(x)=='rna-evidence.sample')]
    }
    temp <- vector(length=10, mode='character')
    names(temp) <- c('description', 'equivocal', 'in-frame', 'other-gene',
                     'pos1', 'pos2', 'status', 'supporting-read-pairs',
                     'targeted-gene', 'type')
    temp['description'] <- x['description']
    temp['equivocal'] <- x['equivocal']
    temp['in-frame'] <- x['in-frame']
    temp['other-gene'] <- x['other-gene']
    temp['pos1'] <- x['pos1']
    temp['pos2'] <- x['pos2']
    temp['status'] <- x['status']
    temp['supporting-read-pairs'] <- x['supporting-read-pairs']
    temp['targeted-gene'] <- x['targeted-gene']
    temp['type'] <- x['type']
    temp
  })
  if(length(reas)==0) {
    rea <- rea
  } else if(is.null(reas$rearrangement)) {
    rea <- rea
  } else {
    reasDf <- data.frame(do.call(rbind,reas), stringsAsFactors = F)
    reasDf <- cbind(reasDf,
                    Sample=gsub('.xml', '', files[i]),
                    disease=toupper(as.vector(xmlData$`variant-report`$.attrs['disease'])),
                    disease_ontology=as.vector(xmlData$`variant-report`$.attrs['disease-ontology']),
                    gender=as.vector(xmlData$`variant-report`$.attrs['gender']),
                    pathology_diagnosis=as.vector(xmlData$`variant-report`$.attrs['pathology-diagnosis']),
                    percent_tumor_nuclei=as.vector(xmlData$`variant-report`$.attrs['percent-tumor-nuclei']),
                    purity_assessment=as.vector(xmlData$`variant-report`$.attrs['purity-assessment']),
                    tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin']),
                    study=as.vector(xmlData$`variant-report`$.attrs['study']),
                    test_type=as.vector(xmlData$`variant-report`$.attrs['test-type'])
    )
    rea <- rbind(rea,reasDf)
  }
  # biomarker
  if(any(names(xmlData$`variant-report`)=='biomarkers')) {
    if(!is.null(xmlData$`variant-report`$biomarkers)) {
      bioDf <- data.frame(Sample=gsub('.xml', '', files[i]),
                          microsatellite.instability=NA,
                          tumor.mutation.burden.score=NA,
                          tumor.mutation.burden.status=NA,
                          tumor.mutation.burden.unit=NA,
                          disease=toupper(as.vector(xmlData$`variant-report`$.attrs['disease'])),
                          disease_ontology=as.vector(xmlData$`variant-report`$.attrs['disease-ontology']),
                          gender=as.vector(xmlData$`variant-report`$.attrs['gender']),
                          pathology_diagnosis=as.vector(xmlData$`variant-report`$.attrs['pathology-diagnosis']),
                          percent_tumor_nuclei=as.vector(xmlData$`variant-report`$.attrs['percent-tumor-nuclei']),
                          purity_assessment=as.vector(xmlData$`variant-report`$.attrs['purity-assessment']),
                          tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin']),
                          study=as.vector(xmlData$`variant-report`$.attrs['study']),
                          test_type=as.vector(xmlData$`variant-report`$.attrs['test-type']),
                          stringsAsFactors = F
      )
      if(!is.null(xmlData$`variant-report`$biomarkers$`microsatellite-instability`)) {
        bioDf$microsatellite.instability <- xmlData$`variant-report`$biomarkers$`microsatellite-instability`['status']
      }
      if(!is.null(xmlData$`variant-report`$biomarkers$`tumor-mutation-burden`)) {
        bioDf$tumor.mutation.burden.score <- xmlData$`variant-report`$biomarkers$`tumor-mutation-burden`['score']
        bioDf$tumor.mutation.burden.status <- xmlData$`variant-report`$biomarkers$`tumor-mutation-burden`['status']
        bioDf$tumor.mutation.burden.unit <- xmlData$`variant-report`$biomarkers$`tumor-mutation-burden`['unit']
      }
      bio <- rbind(bio,bioDf)
    } else {
      bio <- bio
    }
  }
  #sample info
  if(any(names(xmlData$`variant-report`)=='samples')) {
    if(!is.null(xmlData$`variant-report`$samples)) {
      samplesDf <- data.frame(Sample=gsub('.xml', '', files[i]),
                              baitSet=NA,
                              meanExonDepth=NA,
                              name=NA,
                              nucleicAcidType=NA,
                              disease=toupper(as.vector(xmlData$`variant-report`$.attrs['disease'])),
                              disease_ontology=as.vector(xmlData$`variant-report`$.attrs['disease-ontology']),
                              gender=as.vector(xmlData$`variant-report`$.attrs['gender']),
                              pathology_diagnosis=as.vector(xmlData$`variant-report`$.attrs['pathology-diagnosis']),
                              percent_tumor_nuclei=as.vector(xmlData$`variant-report`$.attrs['percent-tumor-nuclei']),
                              purity_assessment=as.vector(xmlData$`variant-report`$.attrs['purity-assessment']),
                              tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin']),
                              study=as.vector(xmlData$`variant-report`$.attrs['study']),
                              test_type=as.vector(xmlData$`variant-report`$.attrs['test-type']),
                              stringsAsFactors = F
      )
      if(!is.null(xmlData$`variant-report`$samples$sample['bait-set'])) {
        samplesDf$baitSet <- xmlData$`variant-report`$samples$sample['bait-set']
      }
      if(!is.null(xmlData$`variant-report`$samples$sample['mean-exon-depth'])) {
        samplesDf$meanExonDepth <- xmlData$`variant-report`$samples$sample['mean-exon-depth']
      }
      if(!is.null(xmlData$`variant-report`$samples$sample['name'])) {
        samplesDf$name <- xmlData$`variant-report`$samples$sample['name']
      }
      if(!is.null(xmlData$`variant-report`$samples$sample['nucleic-acid-type'])) {
        samplesDf$nucleicAcidType <- xmlData$`variant-report`$samples$sample['nucleic-acid-type']
      }
      samples <- rbind(samples,samplesDf)
    } else {
      samples <- samples
    }
  }
  # non human content
  if(any(names(xmlData$`variant-report`)=='non-human-content')) {
    if(!is.null(xmlData$`variant-report`$`non-human-content`)) {
      t <- unlist(xmlData$`variant-report`$`non-human-content`)
      nhcDf <- data.frame(Sample=gsub('.xml', '', files[i]),
                          nonHumanDNAEvidenceSample = NA,
                          nonHumananOrganism = NA,
                          nonHumanReadsPerMillion = NA,
                          nonHumanStatus = NA,
                          disease=toupper(as.vector(xmlData$`variant-report`$.attrs['disease'])),
                          disease_ontology=as.vector(xmlData$`variant-report`$.attrs['disease-ontology']),
                          gender=as.vector(xmlData$`variant-report`$.attrs['gender']),
                          pathology_diagnosis=as.vector(xmlData$`variant-report`$.attrs['pathology-diagnosis']),
                          percent_tumor_nuclei=as.vector(xmlData$`variant-report`$.attrs['percent-tumor-nuclei']),
                          purity_assessment=as.vector(xmlData$`variant-report`$.attrs['purity-assessment']),
                          tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin']),
                          study=as.vector(xmlData$`variant-report`$.attrs['study']),
                          test_type=as.vector(xmlData$`variant-report`$.attrs['test-type']),
                          stringsAsFactors = F
      )
      nhcDf <- nhcDf[rep(seq_len(nrow(nhcDf)), each=length(t)/4),]
      nhcDf$nonHumanDNAEvidenceSample <- t[names(t)=='non-human.dna-evidence.sample']
      nhcDf$nonHumananOrganism <- t[names(t)=='non-human..attrs.organism']
      nhcDf$nonHumanReadsPerMillion <- t[names(t)=='non-human..attrs.reads-per-million']
      nhcDf$nonHumanStatus <- t[names(t)=='non-human..attrs.status']
      nhc <- rbind(nhc,nhcDf)
    } else {
      nhc <- nhc
    }
  }
}
# do some reformatting
var$depth <- as.numeric(var$depth)
var$percent.reads <- as.numeric(var$percent.reads)
var$subclonal <- toupper(var$subclonal)
rownames(var) <- seq(1:dim(var)[1])
cnv$copy.number <- as.numeric(cnv$copy.number)
cnv$equivocal <- toupper(cnv$equivocal)
cnv$ratio <- as.numeric(cnv$ratio)
rownames(cnv) <- seq(1:dim(cnv)[1])
rea$supporting.read.pairs <- as.numeric(rea$supporting.read.pairs)
rownames(rea) <- seq(1:dim(rea)[1])
bio$tumor.mutation.burden.score <- as.numeric(bio$tumor.mutation.burden.score)
rownames(bio) <- seq(1:dim(bio)[1])
# add avera tmb
tmbDat <- read.csv2('c://Users//Tobias//Documents//fmtmb.csv', header=T, stringsAsFactors = F, sep=' ')
bio$tmb <- round(as.numeric(tmbDat$TMB[match(bio$Sample, tmbDat$Sample)]),2)
samples$meanExonDepth <- as.numeric(samples$meanExonDepth)
# add avera tmb
samples$tmb <- round(as.numeric(tmbDat$TMB[match(samples$Sample, tmbDat$Sample)]),2)
nhc$nonHumanReadsPerMillion <- as.numeric(nhc$nonHumanReadsPerMillion)
## add CCD ID
ccd2trf <- read.csv2('c://Users//Tobias//Documents//ccdfun.csv', sep=',', header=T, stringsAsFactors = F)
var$CCDID <- ccd2trf$Research[match(var$Sample, ccd2trf$Barcode)]
cnv$CCDID <- ccd2trf$Research[match(cnv$Sample, ccd2trf$Barcode)]
rea$CCDID <- ccd2trf$Research[match(rea$Sample, ccd2trf$Barcode)]
# querry variants myvariant to annotate CADD & ExAc
library(myvariant)
getHGVS <- function(genes, peffect) {
  pos <- subset(var, gene==genes & protein.effect==peffect)$position[1]
  cds <- subset(var, gene==genes & protein.effect==peffect)$cds.effect[1]
  strand <- subset(var, gene==genes & protein.effect==peffect)$strand[1]
  if(strand=='-') {
    cds <- paste(lapply(unlist(strsplit(cds, '')), function(x) {
      r <- x
      if(x=='T') {r='A'}
      if(x=='A') {r='T'}
      if(x=='G') {r='C'}
      if(x=='C') {r='G'}
      r
    }), collapse = ''
    )
  }
  hgvs <- paste0(gsub(':', ':g.', pos), gsub('^[0-9_+-]*', '', cds))
  hgvs
}
hgvs <- apply(var, 1, function(x) {
  getHGVS(x['gene'], x['protein.effect'])
})
#hgvs <- paste0(gsub(':', ':g.', var$position), gsub('^[0-9_+-]*', '', var$cds.effect))
var <- cbind(var, hgvs=hgvs)
head(var)
hgvsUq <- unique(hgvs)
x <- getVariants(hgvsUq, fields = c('cadd.phred',
                                    'exac_nontcga.af',
                                    'snpeff.ann',
                                    'dbsnp.rsid',
                                    'wellderly.adviser_score',
                                    'emv.egl_classification',
                                    'mutdb.mutpred_score',
                                    'docm.domain',
                                    'cosmic.cosmic_id',
                                    'cosmic.mut_freq'#,
                                    #'snpedia.text'
                                    #'cgi',
                                    #'civic'
))
res <- x[-which(x@listData$notfound),]
#res$dbnsfp.mutationtaster.pred <- unlist(lapply(res$dbnsfp.mutationtaster.pred, paste, collapse='/'))
#res$clinvar.rcv <- grepl('pathogenic', unlist(lapply(lapply(res$clinvar.rcv, function(x) x['clinical_significance']), paste, collapse=',')))
res$snpeff.ann <- unlist(lapply(res$snpeff.ann, function(x) {
  if(!is.null(dim(x))) {
    f <- gsub('\\..', '', x$feature_id)
    return(x$putative_impact[which(f %in% var$transcript)][1])
  } else {
    return(x$putative_impact)
  }
}))
res <- data.frame(res)
colnames(res)[c(5,6,8,10,12,13,14,16,18,20)] <- c('COSMIC70', 'COSMIC70_MUTFREQ', 'EFFECT', 'CADD', 'RSID', 'EXAC', 'DOMAIN', 'EMV', 'MUTPRED','WELLDERLY_SCORE')
df <- cbind(var,
            RSID=res$RSID[match(var$hgvs, res$query)],
            EFFECT=res$EFFECT[match(var$hgvs, res$query)],
            CADD=res$CADD[match(var$hgvs, res$query)],
            EXAC=res$EXAC[match(var$hgvs, res$query)],
            COSMIC70=res$COSMIC70[match(var$hgvs, res$query)],
            COSMIC70_MUTFREQ=res$COSMIC70_MUTFREQ[match(var$hgvs, res$query)],
            EMV=res$EMV[match(var$hgvs, res$query)],
            MUTPRED=res$MUTPRED[match(var$hgvs, res$query)],
            WELLDERLY_SCORE=res$WELLDERLY_SCORE[match(var$hgvs, res$query)],
            DOMAIN=res$DOMAIN[match(var$hgvs, res$query)]
)
var <- df

####
###
alldb <- unique(sapply(strsplit(colnames(y), '\\.'), '[[', 1))



###

## annotate with gavin
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1141-7
################################
## GAVIN Applied Rule Guide r0.3
################################
##
## Variant can be interpreted by following the columns from left to right.
## This classification scheme was implemented in Java, and used to
## benchmark GAVIN in the paper (Van der Velde et al., using r0.2), see:
## https://github.com/molgenis/molgenis and https://github.com/molgenis/gavin
##
## Genome-wide rules are used if the gene-specific rules fail to classify.
## These rules are applied as follows:
## 1) If impact equals MODIFIER -> benign,
## 2) if MAF greater than 0.003456145 -> benign,
## 3) if CADD greater than 15 -> pathogenic,
## 4) if CADD less than 15 -> benign.
##
## Explanation of the gene calibration categories:
## C1 = CADD scores highly significantly predictive for pathogenicity (pval < 0.01).
## C2 = CADD scores significantly predictive for pathogenicity (pval < 0.05).
## C3 = CADD scores may be predictive for pathogenicity (pval > 0.05 but with few samples).
## C4 = CADD scores less predictive for pathogenicity (pval > 0.05 with enough samples).
## C5 = CADD scores less predictive for pathogenicity (population CADD > pathogenic CADD).
## I1 = HIGH impact unique for, thus predictive, for pathogenic variants.
## I2 = MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.
## I3 = LOW or MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.
## T1 = Too few exac variants after filtering with pathogenic 95th percentile MAF.
## T2 = Too few exac variants after filtering with impact distribution.
## N1 = Too few ClinVar variants for calibration at this time.
## N2 = Too few ExAC variants found for calibration.
##
## For C1 and C2, CADD score thresholds are means of stratified benign and pathogenic variants.
## For C3, C4 and C5, CADD score thresholds are 95th sensitivity/specificity percentiles of stratified benign and pathogenic variants.
##
#gavinRules <- read.csv2('~/Documents/GAVIN_ruleguide_r0.3.tsv', skip=33, sep='\t', dec='.', stringsAsFactors = F, na.strings = 'n/a')
gavinRules <- read.csv2('C:/Users/Tobias/Documents/GAVIN_ruleguide_r0.3.tsv', skip=33, sep='\t', dec='.', stringsAsFactors = F, na.strings = 'n/a')
gavin <- function(gene, cadd, exac, effect) {
  pred <- NA
  cat <- NA
  g <- 0
  ## calibrations
  cal <- list(C1 = 'CADD scores highly significantly predictive for pathogenicity (pval < 0.01).',
              C2 = 'CADD scores significantly predictive for pathogenicity (pval < 0.05).',
              C3 = 'CADD scores may be predictive for pathogenicity (pval > 0.05 but with few samples).',
              C4 = 'CADD scores less predictive for pathogenicity (pval > 0.05 with enough samples).',
              C5 = 'CADD scores less predictive for pathogenicity (population CADD > pathogenic CADD).',
              I1 = 'HIGH impact unique for, thus predictive, for pathogenic variants.',
              I2 = 'MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.',
              I3 = 'LOW or MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.',
              T1 = 'Too few exac variants after filtering with pathogenic 95th percentile MAF.',
              T2 = 'Too few exac variants after filtering with impact distribution.',
              N1 = 'Too few ClinVar variants for calibration at this time.',
              N2 = 'Too few ExAC variants found for calibration.'
  )
  if(gene %in% gavinRules$Gene) {
    r <- subset(gavinRules, Gene==gene)
    if (!is.na(cadd) & !is.na(r$PathogenicIfCADDScoreGreaterThan) & cadd >= r$PathogenicIfCADDScoreGreaterThan) {
      pred <- 'pathogenic'
      cat <- as.vector(unlist(cal[r$CalibrationCategory]))
      g <- 1
    } else if (!is.na(cadd) & !is.na(r$PathogenicIfCADDScoreGreaterThan) & cadd < r$BenignIfCADDScoreLessThan) {
      pred <- 'benign'
      cat <- as.vector(unlist(cal[r$CalibrationCategory]))
      g <- 1
    } else if (!is.na(exac) & !is.na(r$BenignIfMAFGreaterThan) & exac >= r$BenignIfMAFGreaterThan) {
      pred <- 'benign'
      cat <- as.vector(unlist(cal[r$CalibrationCategory]))
      g <- 1
    } else if (!is.na(effect) & !is.na(r$PathogenicIfImpactEquals) & effect==r$PathogenicIfImpactEquals) {
      pred <- 'pathogenic'
      cat <- as.vector(unlist(cal[r$CalibrationCategory]))
      g <- 1
    }
  }
  # apply global fallback
  if (g==0 & !is.na(effect) & effect=='MODIFIER') {
    pred <- 'benign'
    cat <- 'Genome wide rule'
  } else if (g==0 &  !is.na(exac) & exac >= 0.003456145) {
    pred <- 'benign'
    cat <- 'Genome wide rule'
  } else if (g==0 &  !is.na(cadd) & cadd >= 15) {
    pred <- 'pathogenic'
    cat <- 'Genome wide rule'
  } else if (g==0 & !is.na(cadd) & cadd < 15) {
    pred <- 'benign'
    cat <- 'Genome wide rule'
  } else if (g==0) {
    pred <- NA
    cat <- NA
  }
  return(list=c(pred,cat))
}
g <- t(apply(var, 1, function(x) {
  gavin(x['gene'], x['CADD'], x['EXAC'], x['EFFECT'])
}))
colnames(g) <- c('GARVIN_PRED', 'GARVIN_CAT')
var <- cbind(var, g)

## ACMG
## use data from pathscan
db <- read.csv2('C:\\data\\db\\pathscan\\complete_db.txt',
                stringsAsFactors = F,
                sep='\t')

evalACMG <- function(fmDatR) {
  chr <- as.numeric(gsub('chr', '', unlist(strsplit(fmDatR$position, ':'))[1]))
  pos <- as.numeric(unlist(strsplit(fmDatR$position, ':'))[2])
  ref <- str_sub(unlist(strsplit(as.vector(fmDatR$hgvs), '>'))[1], -1)
  alt <- unlist(strsplit(as.vector(fmDatR$hgvs), '>'))[2]

  r <- subset(db, chromosome==chr & position==pos & ref==ref & alt==alt)
  if(dim(r)[1]!=0) {
    r[1,7:9]
  } else {
    data.frame(Gene_reviews=NA, OMIM=NA, clinical_sig_code=NA)
  }

}

ac <- NULL
for(i in 1:dim(var)[1]) {
  ac <- rbind(ac, evalACMG(var[i,]))
}
var <- cbind(var, ac)

# variant reocurrance
var$N <- as.numeric(table(var$hgvs)[match(var$hgvs, names(table(var$hgvs)))])
var$freq <- var$N / length(samples$Sample)
cnv$N <- as.numeric(table(cnv$gene)[match(cnv$gene, names(table(cnv$gene)))])
cnv$freq <- cnv$N / length(samples$Sample)

varSum <- dplyr::arrange(var[match(unique(var[var$freq >= 0.005,]$hgvs), var$hgvs),c(5,7,8,9,12,24,25,26,27,28,30,31,33,34)], desc(N))
x1 <- subset(varSum, status=='known')
x2 <- subset(varSum, status=='unknown')

WriteXLS(c('x1', 'x2'), 'C:/data/topVars.xls')

#####
# save data
# export to xls
WriteXLS(c('var', 'cnv', "rea", 'bio', 'samples', 'nhc'),
         paste(out_path, '/', 'fm_03012017', '.xls', sep='')
)
# save data as R object
save(list = c('var','cnv','rea', 'bio', 'samples', 'nhc'), file='C:/data/fm_03012017.Rdata')
