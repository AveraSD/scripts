require(XML)
require(WriteXLS)

path <- '~/averalajolla/fm/'
files <- dir(path, pattern = '.xml')
out_path = '~/Downloads'

var <- NULL
cnv <- NULL
rea <- NULL

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
    temp <- vector(length=11)
    names(temp) <- c('cds-effect', 'depth', 'functional-effect', 'gene', 'percent-reads',
                    'position', 'protein-effect', 'status', 'strand', 'subclonal', 'transcript')
    temp[match(names(x), names(temp))] <- x
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
                   tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin'])
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
                   tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin'])
                   )
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
                    tissue_of_origin=as.vector(xmlData$`variant-report`$.attrs['tissue-of-origin'])
                    )
    rea <- rbind(rea,reasDf)
  }
}

# do some reformatting
var$depth <- as.numeric(var$depth)
var$percent.reads <- as.numeric(var$percent.reads)
var$subclonal <- toupper(var$subclonal)

cnv$copy.number <- as.numeric(cnv$copy.number)
cnv$equivocal <- toupper(cnv$equivocal)
cnv$ratio <- as.numeric(cnv$ratio)

# export to xls
WriteXLS(c('var', 'cnv', "rea"),
         paste(out_path, '/', 'fm_one', '.xls', sep='')
)
