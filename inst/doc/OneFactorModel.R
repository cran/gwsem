## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
is_CRAN <- !identical(Sys.getenv("NOT_CRAN"), "true")
knitr::opts_chunk$set(eval = !is_CRAN)

## -----------------------------------------------------------------------------
library(gwsem)

## -----------------------------------------------------------------------------
location <- 'https://jpritikin.github.io/gwsem/gwsemOneFacExample'
phenoData <- read.table(file.path(location, "oneFacphenoData.txt"), header=TRUE)
head(phenoData)

## -----------------------------------------------------------------------------
                                                                                     # You must tell GW-SEM:
addFac <- buildOneFac(phenoData = phenoData,                                         # what the data object is (which you read in above)
                     itemNames = c("tobacco", "cannabis", "alcohol"),                # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                    # what covariates that you want to include in the analysis
                     fitfun = "WLS")                                                 # and the fit function that you would like to use (WLS is much faster than ML)


## -----------------------------------------------------------------------------
addFacFit <- mxRun(addFac)
summary(addFacFit)

## -----------------------------------------------------------------------------
library(curl)
curl_download(file.path(location, 'example.pgen'),
              file.path(tempdir(),'example.pgen'))
curl_download(file.path(location, 'example.pvar'),
              file.path(tempdir(),'example.pvar'))

GWAS(model = addFac,                                                                 # what model object you would like to fit
	snpData = file.path(tempdir(), 'example.pgen'),                                                        # that path to the snpData file.
	out=file.path(tempdir(), "latFac.log"),                                                                # the file that you would like to save the full results into
	SNP=1:200)                                                                       # the index of the snps (how many) you would like to fit

## -----------------------------------------------------------------------------
library(gwsem)

phenoData <- read.table(file.path(location, "oneFacphenoData.txt"), header=TRUE)
head(phenoData)

addFac <- buildOneFac(phenoData = phenoData,                                         # what the data object is (which you read in above)
                     itemNames = c("tobacco", "cannabis", "alcohol"),                # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                    # what covariates that you want to include in the analysis
                     fitfun = "WLS")                                                 # and the fit function that you would like to use (WLS is much faster than ML)

GWAS(model = addFac,                                                                 # what model object you would like to fit
	snpData = file.path(tempdir(), 'example.pgen'),                                                    # that path to the snpData file.
	out=file.path(tempdir(), "latFac.log"))                                                            # the file that you would like to save the full results 

## -----------------------------------------------------------------------------
FullResult <- read.delim(file.path(tempdir(), "latFac.log"))

## -----------------------------------------------------------------------------
succinct <- loadResults(path = file.path(tempdir(), "latFac.log"), focus =  "snp_to_F")

## -----------------------------------------------------------------------------
plot(succinct)

