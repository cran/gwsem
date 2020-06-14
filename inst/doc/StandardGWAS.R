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
location <- 'https://jpritikin.github.io/gwsem/gwsemItemExample/'

# Read the phenotype data into R and look at the data
phenoData <- read.table(file.path(location, "itemPhenoData.txt"),
                        header=TRUE)
head(phenoData)

table(phenoData$tobacco)
table(phenoData$cannabis)
table(phenoData$alcohol)

## -----------------------------------------------------------------------------
phenoData$tobacco  <- mxFactor(phenoData$tobacco  , levels = 0:2)
phenoData$cannabis <- mxFactor(phenoData$cannabis , levels = 0:3)
phenoData$alcohol  <- mxFactor(phenoData$alcohol  , levels = 0:4)

## -----------------------------------------------------------------------------
tob <- buildItem(phenoData = phenoData,                              # the phenotypic data object
                     depVar = c("tobacco"),                          # the outcome or dependent variable 
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),    # the necessary covariates
                     fitfun = "WLS")                                 # the fit function (WLS is much faster than ML)

## -----------------------------------------------------------------------------
tobFit <- mxRun(tob)
summary(tobFit)

## -----------------------------------------------------------------------------
library(curl)
curl_download(file.path(location, 'example.pgen'),
              file.path(tempdir(),'example.pgen'))
curl_download(file.path(location, 'example.pvar'),
              file.path(tempdir(),'example.pvar'))

GWAS(model = tob,                       # where the model is
	snpData = file.path(tempdir(),'example.pgen'),       # where the snpData 
	out=file.path(tempdir(),"tob.log"))   # where you want to save the results

## -----------------------------------------------------------------------------
TobFullResult <- read.delim(file.path(tempdir(), "tob.log"))  

## -----------------------------------------------------------------------------
succinct <- loadResults(path = file.path(tempdir(), "tob.log"),
                        focus =  "snp_to_tobacco")

## -----------------------------------------------------------------------------
plot(succinct)

## -----------------------------------------------------------------------------
library(gwsem)

# Read the phenotype data into R and look at the data
phenoData <- read.table(file.path(location, "itemPhenoData.txt"), header=TRUE)
phenoData$tobacco  <- mxFactor(phenoData$tobacco  , levels = 0:2)
phenoData$cannabis <- mxFactor(phenoData$cannabis , levels = 0:3)
phenoData$alcohol  <- mxFactor(phenoData$alcohol  , levels = 0:4)

tob <- buildItem(phenoData = phenoData, depVar = c("tobacco"), 
                 covariates=c('pc1','pc2','pc3','pc4','pc5'),
                 fitfun = "WLS") 

GWAS(model = tob,
     snpData = file.path(tempdir(), 'example.pgen'),
     out= file.path(tempdir(), "tob.log"))

## -----------------------------------------------------------------------------
multi   <- buildItem(phenoData,                               # the data object
              depVar = c("tobacco", "cannabis", "alcohol"),   # the dependent variables
              covariates=c('pc1','pc2','pc3','pc4','pc5'),    # the covariates
              fitfun = "WLS")                                 # the fit function

## -----------------------------------------------------------------------------
multiFit <- mxRun(multi)
summary(multiFit)

## -----------------------------------------------------------------------------
GWAS(model = multi,
     snpData = file.path(tempdir(), 'example.pgen'),
     out=file.path(tempdir(), "multi.log"))

## -----------------------------------------------------------------------------
succinctTob <- loadResults(path = file.path(tempdir(), "multi.log"),
                           focus =  "snp_to_tobacco")
succinctCan <- loadResults(path = file.path(tempdir(), "multi.log"),
                           focus =  "snp_to_cannabis")
succinctAlc <- loadResults(path = file.path(tempdir(), "multi.log"),
                           focus =  "snp_to_alcohol")

## -----------------------------------------------------------------------------
plot(succinctTob)
plot(succinctCan)
plot(succinctAlc)

