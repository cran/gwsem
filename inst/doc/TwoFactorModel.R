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
location <- 'https://jpritikin.github.io/gwsem/gwsemTwoFacExample'
TwoDat <- read.table(file.path(location, "phenoTwoData.txt"), header=TRUE)
head(TwoDat)

## -----------------------------------------------------------------------------

twoFac <- buildTwoFac(phenoData = TwoDat,                                         # what the data object is (which you read in above)
                     F1itemNames = c("A1", "A2", "A3"),                # what the items of the latent factor are
                     F2itemNames = c("B1", "B2", "B3"),                # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                    # what covariates that you want to include in the analysis
                     fitfun = "WLS",
					 exogenous = T)                                                 # and the fit function that you would like to use (WLS is much faster than ML)


## -----------------------------------------------------------------------------
twoFacFit <- mxRun(twoFac)
summary(twoFacFit)

## -----------------------------------------------------------------------------
library(curl)
curl_download(file.path(location, 'example.pgen'),
              file.path(tempdir(),'example.pgen'))
curl_download(file.path(location, 'example.pvar'),
              file.path(tempdir(),'example.pvar'))

GWAS(model = twoFac,                                                                 # what model object you would like to fit
	snpData = file.path(tempdir(),'example.pgen'),                                                        # that path to the snpData file.
	out=file.path(tempdir(), "twoFac.log"))                                                                # the file that you would like to save the full results into
	                                                                       # the index of the snps (how many) you would like to fit


## -----------------------------------------------------------------------------
library(gwsem)

TwoDat <- read.table(file.path(location, "phenoTwoData.txt"), header=TRUE)
head(TwoDat)

twoFac <- buildTwoFac(phenoData = TwoDat,                              # data object
                     F1itemNames = c("A1", "A2", "A3"),                # items of the first latent factor
                     F2itemNames = c("B1", "B2", "B3"),                # items of the second latent factor 
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),      # covariates 
                     fitfun = "WLS",
		     exogenous = T)                                                 


GWAS(model = twoFac,
     snpData = file.path(tempdir(), 'example.pgen'),
     out=file.path(tempdir(), "twoFac.log"))


## -----------------------------------------------------------------------------
TwoResult <- read.delim(file.path(tempdir(), "twoFac.log"))

## -----------------------------------------------------------------------------
succinct1 <- loadResults(path = file.path(tempdir(), "twoFac.log"), focus =  "snp_to_F1")
succinct2 <- loadResults(path = file.path(tempdir(), "twoFac.log"), focus =  "snp_to_F2")

## -----------------------------------------------------------------------------
plot(succinct1)
plot(succinct2)

## -----------------------------------------------------------------------------
head(succinct1[order(succinct1$Z, decreasing = T),])

head(succinct2[order(succinct2$Z, decreasing = T),])


## -----------------------------------------------------------------------------
facCov <- loadResults(path = file.path(tempdir(), "twoFac.log"), focus =  "facCov", .retainSE = T)

subset(facCov, SNP == "snp141")

subset(facCov, SNP == "snp50")

subset(facCov, SNP == "snp719")

head(facCov)


