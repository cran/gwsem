## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
is_CRAN <- !identical(Sys.getenv("NOT_CRAN"), "true")
knitr::opts_chunk$set(eval = !is_CRAN)

## -----------------------------------------------------------------------------
library(gwsem)

location <- 'https://jpritikin.github.io/gwsem/gwsemGxEexample'

# Read the phenotype data into R and look at the data
gxeData <- read.table(file.path(location, "gxeData.txt"), header=TRUE)

## -----------------------------------------------------------------------------
gxe <- buildItem(phenoData = gxeData, depVar  = "phe",
                  covariates = c("mod", "snp_mod" ,"pc1", "pc2", "pc3", "pc4","pc5"),
                  fitfun = "WLS", exogenous = T, gxe = "mod")

## -----------------------------------------------------------------------------
gxeFit <- mxRun(gxe)
summary(gxeFit)

## -----------------------------------------------------------------------------
library(curl)
curl_download(file.path(location, 'example.pgen'),
              file.path(tempdir(),'example.pgen'))
curl_download(file.path(location, 'example.pvar'),
              file.path(tempdir(),'example.pvar'))

GWAS(model = gxe,                            # the model object
	snpData = file.path(tempdir(), 'example.pgen'),            # the path to the snpData
	out=file.path(tempdir(), "gxe.log"))                       # the results file name


## -----------------------------------------------------------------------------
gxeResult <- read.delim(file.path(tempdir(), "gxe.log"))                            

## -----------------------------------------------------------------------------
succinctCond <- loadResults(path = file.path(tempdir(), "gxe.log"), focus =  "snp_to_phe")
succinctInt  <- loadResults(path = file.path(tempdir(), "gxe.log"), focus =  "snp_mod_to_phe")

head(succinctCond[order(succinctCond$Z, decreasing = T),])
head(succinctInt[order(succinctInt$Z, decreasing = T),])


## -----------------------------------------------------------------------------
margLow  <- loadResults(path = file.path(tempdir(), "gxe.log"), focus =  "snp_mod_to_phe", moderatorLevel= -2)
margHigh <- loadResults(path = file.path(tempdir(), "gxe.log"), focus =  "snp_mod_to_phe", moderatorLevel=  2)


## -----------------------------------------------------------------------------
# Manhattan Plots for directly estimated parameters
plot(succinctCond)  # To plot p-values of the conditional effect
plot(succinctInt)   # To plot the interaction coefficient

# Manhattan Plots for marginal effects
plot(margLow)    # To plot p-values for a low level (-2SD) of the moderator
plot(margHigh)   # To plot p-values for a high level (+2SD) of the moderator


