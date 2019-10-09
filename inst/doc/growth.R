## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE------------------------------------------------
library(gwsem)  # load gwsem

## ------------------------------------------------------------------------
manifests<-c("t1","t2","t3","t4","snp")
latents<-c("int","slope")
path <- list(mxPath(from="int",to=c("t1","t2","t3","t4"), free=c(FALSE,FALSE,FALSE,FALSE), value=c(1.0,1.0,1.0,1.0) , arrows=1, label=c("int__t1","int__t2","int__t3","int__t4") ),
                 mxPath(from="slope",to=c("t1","t2","t3","t4"), free=c(FALSE,FALSE,FALSE,FALSE), value=c(0.0,1.0,2.0,3.0) , arrows=1, label=c("slope__t1","slope__t2","slope__t3","slope__t4") ),
                 mxPath(from="one",to=c("int","slope"), free=c(TRUE,TRUE), value=c(0.0,0.0) , arrows=1, label=c("const__int","const__slope") ),
                 mxPath(from="snp",to=c("slope","int"), free=c(TRUE,TRUE), value=c(1.0,0.0) , arrows=1, label=c("snp__slope","snp__int") ),
                 mxPath(from="int",to=c("int"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_int") ),
                 mxPath(from="slope",to=c("slope","int"), free=c(TRUE,TRUE), value=c(1.0,0.1) , arrows=2, label=c("VAR_slope","COV_slope_int") ),
                 mxPath(from="t1",to=c("t1"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_err") ),
                 mxPath(from="t2",to=c("t2"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_err") ),
                 mxPath(from="t3",to=c("t3"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_err") ),
                 mxPath(from="t4",to=c("t4"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_err") ),
                 mxPath(from="snp",to=c("snp"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("snp_res") ),
                 mxPath(from="one",to=c("t1","t2","t3","t4","snp"), free=F, value=0, arrows=1))
model <- mxModel("lgc", 
                 type="RAM",
                 manifestVars = manifests,
                 latentVars = latents,
                 path)

## ---- results='hidden', echo=FALSE---------------------------------------
fam <- read.table(file.path(system.file("extdata", package = "gwsem"), "example.fam"))
N <- nrow(fam)

## ------------------------------------------------------------------------
sim1 <- mxGenerateData(model, N)
head(sim1)

## ------------------------------------------------------------------------
for (ii in 1:5) {
  sim1[[paste0('pc', ii)]] <- rnorm(N)
}

## ------------------------------------------------------------------------
m2 <- mxModel("lgc", type="RAM",
        manifestVars = model$manifestVars,
        latentVars = c(model$latentVars, paste0('pc', 1:5)),
        path,
        mxExpectationRAM(M="M"),
        mxFitFunctionWLS(allContinuousMethod="marginals"),
        mxData(observed=sim1, type="raw", minVariance=0.1, warnNPDacov=FALSE))

## ------------------------------------------------------------------------
m2 <- setupCovariates(m2, paste0('pc', 1:5))

## ------------------------------------------------------------------------
tdir <- tempdir()
dir <- system.file("extdata", package = "gwsem")
snp1 <- GWAS(m2, file.path(dir,"example.pgen"), file.path(tdir, "out.log"), SNP=1)
summary(snp1)

## ------------------------------------------------------------------------
GWAS(m2, file.path(dir,"example.pgen"), file.path(tdir, "out.log"))
got <- loadResults(file.path(tdir, "out.log"), 'snp__slope')

## ------------------------------------------------------------------------
head(got[is.na(got$P),])


## ---- message=FALSE------------------------------------------------------
got$P[is.na(got$P)] <- 1

library(qqman)
manhattan(got)

