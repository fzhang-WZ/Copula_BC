# You need to run these every time you start R.
library(JuliaCall)
library(zoo)
library(quantreg)
library(RColorBrewer)
julia_setup(JULIA_HOME = '/Applications/Julia-1.11.app/Contents/Resources/julia/bin')
julia_library("BulkAndTails")

## loading the data
library(quantreg)
library(RColorBrewer)


#########################################################################################################
#########################################################################################################
#########################################################################################################
data = read.csv("~/Desktop/argonne07.csv")
data$obs_Spd60m[1] <- 5


###############################################################################
####################### box-cox spline #############################
###############################################################################
library(pbs)

diurnalmatrix = pbs::pbs(hour %% 24,df=3, Boundary.knots = c(0,24))
seasonalmatrix = pbs::pbs(day %% 365.25,df=4, Boundary.knots = c(0,365.25))

bigbasismatrix = data.frame(cbind(diurnalmatrix, seasonalmatrix))
phimatrix = cbind(rep(1.0, length(hour)), diurnalmatrix, seasonalmatrix)
taumatrix = cbind(rep(1.0, length(hour)), diurnalmatrix, seasonalmatrix)

## example initial guess

regcoef.bc.spline.obs10 = lm(bc.obs10 ~ seasonalmatrix+diurnalmatrix)$coef
regcoef.bc.spline.obs60 = lm(bc.obs60 ~ seasonalmatrix+diurnalmatrix)$coef
regcoef.bc.spline.era10 = lm(bc.era10 ~ seasonalmatrix+diurnalmatrix)$coef
regcoef.bc.spline.era60 = lm(bc.era60 ~ seasonalmatrix+diurnalmatrix)$coef

phi0init.bc.spline.obs10 = regcoef.bc.spline.obs10
phi1init.bc.spline.obs10 = regcoef.bc.spline.obs10
phi0init.bc.spline.era10 = regcoef.bc.spline.era10
phi1init.bc.spline.era10 = regcoef.bc.spline.era10
phi0init.bc.spline.obs60 = regcoef.bc.spline.obs60
phi1init.bc.spline.obs60 = regcoef.bc.spline.obs60
phi0init.bc.spline.era60 = regcoef.bc.spline.era60
phi1init.bc.spline.era60 = regcoef.bc.spline.era60

kappa0init = 0.0
kappa1init = 0.0
nuinit = 1.0
tau0init = rep(0,3)
tau1init = rep(0,3)

initguess.bc.spline.obs10 = c(kappa0init, tau0init, phi0init.bc.spline.obs10, kappa1init, tau1init, phi1init.bc.spline.obs10, nuinit)
initguess.bc.spline.era10 = c(kappa0init, tau0init, phi0init.bc.spline.era10, kappa1init, tau1init, phi1init.bc.spline.era10, nuinit)
initguess.bc.spline.obs60 = c(kappa0init, tau0init, phi0init.bc.spline.obs60, kappa1init, tau1init, phi1init.bc.spline.obs60, nuinit)
initguess.bc.spline.era60 = c(kappa0init, tau0init, phi0init.bc.spline.era60, kappa1init, tau1init, phi1init.bc.spline.era60, nuinit)


## fit
# note this takes 61 iterations and converges after about 5 mins on my laptop
# this will return an answer even if the software hasn't converged, see status below

bc.spline.obs10 = julia_call("fitbats_covariates", bc.obs10, taumatrix, phimatrix, initguess.bc.spline.obs10, print_level=5L)
bc.spline.era10 = julia_call("fitbats_covariates", bc.era10, taumatrix, phimatrix, initguess.bc.spline.era10, print_level=5L)
bc.spline.obs60 = julia_call("fitbats_covariates", bc.obs60, taumatrix, phimatrix, initguess.bc.spline.obs60, print_level=5L)
bc.spline.era60 = julia_call("fitbats_covariates", bc.era60, taumatrix, phimatrix, initguess.bc.spline.era60, print_level=5L)

## plotting at various seasons

seasonlist = list("DJF" = c(331:365,1:60), "MAM" = 61:150, "JJA" = 151:240, "SON" = 241:330)
whichdays = sapply(seasonlist, function(x){floor(median(x))}) #picking the median for each season as representative

quantvals = c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
quantilemat_bc_spline_obs10 = matrix(NA, 24, length(quantvals))
quantilemat_qr = matrix(NA, 24, length(quantvals))

library(quantreg)
library(RColorBrewer)
for(k in 1:length(seasonlist)) {
  
  thisday = whichdays[k]
  season = names(seasonlist)[k]
  
  smalldiurnalmatrix = pbs::pbs(1:24 %% 24,df=3, Boundary.knots = c(0,24))
  smallseasonalmatrix = matrix(rep(seasonalmatrix[which(day == thisday)[1],],24),nrow=24,byrow = TRUE)
  smallphimatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  smalltaumatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  
  # extracting the MLE coefficients
  phi0vec = c(smallphimatrix %*% bc.spline.obs10$phi0)
  phi1vec = c(smallphimatrix %*% bc.spline.obs10$phi1)
  tau0vec = c(exp(smalltaumatrix %*% bc.spline.obs10$tau0))
  tau1vec = c(exp(smalltaumatrix %*% bc.spline.obs10$tau1))
  nu = bc.spline.obs10$nu
  kappa0 = bc.spline.obs10$kappa0
  kappa1 = bc.spline.obs10$kappa1
  
  # getting some empirical hourly max/med/min for the season
  
  minvals = rep(NA,24)
  maxvals = rep(NA,24)
  medvals = rep(NA,24)
  
  for (i in 1:24){
    subsdat = bc.obs10[which(((day %in% seasonlist[[k]])) & hour == i)]
    maxvals[i] = max(subsdat)
    minvals[i] = min(subsdat)
    medvals[i] = median(subsdat)
  }
  
  # quantiles from the MLE at 24 hours at the selected day
  
  for (j in 1:length(quantvals)){
    thisquantile = quantvals[j]
    
    for(i in 1:24){
      quantilemat_bc_spline_obs10[i,j] = julia_call("batsquantile", thisquantile, c(kappa0,tau0vec[i],phi0vec[i],kappa1,tau1vec[i],phi1vec[i],nu))
    }
  }
  
  pdf(paste0("~/Downloads/","pix/quantiles",season,"_bc_spline.pdf"))
  
  name = paste0("OBS (", season, ")") #which datatype was fit?
  
  colz <-
    c("blue",
      "blue",
      "blue",
      "green3",
      "green3",
      "green3",
      "red",
      "red",
      "red")
  matplot(
    quantilemat,
    type = "l",
    col = colz,
    lty = c(1),
    main = name,
    axes = F,
    xlab = "",
    ylab = "Box-Cox Wind Speed",
    ylim = c(ylimlower, ylimupper),
    cex.main = 2.3,
    cex.lab = 1.5,
    lwd = 1.5
  )
  
  axis(
    side = 1,
    at = 1:24,
    las = 2,
    cex.axis = 1.5
  )
  axis(2, cex.axis = 1.5)
  
  title(xlab = "Hour", cex.lab=1.5, line=4)
  
  lines(maxvals)
  lines(minvals)
  lines(medvals)
  dev.off()
  
  # densities
  
  mypdf = matrix(NA,length(pdfvals),24)
  
  for (j in 1:length(pdfvals)){
    thisquantile = pdfvals[j]
    
    for(i in 1:24){
      mypdf[j,i] = julia_call("batspdf", thisquantile, c(kappa0,tau0vec[i],phi0vec[i],kappa1,tau1vec[i],phi1vec[i],nu))
    }
  }
  
  pdf(paste0("~/Downloads/","pix/density",season,"_bc_spline.pdf"))
  # colz <- c("blue", "green3", "red", "black")
  # # Create a color ramp function using the list of colors
  # rampFunc <- colorRampPalette(colz)
  # interpColors <- rampFunc(50)
  
  colz <- colorRampPalette(brewer.pal(8, "Greys")[3:8])(24)
  
  matplot(
    pdfvals,
    mypdf,
    type = "l",
    col = colz,
    lty = 1,
    main = name,
    yaxs = "i",
    ylim = c(0, max(mypdf)),
    xlab = "Box-Cox Wind Speed",
    ylab = "Density",
    cex.axis = 1.5,
    cex.main = 2.3,
    cex.lab = 1.5,
    lwd = 1.5
  )
  # matplot(
  # kdematx,
  # kdematy,
  # type = "l",
  # col = colz,
  # lty = "33",
  # add = TRUE,
  # lwd = 1.5
  # )
  
  dev.off()
}


quantilemat = matrix(NA, 24, length(quantvals))

library(quantreg)
library(RColorBrewer)
for(k in 1:length(seasonlist)) {
  
  thisday = whichdays[k]
  season = names(seasonlist)[k]
  
  smalldiurnalmatrix = pbs::pbs(1:24 %% 24,df=3, Boundary.knots = c(0,24))
  smallseasonalmatrix = matrix(rep(seasonalmatrix[which(day == thisday)[1],],24),nrow=24,byrow = TRUE)
  smallphimatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  smalltaumatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  
  
  # extracting the MLE coefficients
  phi0vec = c(smallphimatrix %*% mle.era10.spline$phi0)
  phi1vec = c(smallphimatrix %*% mle.era10.spline$phi1)
  tau0vec = c(exp(smalltaumatrix %*% mle.era10.spline$tau0))
  tau1vec = c(exp(smalltaumatrix %*% mle.era10.spline$tau1))
  nu = mle.era10.spline$nu
  kappa0 = mle.era10.spline$kappa0
  kappa1 = mle.era10.spline$kappa1
  
  # getting some empirical hourly max/med/min for the season
  
  minvals = rep(NA,24)
  maxvals = rep(NA,24)
  medvals = rep(NA,24)
  
  for (i in 1:24){
    subsdat = data$era10[which(((day %in% seasonlist[[k]])) & hour == i)]
    maxvals[i] = max(subsdat)
    minvals[i] = min(subsdat)
    medvals[i] = median(subsdat)
  }
  
  # quantiles from the MLE at 24 hours at the selected day
  
  for (j in 1:length(quantvals)){
    thisquantile = quantvals[j]
    
    for(i in 1:24){
      quantilemat[i,j] = ((julia_call("batsquantile", thisquantile, c(kappa0,tau0vec[i],phi0vec[i],kappa1,tau1vec[i],phi1vec[i],nu)))*0.3+1)^(10/3)
    }
  }
  
  pdf(paste0("~/Downloads/","pix/quantiles",season,"_bc_ws_spline_era10.pdf"))
  
  name = paste0("ERA10 (", season, ")") #which datatype was fit?
  
  colz <-
    c("blue",
      "blue",
      "blue",
      "green3",
      "green3",
      "green3",
      "red",
      "red",
      "red")
  matplot(
    quantilemat,
    type = "l",
    col = colz,
    lty = c(1),
    main = name,
    axes = F,
    xlab = "",
    ylab = "Wind Speed",
    ylim = c((ylimlower10*0.3+1)^(10/3), (ylimupper10*0.3+1)^(10/3)),
    cex.main = 2.3,
    cex.lab = 1.5,
    lwd = 1.5
  )
  
  axis(
    side = 1,
    at = 1:24,
    las = 2,
    cex.axis = 1.5
  )
  axis(2, cex.axis = 1.5)
  
  title(xlab = "Hour", cex.lab=1.5, line=4)
  
  lines(maxvals)
  lines(minvals)
  lines(medvals)
  dev.off()
}


## example initial guess
regcoef.era60 = lm(bc.era60 ~ diurnalmatrix+seasonalmatrix)$coef

phi0init.era60 = regcoef.era60
phi1init.era60 = regcoef.era60
kappa0init = 0.0
kappa1init = 0.0
nuinit = 1.0
tau0init = rep(0,7+1)
tau1init = rep(0,7+1)

initguess.era60 = c(kappa0init, tau0init, phi0init.era60, kappa1init, tau1init, phi1init.era60, nuinit)

## fit
# this will return an answer even if the software hasn't converged, see status below

mle.era60.spline = julia_call("fitbats_covariates", bc.era60, taumatrix, phimatrix, initguess.era60, print_level=5L)


quantilemat = matrix(NA, 24, length(quantvals))

library(quantreg)
library(RColorBrewer)
for(k in 1:length(seasonlist)) {
  
  thisday = whichdays[k]
  season = names(seasonlist)[k]
  
  smalldiurnalmatrix = pbs::pbs(1:24 %% 24,df=3, Boundary.knots = c(0,24))
  smallseasonalmatrix = matrix(rep(seasonalmatrix[which(day == thisday)[1],],24),nrow=24,byrow = TRUE)
  smallphimatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  smalltaumatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  
  
  # extracting the MLE coefficients
  phi0vec = c(smallphimatrix %*% mle.era60.spline$phi0)
  phi1vec = c(smallphimatrix %*% mle.era60.spline$phi1)
  tau0vec = c(exp(smalltaumatrix %*% mle.era60.spline$tau0))
  tau1vec = c(exp(smalltaumatrix %*% mle.era60.spline$tau1))
  nu = mle.era10.spline$nu
  kappa0 = mle.era60.spline$kappa0
  kappa1 = mle.era60.spline$kappa1
  
  # getting some empirical hourly max/med/min for the season
  
  minvals = rep(NA,24)
  maxvals = rep(NA,24)
  medvals = rep(NA,24)
  
  for (i in 1:24){
    subsdat = na.approx(data$era_Spd60m)[which(((day %in% seasonlist[[k]])) & hour == i)]
    maxvals[i] = max(subsdat)
    minvals[i] = min(subsdat)
    medvals[i] = median(subsdat)
  }
  
  # quantiles from the MLE at 24 hours at the selected day
  
  for (j in 1:length(quantvals)){
    thisquantile = quantvals[j]
    
    for(i in 1:24){
      quantilemat[i,j] =  ((julia_call("batsquantile", thisquantile, c(kappa0,tau0vec[i],phi0vec[i],kappa1,tau1vec[i],phi1vec[i],nu)))*0.6+1)^(10/6)
    }
  }
  
  pdf(paste0("~/Downloads/","pix/quantiles",season,"_bc_ws_spline60_era.pdf"))
  
  name = paste0("ERA60 (", season, ")") #which datatype was fit?
  
  colz <-
    c("blue",
      "blue",
      "blue",
      "green3",
      "green3",
      "green3",
      "red",
      "red",
      "red")
  matplot(
    quantilemat,
    type = "l",
    col = colz,
    lty = c(1),
    main = name,
    axes = F,
    xlab = "",
    ylab = "Wind Speed",
    ylim = c((ylimlower60*0.6+1)^(10/6), 30),
    cex.main = 2.3,
    cex.lab = 1.5,
    lwd = 1.5
  )
  
  axis(
    side = 1,
    at = 1:24,
    las = 2,
    cex.axis = 1.5
  )
  axis(2, cex.axis = 1.5)
  
  title(xlab = "Hour", cex.lab=1.5, line=4)
  
  lines(maxvals)
  lines(minvals)
  lines(medvals)
  dev.off()
}


## example initial guess
regcoef.obs60 = lm(bc.obs60 ~ diurnalmatrix+seasonalmatrix)$coef

phi0init.obs60 = regcoef.obs60
phi1init.obs60 = regcoef.obs60
kappa0init = 0.0
kappa1init = 0.0
nuinit = 1.0
tau0init = rep(0,7+1)
tau1init = rep(0,7+1)

initguess.obs60 = c(kappa0init, tau0init, phi0init.obs60, kappa1init, tau1init, phi1init.obs60, nuinit)

## fit
# this will return an answer even if the software hasn't converged, see status below

mle.obs60.spline = julia_call("fitbats_covariates", bc.obs60, taumatrix, phimatrix, initguess.obs60, print_level=5L)


quantilemat = matrix(NA, 24, length(quantvals))

library(quantreg)
library(RColorBrewer)
for(k in 1:length(seasonlist)) {
  
  thisday = whichdays[k]
  season = names(seasonlist)[k]
  
  smalldiurnalmatrix = pbs::pbs(1:24 %% 24,df=3, Boundary.knots = c(0,24))
  smallseasonalmatrix = matrix(rep(seasonalmatrix[which(day == thisday)[1],],24),nrow=24,byrow = TRUE)
  smallphimatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  smalltaumatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  
  
  # extracting the MLE coefficients
  phi0vec = c(smallphimatrix %*% mle.obs60.spline$phi0)
  phi1vec = c(smallphimatrix %*% mle.obs60.spline$phi1)
  tau0vec = c(exp(smalltaumatrix %*% mle.obs60.spline$tau0))
  tau1vec = c(exp(smalltaumatrix %*% mle.obs60.spline$tau1))
  nu = mle.obs60.spline$nu
  kappa0 = mle.obs60.spline$kappa0
  kappa1 = mle.obs60.spline$kappa1
  
  # getting some empirical hourly max/med/min for the season
  
  minvals = rep(NA,24)
  maxvals = rep(NA,24)
  medvals = rep(NA,24)
  
  for (i in 1:24){
    subsdat = dat[which(((day %in% seasonlist[[k]])) & hour == i)]
    maxvals[i] = max(subsdat)
    minvals[i] = min(subsdat)
    medvals[i] = median(subsdat)
  }
  
  # quantiles from the MLE at 24 hours at the selected day
  
  for (j in 1:length(quantvals)){
    thisquantile = quantvals[j]
    
    for(i in 1:24){
      quantilemat[i,j] = ((julia_call("batsquantile", thisquantile, c(kappa0,tau0vec[i],phi0vec[i],kappa1,tau1vec[i],phi1vec[i],nu)))*0.4+1)^(10/4)
    }
  }
  
  pdf(paste0("~/Downloads/","pix/quantiles",season,"_bc_ws_spline60_obs.pdf"))
  
  name = paste0("OBS60 (", season, ")") #which datatype was fit?
  
  colz <-
    c("blue",
      "blue",
      "blue",
      "green3",
      "green3",
      "green3",
      "red",
      "red",
      "red")
  matplot(
    quantilemat,
    type = "l",
    col = colz,
    lty = c(1),
    main = name,
    axes = F,
    xlab = "",
    ylab = "Wind Speed",
    ylim = c((ylimlower60*0.3+1)^(10/3), 20),
    cex.main = 2.3,
    cex.lab = 1.5,
    lwd = 1.5
  )
  
  axis(
    side = 1,
    at = 1:24,
    las = 2,
    cex.axis = 1.5
  )
  axis(2, cex.axis = 1.5)
  
  title(xlab = "Hour", cex.lab=1.5, line=4)
  
  lines(maxvals)
  lines(minvals)
  lines(medvals)
  dev.off()
}


######################################################################################

ylimupper <- max(dat)
ylimlower <- min(dat.obs)
pdfvals <- seq(min(ylimlower), max(ylimupper), , 1000)

ylimupper10 <- max(bc.obs10,bc.era10)
ylimlower10 <- min(bc.obs10,bc.era10)
pdfvals10 <- seq(min(ylimlower10), max(ylimupper10), , 1000)

ylimupper60 <- max(bc.obs60,bc.era60)
ylimlower60 <- min(bc.obs60,bc.era60)
pdfvals60 <- seq(min(ylimlower60), max(ylimupper60), , 1000)

#nsplines = 8
diurnalmatrix = pbs::pbs(hour %% 24,df=3, Boundary.knots = c(0,24))
seasonalmatrix = pbs::pbs(day %% 365.25,df=4, Boundary.knots = c(0,365.25))

bigbasismatrix = data.frame(cbind(diurnalmatrix, seasonalmatrix))
phimatrix = cbind(rep(1.0, length(hour)), diurnalmatrix, seasonalmatrix)
taumatrix = cbind(rep(1.0, length(hour)), diurnalmatrix, seasonalmatrix)

## example initial guess

regcoef.obs = lm(dat.obs ~ diurnalmatrix+seasonalmatrix)$coef

regcoef.obs10 = lm(bc.obs10 ~ diurnalmatrix+seasonalmatrix)$coef
regcoef.obs60 = lm(bc.obs60 ~ diurnalmatrix+seasonalmatrix)$coef
regcoef.era10 = lm(bc.era10 ~ diurnalmatrix+seasonalmatrix)$coef
regcoef.era60 = lm(bc.era60 ~ diurnalmatrix+seasonalmatrix)$coef

phi0init.era10 = regcoef.era10
phi1init.era10 = regcoef.era10
kappa0init = 0.0
kappa1init = 0.0
nuinit = 1.0
tau0init = rep(0,7+1)
tau1init = rep(0,7+1)

initguess.era10 = c(kappa0init, tau0init, phi0init.era10, kappa1init, tau1init, phi1init.era10, nuinit)

## fit
# this will return an answer even if the software hasn't converged, see status below

mle.era10.spline = julia_call("fitbats_covariates", bc.era10, taumatrix, phimatrix, initguess.era10, print_level=5L)

## plotting at various seasons

seasonlist = list("DJF" = c(331:365,1:60), "MAM" = 61:150, "JJA" = 151:240, "SON" = 241:330)
whichdays = sapply(seasonlist, function(x){floor(median(x))}) #picking the median for each season as representative

quantvals = c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
quantilemat_bc_spline_era10 = matrix(NA, 24, length(quantvals))
quantilemat_qr = matrix(NA, 24, length(quantvals))

library(quantreg)
library(RColorBrewer)
for(k in 1:length(seasonlist)) {
  
  thisday = whichdays[k]
  season = names(seasonlist)[k]
  
  smalldiurnalmatrix = pbs::pbs(1:24 %% 24,df=3, Boundary.knots = c(0,24))
  smallseasonalmatrix = matrix(rep(seasonalmatrix[which(day == thisday)[1],],24),nrow=24,byrow = TRUE)
  smallphimatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  smalltaumatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  
  # extracting the MLE coefficients
  phi0vec = c(smallphimatrix %*% mle.era10.spline$phi0)
  phi1vec = c(smallphimatrix %*% mle.era10.spline$phi1)
  tau0vec = c(exp(smalltaumatrix %*% mle.era10.spline$tau0))
  tau1vec = c(exp(smalltaumatrix %*% mle.era10.spline$tau1))
  nu = mle.era10.spline$nu
  kappa0 = mle.era10.spline$kappa0
  kappa1 = mle.era10.spline$kappa1
  
  # getting some empirical hourly max/med/min for the season
  
  minvals = rep(NA,24)
  maxvals = rep(NA,24)
  medvals = rep(NA,24)
  
  for (i in 1:24){
    subsdat = data$era10[which(((day %in% seasonlist[[k]])) & hour == i)]
    maxvals[i] = max(subsdat)
    minvals[i] = min(subsdat)
    medvals[i] = median(subsdat)
  }
  
  # quantiles from the MLE at 24 hours at the selected day
  
  for (j in 1:length(quantvals)){
    thisquantile = quantvals[j]
    
    for(i in 1:24){
      quantilemat_bc_spline_era10[i,j] = julia_call("batsquantile", thisquantile, c(kappa0,tau0vec[i],phi0vec[i],kappa1,tau1vec[i],phi1vec[i],nu))
    }
         qrfit <- rq(bc.era10 ~ . , data=bigbasismatrix, tau = thisquantile)
         bigbasismatrix2 <- cbind(smallbasismatrix)
        qrpred <- predict(qrfit, newdata = data.frame(bigbasismatrix2))
    
        quantilemat_qr[,j] = qrpred
  }
  
  pdf(paste0("~/Downloads/","pix/quantiles",season,"_bc_spline_era10.pdf"))
  
  name = paste0("ERA10 (", season, ")") #which datatype was fit?
  
  colz <-
    c("blue",
      "blue",
      "blue",
      "green3",
      "green3",
      "green3",
      "red",
      "red",
      "red")
  matplot(
    quantilemat,
    type = "l",
    col = colz,
    lty = c(1),
    main = name,
    axes = F,
    xlab = "",
    ylab = "BC Wind Speed",
    ylim = c(ylimlower, ylimupper),
    cex.main = 2.3,
    cex.lab = 1.5,
    lwd = 1.5
  )
  
  axis(
    side = 1,
    at = 1:24,
    las = 2,
    cex.axis = 1.5
  )
  axis(2, cex.axis = 1.5)
  
  title(xlab = "Hour", cex.lab=1.5, line=4)
  
  lines(maxvals)
  lines(minvals)
  lines(medvals)
  dev.off()
  
  # densities
  
  mypdf = matrix(NA,length(pdfvals),24)
  
  for (j in 1:length(pdfvals)){
    thisquantile = pdfvals[j]
    
    for(i in 1:24){
      mypdf[j,i] = julia_call("batspdf", thisquantile, c(kappa0,tau0vec[i],phi0vec[i],kappa1,tau1vec[i],phi1vec[i],nu))
    }
  }
  
  pdf(paste0("~/Downloads/","pix/density",season,"_log_spline.pdf"))
  # colz <- c("blue", "green3", "red", "black")
  # # Create a color ramp function using the list of colors
  # rampFunc <- colorRampPalette(colz)
  # interpColors <- rampFunc(50)
  
  colz <- colorRampPalette(brewer.pal(8, "Greys")[3:8])(24)
  
  matplot(
    pdfvals,
    mypdf,
    type = "l",
    col = colz,
    lty = 1,
    main = name,
    yaxs = "i",
    ylim = c(0, max(mypdf)),
    xlab = "Log Wind Speed",
    ylab = "Density",
    cex.axis = 1.5,
    cex.main = 2.3,
    cex.lab = 1.5,
    lwd = 1.5
  )
  # matplot(
  # kdematx,
  # kdematy,
  # type = "l",
  # col = colz,
  # lty = "33",
  # add = TRUE,
  # lwd = 1.5
  # )
  
  dev.off()
}


quantilemat = matrix(NA, 24, length(quantvals))

library(quantreg)
library(RColorBrewer)
for(k in 1:length(seasonlist)) {
  
  thisday = whichdays[k]
  season = names(seasonlist)[k]
  
  smalldiurnalmatrix = pbs::pbs(1:24 %% 24,df=3, Boundary.knots = c(0,24))
  smallseasonalmatrix = matrix(rep(seasonalmatrix[which(day == thisday)[1],],24),nrow=24,byrow = TRUE)
  smallphimatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  smalltaumatrix = cbind(1.0, smalldiurnalmatrix, smallseasonalmatrix)
  
  # extracting the MLE coefficients
  phi0vec = c(smallphimatrix %*% mle.obs.spline$phi0)
  phi1vec = c(smallphimatrix %*% mle.obs.spline$phi1)
  tau0vec = c(exp(smalltaumatrix %*% mle.obs.spline$tau0))
  tau1vec = c(exp(smalltaumatrix %*% mle.obs.spline$tau1))
  nu = mle.obs.spline$nu
  kappa0 = mle.obs.spline$kappa0
  kappa1 = mle.obs.spline$kappa1
  
  # getting some empirical hourly max/med/min for the season
  
  minvals = rep(NA,24)
  maxvals = rep(NA,24)
  medvals = rep(NA,24)
  
  for (i in 1:24){
    subsdat = dat[which(((day %in% seasonlist[[k]])) & hour == i)]
    maxvals[i] = max(subsdat)
    minvals[i] = min(subsdat)
    medvals[i] = median(subsdat)
  }
  
  # quantiles from the MLE at 24 hours at the selected day
  
  for (j in 1:length(quantvals)){
    thisquantile = quantvals[j]
    
    for(i in 1:24){
      quantilemat[i,j] = exp(julia_call("batsquantile", thisquantile, c(kappa0,tau0vec[i],phi0vec[i],kappa1,tau1vec[i],phi1vec[i],nu)))
    }
  }
  
  pdf(paste0("~/Downloads/","pix/quantiles",season,"_log_ws_spline.pdf"))
  
  name = paste0("OBS (", season, ")") #which datatype was fit?
  
  colz <-
    c("blue",
      "blue",
      "blue",
      "green3",
      "green3",
      "green3",
      "red",
      "red",
      "red")
  matplot(
    quantilemat,
    type = "l",
    col = colz,
    lty = c(1),
    main = name,
    axes = F,
    xlab = "",
    ylab = "Wind Speed",
    ylim = c(exp(ylimlower), exp(ylimupper)),
    cex.main = 2.3,
    cex.lab = 1.5,
    lwd = 1.5
  )
  
  axis(
    side = 1,
    at = 1:24,
    las = 2,
    cex.axis = 1.5
  )
  axis(2, cex.axis = 1.5)
  
  title(xlab = "Hour", cex.lab=1.5, line=4)
  
  lines(maxvals)
  lines(minvals)
  lines(medvals)
  dev.off()
}

