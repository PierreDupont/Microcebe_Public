########################################################
## ---------------- MANDENA MICROCEBE --------------- ##
## ---- CAPTURE-MARK_RECAPTURE OUTPUT PROCESSING ---- ##
########################################################
rm(list = ls())


## ------ LIBRARIES ------

library(coda)
library(nimble)
library(xtable)
library(data.table)
library(ggplot2)
library(dplyr)
library(magrittr)



## ------ WORKING DIRECTORIES ------

source("./functions/ProcessCodaOutput.R")
source("./functions/PlotViolins.R")



## -----------------------------------------------------------------------------
## ------   1. SET-UP -----

##-- Load model input
load("./data/Microcebe_Mandena_data.RData")

##-- Load model MCMC outputs
load("./data/Microcebe_Mandena_output.RData")

##-- Process output for easier use
results <- ProcessCodaOutput(x = nimOutput) 

##-- Get model features
n.chains <- length(nimOutput)
n.iterations <- dim(nimOutput[[1]])[1]
n.months <- nimConstants$n.months
temp <- nimConstants$temp
month <- nimConstants$month

##-- Choose colors for plots
myCols <- hcl.colors(6)



## ------   2. RJ-MCMC PLOTS ------

##---- Get target covariates
covNames <- c("lphi0.dist.f","lphi0.dist.m","lphi0.prot.f","lphi0.prot.m",
              "temp.dist.f", "temp.dist.m", "temp.prot.f", "temp.prot.m",
              "time.dist.f","time.dist.m", "time.prot.f", "time.prot.m",
              "transloc.dist.f","transloc.dist.m", "transloc.prot.f", "transloc.prot.m")

zRJ.wide <- data.table(cbind(rep(1,200000),rep(1,200000),rep(1,200000),rep(1,200000),
                             res$sims.list$z.temp[ ,1,1:2],res$sims.list$z.temp[ ,2,1:2],
                             res$sims.list$z.time[ ,1,1:2],res$sims.list$z.time[ ,2,1:2],
                             res$sims.list$z.transloc[ ,1,1:2],res$sims.list$z.transloc[ ,2,1:2]))
dimnames(zRJ.wide) <- list(NULL, covNames)

betas.wide <- data.table(cbind(res$sims.list$logit.phi0[ ,1,1:2],res$sims.list$logit.phi0[ ,2,1:2],
                               res$sims.list$beta.temp[ ,1,1:2],res$sims.list$beta.temp[ ,2,1:2],
                               res$sims.list$beta.time[ ,1,1:2],res$sims.list$beta.time[ ,2,1:2],
                               res$sims.list$beta.transloc[ ,1,1:2],res$sims.list$beta.transloc[ ,2,1:2]))
dimnames(betas.wide) <- list(NULL, covNames)


##-- List model combinations
mods <- apply(zRJ.wide, 1, function(x){paste(covNames[x == 1], collapse = "+")})

betas.wide$model <- zRJ.wide$model <- gsub("\\(Intercept\\)\\+", "", mods)
betas.wide$chain <- zRJ.wide$chain <- rep(1:n.chains, each = n.iterations)
betas.wide$iteration <- zRJ.wide$iteration <- rep(1:n.iterations, n.chains)

zRJ.df <- melt(zRJ.wide, id.vars = c("iteration", "chain", "model"))
names(zRJ.df) <- c("iteration", "chain", "model", "variable", "value")

betas.df <-  melt(betas.wide, id.vars = c("iteration", "chain", "model"))
names(betas.df) <-  c("iteration", "chain", "model", "variable", "value")

betas.df$value[zRJ.df$value == 0] <- NA

betas.aggr <- data.table(do.call(rbind, lapply(levels(betas.df$variable), function(x){
  tmp <- betas.df[betas.df$variable == x, ]
  out <- cbind("variable" = x,
               "p.inclusion" = mean(!is.na(tmp$value)))
})))

betas.df <- merge(betas.df, betas.aggr)

included <- zRJ.df$value == 1

betas.df <- betas.df[included, ]
zRJ.df <- zRJ.df[included, ]

betas.df <- betas.df[order(betas.df$variable, betas.df$model, betas.df$chain), ]
zRJ.df <- zRJ.df[order(zRJ.df$variable, zRJ.df$model, zRJ.df$chain), ]

myfun1 <- function(x) 1:length(x)

temp <- betas.df %>% group_by(variable, model, chain) %>%
  summarize(iteration.model = myfun1(value))

betas.df$iteration.model <- temp$iteration.model

aggr <- data.frame(table(betas.df$model) / length(betas.df$model))
names(aggr) <- c("model", "weight")
betas.df <- merge(betas.df, aggr)

##-- Identify parameters by sex
betas.df$sex <- "male"
betas.df$sex[grep(pattern = ".f",x = betas.df$variable)] <- "female"

##-- Identify parameters by protection status
betas.df$status <- "protected"
betas.df$status[grep(pattern = ".dist.",x = betas.df$variable)] <- "fragmented"

##-- Use simpler variable names for plotting
betas.df$variable.simple <- betas.df$variable
betas.df$variable.simple <- gsub('\\b.f\\b','',betas.df$variable.simple)
betas.df$variable.simple <- gsub('\\b.m\\b','',betas.df$variable.simple)
betas.df$variable.simple <- gsub('\\b.dist\\b','',betas.df$variable.simple)
betas.df$variable.simple <- gsub('\\b.prot\\b','',betas.df$variable.simple)
unique(betas.df$variable.simple)



##---- MODEL TALLY
aggr <- aggr[order(aggr$weight, decreasing = TRUE),]
aggr$model <- factor(aggr$model, levels = aggr$model)
reduced_aggr <- filter(aggr, weight >= 0.01)


##---- PLOT MODEL WEIGHTS
ggplot(data = reduced_aggr,
       mapping =  aes(x = model, y = weight, alpha = weight)) +
  geom_col(fill = "magenta") +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 1,
    hjust = 1
  )) + ylab("Weight") + xlab("Models")


##---- PLOT COEFFICIENT ESTIMATES (OVERALL)
ggplot(betas.df, aes(value, variable, alpha = as.numeric(p.inclusion))) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "turquoise",
    color = "white") +
  xlim(-5,5) +
  geom_vline(xintercept = 0)

ggplot( betas.df,
        aes( value,
             variable.simple,
             alpha = as.numeric(p.inclusion))) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "turquoise",
    color = "white") +
  xlim(-3,3) +
  geom_vline(xintercept = 0) + 
  facet_wrap(~ status + sex)


##---- PLOT COEFFICIENT ESTIMATES (MODEL-SPECIFIC)
reduced_betas.df <- betas.df[betas.df$model %in% reduced_aggr$model]
ggplot(reduced_betas.df, aes(value, variable, alpha = weight)) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "magenta",
    color = grey(1)) +
  xlim(-5,5) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ model)




## ------   3. MONTHLY SURVIVAL PROBABLILITIES -----
PHI <- array(NA, c(dim(res$sims.list$beta.temp)[1],2,2,2,n.months))
for(d in 1:2){
  for(s in 1:2){
    for(m in 1:n.months){
      PHI[ ,d,s,1,m] <- ilogit(res$sims.list$logit.phi0[ ,d,s] +
                                 res$sims.list$beta.temp[ ,d,s] * res$sims.list$z.temp[ ,d,s] * temp[m] +
                                 res$sims.list$beta.time[ ,d,s] * res$sims.list$z.time[ ,d,s] * month[m])
      PHI[ ,d,s,2,m] <- ilogit(res$sims.list$logit.phi0[ ,d,s] +
                                 res$sims.list$beta.temp[ ,d,s] * res$sims.list$z.temp[ ,d,s] * temp[m] +
                                 res$sims.list$beta.time[ ,d,s] * res$sims.list$z.time[ ,d,s] * month[m] +
                                 res$sims.list$beta.transloc[ ,d,s] * res$sims.list$z.transloc[ ,d,s])
    }#m
  }#ss
}#s
dimnames(PHI) <- list( 
  iteration = 1:dim(res$sims.list$beta.temp)[1],
  status = c("protected","degraded"),
  sex = c("female","male"),
  translocated = c("no","yes"),
  month = 1:n.months)

mean.PHI <- apply(PHI, c(2,3,4,5), mean)
upper.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.975))
lower.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.025))


##-- Identify first month w/ translocation
test <- cbind.data.frame( 
  transloc = nimConstants$transloc,
  sex = nimConstants$sex,
  status = nimConstants$status,
  startM = nimConstants$start.int[ ,1]) %>%
  filter(., transloc == 2)

start <- matrix(NA,2,2)
start[1,1] <- min(test$startM[test$status == 1 & test$sex == 1])
start[1,2] <- min(test$startM[test$status == 1 & test$sex == 2])
start[2,1] <- min(test$startM[test$status == 2 & test$sex == 1])
start[2,2] <- min(test$startM[test$status == 2 & test$sex == 2])

##-- Plot Survival
par(mfrow = c(2,2))
sex <- c("females", "males")
status <- c("degraded", "protected")
for(f in 1:2){
  for(s in 1:2){
    plot(1, type = "n",
         xlim = c(0,n.months+1),
         ylim = c(0, 1),
         axes = F,
         ylab = "Survival prob.",
         xlab = "Months",
         main = paste0(status[f],"-", sex[s]))
    axis(1, at = seq(0,250,50), labels = seq(0,250,50))
    axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
    legend( x = 1, y = 0.2,
            title = "Individual status:",
            legend = c("resident", "translocated"),
            bty = "n",
            fill = myCols[c(2,4)])
    
    polygon(x = c(1:n.months,n.months:1),
            y = c(upper.PHI[f,s,1, ],rev(lower.PHI[f,s,1, ])),
            col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
    points(1:n.months, mean.PHI[f,s,1,1:n.months], type = "l", lwd = 2, col = myCols[2])
    
    polygon(x = c(start[f,s]:n.months,n.months:start[f,s]),
            y = c( upper.PHI[f,s,2,start[f,s]:n.months],
                   rev(lower.PHI[f,s,2,start[f,s]:n.months])),
            col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
    points(start[f,s]:n.months, mean.PHI[f,s,2,start[f,s]:n.months], type = "l", lwd = 2, col = myCols[4])
  }#s
}#f



## ------   4. ANNUAL SURVIVAL PROBABILITIES -----
MPHI <- array(NA, c(dim(res$sims.list$beta.temp)[1],2,n.months))
mean.PHI <- upper.PHI <- lower.PHI <- array(NA, c(2,2,2,n.months-11))
frag <- c("protected","disturbed")
sex <- c("females", "males")
for(d in 1:2){
  for(s in 1:2){
    print(paste0("###  Processing survival for ", sex[s], " in ", frag[d], " fragments!  ###" ))
    for(m in 1:n.months){
      MPHI[ ,1,m] <- ilogit(res$sims.list$logit.phi0[ ,d,s] +
                              res$sims.list$beta.temp[ ,d,s] * res$sims.list$z.temp[ ,d,s] * temp[m] +
                              res$sims.list$beta.time[ ,d,s] * res$sims.list$z.time[ ,d,s] * month[m])
      MPHI[ ,2,m] <- ilogit(res$sims.list$logit.phi0[ ,d,s] +
                              res$sims.list$beta.temp[ ,d,s] * res$sims.list$z.temp[ ,d,s] * temp[m] +
                              res$sims.list$beta.time[ ,d,s] * res$sims.list$z.time[ ,d,s] * month[m] +
                              res$sims.list$beta.transloc[ ,d,s] * res$sims.list$z.transloc[ ,d,s])
    }#m
    for(m in 1:(n.months-11)){
      PHI1 <- apply(MPHI[ ,1,m:(m+11)],1,prod)
      PHI2 <- apply(MPHI[ ,2,m:(m+11)],1,prod)
      
      mean.PHI[d,s,1,m] <- mean(PHI1)
      upper.PHI[d,s,1,m] <- quantile(PHI1,0.975)
      lower.PHI[d,s,1,m] <- quantile(PHI1,0.025)
      
      mean.PHI[d,s,2,m] <- mean(PHI2)
      upper.PHI[d,s,2,m] <- quantile(PHI2,0.975)
      lower.PHI[d,s,2,m] <- quantile(PHI2,0.025)
    }#m
    print(paste0("sex : ",s))
  }#s
  print(paste0("status : ",d))
}#d


##-- Identify first month w/ translocation
test <- cbind.data.frame( 
  transloc = nimConstants$transloc,
  sex = nimConstants$sex,
  status = nimConstants$status,
  startM = nimConstants$start.int[ ,1]) %>%
  filter(., transloc == 2)

start <- matrix(NA,2,2)
start[1,1] <- min(test$startM[test$status == 1 & test$sex == 1])
start[1,2] <- min(test$startM[test$status == 1 & test$sex == 2])
start[2,1] <- min(test$startM[test$status == 2 & test$sex == 1])
start[2,2] <- min(test$startM[test$status == 2 & test$sex == 2])



par(mfrow = c(2,2))
n.months = dim(lower.PHI)[4]

##-- Females protected plot 
par(mar = c(0,4,6,0))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0,0.6), axes = F,
     ylab = "Survival probability",
     xlab = "",
     main = "Females")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,0.6,0.1), labels = seq(0,0.6,0.1))
polygon(x = c(start[2,1]:n.months,n.months:start[2,1]),
        y = c( upper.PHI[2,1,2,start[2,1]:n.months],
               rev(lower.PHI[2,1,2,start[2,1]:n.months])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points( start[2,1]:n.months, mean.PHI[2,1,2,start[2,1]:n.months],
        type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[2,1,1, ],rev(lower.PHI[2,1,1, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[2,1,1, ], type = "l", lwd = 3, col = myCols[4])


##-- Males protected plot 
par(mar = c(0,2,6,2))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 0.6), axes = F,
     ylab = "",
     xlab = "",
     main = "Males")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,0.6,0.1), labels = seq(0,0.6,0.1))
polygon(x = c(start[2,2]:n.months,n.months:start[2,2]),
        y = c(upper.PHI[2,2,2,start[2,2]:n.months],
              rev(lower.PHI[2,2,2,start[2,2]:n.months])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(start[2,2]:n.months, mean.PHI[2,2,2,start[2,2]:n.months],
       type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[2,2,1, ],rev(lower.PHI[2,2,1, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[2,2,1, ], type = "l", lwd = 3, col = myCols[4])


##-- Females fragmented plot 
par(mar = c(4,4,2,0))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 0.6), axes = F,
     ylab = "Survival probability",
     xlab = "",
     main = "")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,0.6,0.1), labels = seq(0,0.6,0.1))
polygon(x = c(start[1,1]:n.months,n.months:start[1,1]),
        y = c(upper.PHI[1,1,2,start[1,1]:n.months],
              rev(lower.PHI[1,1,2,start[1,1]:n.months])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(start[1,1]:n.months, mean.PHI[1,1,2,start[1,1]:n.months],
       type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[1,1,1, ],rev(lower.PHI[1,1,1, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[1,1,1, ], type = "l", lwd = 3, col = myCols[4])


##-- Males fragmented plot 
par(mar = c(4,2,2,2))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 0.6), axes = F,
     ylab = "",
     xlab = "",
     main = "")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,0.6,0.1), labels = seq(0,0.6,0.1))
polygon(x = c(start[1,2]:n.months,n.months:start[1,2]),
        y = c(upper.PHI[1,2,2,start[1,2]:n.months],
              rev(lower.PHI[1,2,2,start[1,2]:n.months])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(start[1,2]:n.months, mean.PHI[1,2,2,start[1,2]:n.months],
       type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[1,2,1, ],rev(lower.PHI[1,2,1, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[1,2,1, ], type = "l", lwd = 3, col = myCols[4])





## ------   5. SURVIVAL TEMPERATURE PLOTS -----
temp <- seq(min(nimConstants$temp), max(nimConstants$temp), by = 0.05)
n.temp <- length(temp)
month <- mean(nimConstants$month)

##-- Plot Survival
par(mfrow = c(2,2))
sex <- c("females", "males")
status <- c("Disturbed", "Protected")
for(d in 1:2){
  for(s in 1:2){
    plot(1, type = "n",
         xlim = c(min(nimConstants$temp),max(nimConstants$temp)+0.2),
         ylim = c(0, 1),
         axes = F,
         ylab = "Survival prob.",
         xlab = "Temperature",
         main = paste0(status[d],"-", sex[s]))
    axis(1,
         at = seq(min(nimConstants$temp)+0.1*0.4385646,
                  max(nimConstants$temp)+0.5*0.4385646,
                  by = 0.4385646),
         labels = seq(20, 29, by = 1))
    axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
    legend( x = 1, y = 0.2,
            title = "Individual status:",
            legend = c("resident", "translocated"),
            bty = "n",
            fill = myCols[c(2,4)])
    
    ## Calculate PHI
    PHI.resid <- PHI.trns <- matrix(NA,dim(res$sims.list$beta.temp)[1],n.temp)
    for(t in 1:n.temp){
      PHI.resid[ ,t] <- ilogit(res$sims.list$logit.phi0[ ,d,s] +
                                 res$sims.list$beta.temp[ ,d,s] * res$sims.list$z.temp[ ,d,s] * temp[t] +
                                 res$sims.list$beta.time[ ,d,s] * res$sims.list$z.time[ ,d,s] * month)
      
      PHI.trns[ ,t] <- ilogit(res$sims.list$logit.phi0[ ,d,s] +
                                res$sims.list$beta.temp[ ,d,s] * res$sims.list$z.temp[ ,d,s] * temp[t] +
                                res$sims.list$beta.time[ ,d,s] * res$sims.list$z.time[ ,d,s] * month +
                                res$sims.list$beta.transloc[ ,d,s] * res$sims.list$z.transloc[ ,d,s])
    }#t
    
    mean.PHI.resid <- apply(PHI.resid, 2, mean)
    upper.PHI.resid  <- apply(PHI.resid, 2, function(x)quantile(x,0.975))
    lower.PHI.resid  <- apply(PHI.resid, 2, function(x)quantile(x,0.025))
    
    mean.PHI.trns <- apply(PHI.trns, 2, mean)
    upper.PHI.trns  <- apply(PHI.trns, 2, function(x)quantile(x,0.975))
    lower.PHI.trns  <- apply(PHI.trns, 2, function(x)quantile(x,0.025))
    
    polygon(x = c(temp,rev(temp)),
            y = c(upper.PHI.trns,rev(lower.PHI.trns)),
            col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
    points(temp, mean.PHI.trns, type = "l", lwd = 2, col = myCols[2])
    
    polygon(x = c(temp,rev(temp)),
            y = c(upper.PHI.resid,rev(lower.PHI.resid)),
            col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
    points(temp, mean.PHI.resid, type = "l", lwd = 2, col = myCols[4])
  }#s
}#f



## ------   6. DETECTION PROBABILITY -----
## Individual detection
temp2 <- seq(min(nimConstants$temp2, na.rm = T), 
             max(nimConstants$temp2, na.rm = T),
             0.05)
n.sites <- nimConstants$n.sites
n.iter <- dim(res$sims.list$lambda0)[1]
P <- array(NA, c(n.iter,n.sites,2,length(temp2)))
for(t in 1:length(temp2)){
  for(s in 1:2){
    for(f in 1:n.sites){
      P[ ,f,s,t] <- 1-exp(-exp( res$sims.list$lambda0[ ,f,s] +
                                  res$sims.list$gamma.temp[ ,f,s] *
                                  res$sims.list$z.det[ ,f,s] *
                                  temp2[t]) * 4) 
    }#f
  }#s
}#t
mean.P <- apply(P, c(2,3,4), function(x)mean(x, na.rm = T))
upper.P <- apply(P, c(2,3,4), function(x)quantile(x, 0.975, na.rm = T))
lower.P <- apply(P, c(2,3,4), function(x)quantile(x, 0.025, na.rm = T))


par(mfrow = c(4,2))
sex <- c("females", "males")
for(f in 1:n.sites){
  plot(1, type = "n",
       xlim = c(min(temp2)-0.05,max(temp2)+0.05),
       ylim = c(0,1), axes = F,
       ylab = "Detection prob.",
       xlab = "Temperature",
       main = fragments[f])
  axis(1,
       at = seq(min(temp2)-0.05,max(temp2),0.5),
       labels = seq(min(temp2)-0.05,max(temp2),0.5))
  axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
  legend( x = 26, y = 1,
          legend = c("females", "males"),
          bty = "n",
          fill = myCols[c(1,3)])
  
  polygon(x = c(temp2,rev(temp2)),
          y = c(upper.P[f,1, ],rev(lower.P[f,1, ])),
          col = adjustcolor(myCols[1],alpha.f = 0.5), border = F)
  polygon(x = c(temp2,rev(temp2)),
          y = c(upper.P[f,2, ],rev(lower.P[f,2, ])),
          col = adjustcolor(myCols[3],alpha.f = 0.5), border = F)
  
  points(temp2, mean.P[f,1, ], type = "l", lwd = 3, col = myCols[1])
  points(temp2, mean.P[f,2, ], type = "l", lwd = 3, col = myCols[3])
}#d



## ------   7. TABLES -----
##---- Calculate "phi0" from "lphi0"
tmp <- betas.df[betas.df$variable.simple == "lphi0", ]
tmp$value <- ilogit(tmp$value)
tmp$variable.simple <- "phi0"
betas.df <- rbind(betas.df,tmp)


##---- Summarize parameter estimates
getCleanEstimates <- function(x,
                              moment = "mean",
                              quantiles = c(0.025,0.975)){
  if(moment == "mean"){
    paste0(format(round(mean(x),digits = 2),nsmall = 2)," (",
           format(round(quantile(x, probs = quantiles[1]), digits = 2),nsmall = 2), "-",
           format(round(quantile(x, probs = quantiles[2]), digits = 2),nsmall = 2), ")")
  } else {
    paste0(format(round(median(x), digits = 2),nsmall = 2)," (",
           format(round(quantile(x, probs = quantiles[1]), digits = 2),nsmall = 2), "-",
           format(round(quantile(x, probs = quantiles[2]), digits = 2),nsmall = 2), ")")
  }# else
}

summary <- betas.df %>% 
  group_by(variable.simple,sex,status) %>%
  summarize( estimate = getCleanEstimates(value))


##---- Create empty table
survival <-  matrix(NA, nrow = 5, ncol = 4)
rownames(survival) <- c("","phi0", "beta.time", "beta.temp", "beta.transloc")
colnames(survival) <- c("protected", "protected", "fragmented", "fragmented")
survival[1, ] <- c("female", "male", "female", "male")
for(s in c("female","male")){
  for(d in c("protected","fragmented")){
    
    thisCol <- colnames(survival) == d & survival[1, ] == s
    tmp <- filter(summary, sex == s, status == d)
    
    survival[2,thisCol] <- tmp$estimate[tmp$variable.simple == "phi0"]
    survival[3,thisCol] <- tmp$estimate[tmp$variable.simple == "time"]
    survival[4,thisCol] <- tmp$estimate[tmp$variable.simple == "temp"]
    survival[5,thisCol] <- tmp$estimate[tmp$variable.simple == "transloc"]
  }
}

##----- Export as .csv
# write.csv( survival,
#            file =  file.path(analysisDir, modelName, "Params.csv"))




## -----------------------------------------------------------------------------

