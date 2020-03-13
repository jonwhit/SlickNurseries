#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Basic non-parameric Bootstrap
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
library(scales) # For plotting
library(tidyverse)
# Setwd
setwd("/Users/jonathan/Google Drive/MS_SlickNurseries/_CURRENT-DRAFT/_SUBMISSION/InitialSubmission_NEE/Data")
dat.raw <- read.csv("SlickNurseries_*.csv")
head(dat.raw)

dat.raw <- dat.raw %>% 
  select(-transect) %>%
  arrange(sample, type)

head(dat.raw)

#########################################
samples = unique(dat.raw$sample)
samples
results = NULL

# choose sample
nboots = 10000
for(i in 1:length(samples)){
  dat = subset(dat.raw,sample==samples[i]) 
  dat$type <- as.factor(dat$type)
  types = unique(dat$type) 
  
  #-----------------------------
  # Do non-parametric bootsrap
  #-----------------------------
  
  boots  =  NULL # objects
  mean_diff <- vector(,nboots)
  median_diff <- vector(,nboots)
  for(k in 1:nboots){
    # Subset type j
    mu.temp = NULL # random data object
    dat.samp = dat
    dat.samp$type <- sample(dat.samp$type,size = length(dat.samp$type),replace = F)
    mean_diff[k] <- diff(tapply(dat.samp$density,dat.samp$type,mean))
    median_diff[k] <- diff(tapply(dat.samp$density,dat.samp$type,median))
    
    for(j in 1:length(types)){
      # resample vector by type with replace and take mean
      mu.temp =  c(mu.temp,mean(sample(dat[dat$type==types[j],3],replace=TRUE)))
    }
    boots = rbind(boots,as.numeric(mu.temp))
  }
  boots= data.frame(boots)
  colnames(boots) = paste(types)
  
  # Ratio Calculation
  boots$Ratio <- (boots$Inside/mean(boots$Outside))
  
  # get stats 
  median = apply(boots,2,median, na.rm=TRUE)
  mean = apply(boots,2,mean, na.rm=TRUE)
  sd = apply(boots,2,sd, na.rm=TRUE)
  se = sd/sqrt(nboots)
  lci = apply(boots,2,quantile,c(0.025), na.rm=TRUE)
  uci = apply(boots,2,quantile,c(0.975), na.rm=TRUE)
  
  org_diff_mean <- diff(tapply(dat$density,dat$type,mean))
  org_diff_median <- diff(tapply(dat$density,dat$type,median))
  
  out = data.frame(sample=samples[i],type=c("Inside","Outside","Ratio"), mean,median,sd,se,lci,uci,`p(mean_out<in)` = c(mean(mean_diff > org_diff_mean)))
  
  results = rbind(results,out)
  write.csv(boots, paste0(getwd(),"/",samples[i],"_boots.csv"),row.names=F)
}

# confirm with R bootstrap package
bootmean <- function(data,index){
  d <- data[index,]
  tapply(d$density,d$type,mean)
}

bootres <- boot::boot(dat,bootmean,R=10000)
boot::boot.ci(bootres,index=c(1,1))
boot::boot.ci(bootres,index=c(2,2))

##############################################
#Save Results
write.csv(results,"boot_results.csv",row.names = F)

#Results Output Metadata
#sample = taxon/group
#type = habitat (Inside=Slick; Outside=Ambient; Ratio=Slick/Ambient)
#mean = mean density (#/m3)
#median = median density (#/m3)
#sd = standard deviation of mean density
#se = standard error
#lci = lower confidence interval
#uci = upper confidence interval
#p.mean_out.in = empirical Probability [P(d ̅slick >d̅ ambient ]) that the mean density (d ̅)
  #inside slicks is greater than in ambient (outside) waters.
