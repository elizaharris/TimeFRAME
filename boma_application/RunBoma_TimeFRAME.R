
rm(list = ls())
require(lubridate)
library(TimeFRAME)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
knitr::opts_chunk$set(
  fig.width = 8,
  fig.height = 5,
  fig.align = "center",
  out.width = "100%",
  comment = "#>"
)
set.seed(1)

# Load the data 
bomadata = read.csv("Boma_Kenya_Short.csv") 
bomadata$dt = dmy(bomadata$date)
iso = bomadata[,c("d15N","d15NSP","d18O")]
sd = iso*0+1 # Set stdevs to 1 for now...
r = which(rowSums(is.na(iso)) == 0 & rowSums(is.na(sd)) == 0)
t = c(1:length(r)) # Create numbering to use instead of time 

# Basic plot of the data
par(mfrow=c(3,1),mar=c(2,4,1,1))
plot(t,iso$d15N[r],col="blue",ylab="d15N")
lines(t,iso$d15N[r],col="blue")
plot(t,iso$d15NSP[r],col="blue",ylab="d15NSP")
lines(t,iso$d15NSP[r],col="blue")
plot(t,iso$d18O[r],col="blue",ylab="d18O")
lines(t,iso$d18O[r],col="blue")

# Load the model
model <- frame_model(n2o_sources, frac = n2o_frac)

# Independent fits, d15N relative to soil (best case)
iso_rel = iso
iso_rel$d15N = iso$d15N - bomadata$d15N_soil
fit.independent.rel <- fit_frame(model, iso_rel[r,], sd=sd[r,], t = r)
autoplot(fit.independent.rel) # autoplot to look

# Independent fits, no fD, no d15N (to match Stephen's model as much as possible)
model <- frame_model(n2o_sources[1:2,c(2,3,5,6)], frac = n2o_frac[c(2,3,5,6)])
fit.independent.nofD <- fit_frame(model, iso[r,2:3], sd=sd[r,2:3], t = t)
autoplot(fit.independent.nofD) # autoplot to look

# Reorganise the final data
models = c("ind_rel","ind_nofD")
pathways = c("Ni","bD","fD","Red")
bomares = matrix(nrow=dim(iso)[1],ncol=length(models)*length(pathways)*2)*NaN # 4 pathways * 2 values (mean, sd) * 4 models tested
n1 = rep(pathways, times=length(models)*2)
n2 = rep(c("mean","sd"),each=length(pathways))
n3 = rep(models, each=length(pathways)*2)
colnames(bomares) = paste(n1, n2, n3, sep="_")
for (m in models){
  if (m=="ind_rel"){ data = fit.independent.rel; thisr = r }
  if (m=="ind_nofD"){ data = fit.independent.nofD; thisr = r }
  data = as.data.frame(data)
  for (f in pathways){ # For each param
    if (sum(data$source==f)>0){
      bomares[thisr,paste(f, "mean", m, sep="_")] = data$mean[data$source==f]
      bomares[thisr,paste(f, "sd", m, sep="_")] = data$sd[data$source==f]
    }
  }
}
write.csv(bomares, "Boma_Kenya_Short_TFoutput.csv", row.names=FALSE)

# Look at data and also compare to Stephen
stephen = bomadata[,c("dt","Reduction_Percent_RedMix_SP_18O","Nitrification_Percent_RedMix_SP_18O","Denitrification_Percent_RedMix_SP_18O","Reduction_Percent_MixRed_SP_18O","Nitrification_Percent_MixRed_SP_18O","Denitrification_Percent_MixRed_SP_18O")]
colnames(stephen) = c("dt","Red_RM","Ni_RM","bD_RM","Red_MR","Ni_MR","bD_MR")
# Stephen's "Red" term is the frac remaining and mine is the frac reduced, so adjust...
stephen$Red_MR = 100-stephen$Red_MR
stephen$Red_RM = 100-stephen$Red_RM

par(mfrow=c(3,2),mar=c(2,4,1,1))
for (f in c("Ni","bD","Red")){
  m = "ind_nofD" # TimeFRAME ind, no fungal denit, d18O and SP only (to match Stephen's)
  plot(stephen[,paste(f, "MR", sep="_")],bomares[,paste(f, "mean", m, sep="_")]*100, type="p", pch=1, col="red",
       xlab=paste(f, "stephen", sep="_"),ylab=paste(f, "TimeFRAME", sep="_"))
  points(stephen[,paste(f, "RM", sep="_")],bomares[,paste(f, "mean", m, sep="_")]*100, pch=1, col="blue")
  lines(stephen[,paste(f, "RM", sep="_")],stephen[,paste(f, "RM", sep="_")], type="l", col="black")
  if (f == "Ni"){
    legend("topright", legend=c(m,"MR","RM"),col=c("grey","red", "blue"), pch=1, cex=0.8)}
  
  m = "ind_rel" # TimeFRAME ind, d15N relative to soil
  plot(stephen[,paste(f, "MR", sep="_")],bomares[,paste(f, "mean", m, sep="_")]*100, type="p", pch=1, col="red",
       xlab=paste(f, "stephen", sep="_"),ylab=paste(f, "TimeFRAME", sep="_"))
  points(stephen[,paste(f, "RM", sep="_")],bomares[,paste(f, "mean", m, sep="_")]*100, pch=1, col="blue")
  lines(stephen[,paste(f, "RM", sep="_")],stephen[,paste(f, "RM", sep="_")], type="l", col="black")
  if (f == "Ni"){
    legend("topright", legend=c(m,"MR","RM"),col=c("grey","red", "blue"), pch=1, cex=0.8)}
}

# Look at the d15N, SP and d18O predictions from each model
iso_pred = bomadata[,c("d15NSP","d18O","d15N")]
for (m in models){
  # SP
  Ni = bomares[,paste("Ni_mean", m, sep="_")]*n2o_sources["Ni","d15NSP"] 
  bD = bomares[,paste("bD_mean", m, sep="_")]*n2o_sources["bD","d15NSP"]  
  fD = bomares[,paste("fD_mean", m, sep="_")]*n2o_sources["fD","d15NSP"] 
  sources = rowSums(cbind(Ni,bD,fD),na.rm=T) # Combined sources before reduction
  Red = bomares[,paste("Red_mean", m, sep="_")] 
  iso_pred[,paste("d15NSP",m,sep="_")] = sources + n2o_frac["Red","d15NSP"]*log(1-Red)
  # d15N
  Ni = bomares[,paste("Ni_mean", m, sep="_")]*n2o_sources["Ni","d15N"] 
  bD = bomares[,paste("bD_mean", m, sep="_")]*n2o_sources["bD","d15N"]  
  fD = bomares[,paste("fD_mean", m, sep="_")]*n2o_sources["fD","d15N"] 
  sources = rowSums(cbind(Ni,bD,fD),na.rm=T) # Combined sources before reduction
  Red = bomares[,paste("Red_mean", m, sep="_")] 
  iso_pred[,paste("d15N",m,sep="_")] = sources + n2o_frac["Red","d15N"]*log(1-Red)
  # d18O
  Ni = bomares[,paste("Ni_mean", m, sep="_")]*n2o_sources["Ni","d18O"] 
  bD = bomares[,paste("bD_mean", m, sep="_")]*n2o_sources["bD","d18O"]  
  fD = bomares[,paste("fD_mean", m, sep="_")]*n2o_sources["fD","d18O"] 
  sources = rowSums(cbind(Ni,bD,fD),na.rm=T) # Combined sources before reduction
  Red = bomares[,paste("Red_mean", m, sep="_")] 
  iso_pred[,paste("d18O",m,sep="_")] = sources + n2o_frac["Red","d18O"]*log(1-Red)
}
for (m in c("MR","RM")){
  # SP
  Ni = stephen[,paste("Ni", m, sep="_")]/100*34.45 # Endmembers that Stephen used
  bD = stephen[,paste("bD", m, sep="_")]/100*-3.9
  sources = rowSums(cbind(Ni,bD,fD),na.rm=T) # Combined sources before reduction
  Red = stephen[,paste("Red", m, sep="_")]/100
  iso_pred[,paste("d15NSP",m,sep="_")] = sources + -5.9*log(1-Red)
  # d18O
  Ni = stephen[,paste("Ni", m, sep="_")]/100*30.15
  bD = stephen[,paste("bD", m, sep="_")]/100*17.55 
  sources = rowSums(cbind(Ni,bD,fD),na.rm=T) # Combined sources before reduction
  Red = stephen[,paste("Red", m, sep="_")]/100
  iso_pred[,paste("d18O",m,sep="_")] = sources + n2o_frac["Red","d18O"]*log(1-Red)
}

all_SP = iso_pred[,grepl("d15NSP_", colnames(iso_pred))]
RMSE_SP = ( colSums((all_SP - iso_pred$d15NSP)^2,na.rm=T) / colSums(!is.na(all_SP)) )^0.5
print(RMSE_SP)

all_d15N = iso_pred[,grepl("d15N_", colnames(iso_pred))]
RMSE_d15N = ( colSums((all_d15N - iso_pred$d15N)^2,na.rm=T) / colSums(!is.na(all_d15N)) )^0.5
print(RMSE_d15N)

all_d18O = iso_pred[,grepl("d18O_", colnames(iso_pred))]
RMSE_d18O = ( colSums((all_d18O - iso_pred$d18O)^2,na.rm=T) / colSums(!is.na(all_d18O)) )^0.5
print(RMSE_d18O)

# Basic plot of Stephen's data
par(mfrow=c(3,1),mar=c(2,4,1,1))
plot(t,stephen$Ni_MR[r],col="blue",ylab="Ni",ylim=c(0,100))
plot(t,stephen$bD_MR[r],col="blue",ylab="bD",ylim=c(0,100))
plot(t,stephen$Red_MR[r],col="blue",ylab="Red",ylim=c(0,100))

