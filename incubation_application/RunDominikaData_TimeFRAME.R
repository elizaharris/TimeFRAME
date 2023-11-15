
rm(list = ls())
require(lubridate)
library(TimeFRAME)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
knitr::opts_chunk$set(
  fig.width = 8,
  fig.height = 5,
  fig.align = "center",
  out.width = "100%",
  comment = "#>"
)
set.seed(1)

# Load the data 
data = read.csv("data_forTimeFrame.csv") 
data$dt = dmy_hm(paste(data$Date,data$Time,sep=" "))
for (c in colnames(data)){
  if (is.numeric(data[1,c])){
    data[data[,c]<(-1000),c] = NaN
}}
iso = data[,c("d15N","SP","d18O")]
sd = data[,c("d15N_sd","SP_sd","d18O_sd")]
L1 = which(data$Experiment=="L1_dryer") # The different experiments from Dominika
L2 = which(data$Experiment=="L2_wetter")

# Create relative times as numeric for FRAME input
data$dt_diff = -1
data$dt_diff[L1] = as.numeric(data$dt[L1] - min(data$dt[L1]))/60/60/24 # t diff in days
data$dt_diff[L2] = as.numeric(data$dt[L2] - min(data$dt[L2]))/60/60/24 # t diff in days

# Load the model
sources = read.csv("endmembers.csv",row.names = 1) # Use 3DIM/FRAME endmembers
sources = sources[,c(1,3,5,2,4,6)] # reorganise to iso then sd
rownames(sources)[rownames(sources)=="bd"] = "bD"
model <- frame_model(sources = sources[1:4,], frac = sources[5,])

# Independent fits
t1 = Sys.time()
fit.independent.L1 <- fit_frame(model, iso[L1,], sd=sd[L1,], t = data$dt_diff[L1])
autoplot(fit.independent.L1) # autoplot to look
head(coef(fit.independent.L1)) # look at head
head(as.data.frame(fit.independent.L1))
fit.independent.L2 <- fit_frame(model, iso[L2,], sd=sd[L2,], t = data$dt_diff[L2])
autoplot(fit.independent.L2) # autoplot to look
t2 = Sys.time()
time.ind = t2 - t1; print(time.ind)

# Spline fits (need to tune M, M.r)
t1 = Sys.time()
fit.splines.L1 <- fit_glm(model, iso[L1,], sd=sd[L1,], t = data$dt_diff[L1], M = 6, M.r = 3)
autoplot(fit.splines.L1)
fit.splines.L2 <- fit_glm(model, iso[L2,], sd=sd[L2,], t = data$dt_diff[L2], M = 6, M.r = 3)
autoplot(fit.splines.L2)
t2 = Sys.time()
time.spline = t2 - t1; print(time.ind)

# HDGP fit, est. correlation length (tune rho, rho.r)
t1 = Sys.time()
fit.hdgp.L1 <- fit_dgp(model, iso[L1,], sd=sd[L1,], t = data$dt_diff[L1], rho = 0.2, rho.r = 0.5, estim.rho = TRUE)
autoplot(fit.hdgp.L1)
fit.hdgp.L2 <- fit_dgp(model, iso[L2,], sd=sd[L2,], t = data$dt_diff[L2], rho = 0.2, rho.r = 0.5, estim.rho = TRUE)
autoplot(fit.hdgp.L2)
t2 = Sys.time()
time.dgp = t2 - t1; print(time.ind)

# Reorganise the final data
models = c("ind","spline","dgp")
pathways = rownames(sources)
pathres = matrix(nrow=dim(iso)[1],ncol=length(models)*length(pathways)*2)*NaN # n pathways * 2 values (mean, sd) * n models tested
n1 = rep(pathways, times=length(models)*2)
n2 = rep(c("mean","sd"),each=length(pathways))
n3 = rep(models, each=length(pathways)*2)
colnames(pathres) = paste(n1, n2, n3, sep="_")
for (m in models){
  if (m=="ind"){ data_L1 = fit.independent.L1; data_L2 = fit.independent.L2 }
  if (m=="spline"){ data_L1 = fit.splines.L1; data_L2 = fit.splines.L2 }
  if (m=="dgp"){ data_L1 = fit.hdgp.L1; data_L2 = fit.hdgp.L2 }
  data_L1 = as.data.frame(data_L1)
  data_L2 = as.data.frame(data_L2)
  for (f in pathways){ # For each param
      pathres[L1,paste(f, "mean", m, sep="_")] = data_L1$mean[data_L1$source==f]
      pathres[L2,paste(f, "mean", m, sep="_")] = data_L2$mean[data_L1$source==f]
      pathres[L1,paste(f, "sd", m, sep="_")] = data_L1$sd[data_L1$source==f]
      pathres[L2,paste(f, "sd", m, sep="_")] = data_L2$sd[data_L1$source==f]
  }
}
write.csv(pathres, "timeframe_pathways.csv", row.names=FALSE)

# Look at data and compare to other approaches

# First, split out all the datasets into the same format
tf_ind = pathres[,paste(n1, n2, "ind", sep="_")]; colnames(tf_ind) = paste(n1, n2, sep="_")
tf_spline = pathres[,paste(n1, n2, "spline", sep="_")]; colnames(tf_spline) = paste(n1, n2, sep="_")
tf_dgp = pathres[,paste(n1, n2, "dgp", sep="_")]; colnames(tf_dgp) = paste(n1, n2, sep="_")
# Create temp data and fill for the other estimates
data_temp = data.frame(matrix(nrow=dim(data)[1],ncol=length(pathways)*2))
n1 = rep(pathways, times=2)
n2 = rep(c("mean","sd"),each=length(pathways))
colnames(data_temp) = paste(n1, n2, sep="_")
FRAME = data_temp
FRAME[,] = data[paste("FRAME", c(pathways[1:4],"r"), "N2O", n2, sep="_")]
SPOmap.RM = data_temp; SPOmap.MR = data_temp;  
SPOmap.RM[,c("bD_mean","red_mean")] = data[,c("f_bD_SP.O.map.RM","r_N2O_SP.O.map.RM")]
SPOmap.MR[,c("bD_mean","red_mean")] = data[,c("f_bD_SP.O.map.MR","r_N2O_SP.O.map.MR")]
SPOmap.RM.mean = data_temp; SPOmap.MR.mean = data_temp; 
SPOmap.RM.mean[,c("bD_mean","red_mean")] = data[,c("f_bD_SP.O.map.RM_mean","r_N2O_SP.O.map.RM_mean")]
SPOmap.MR.mean[,c("bD_mean","red_mean")] = data[,c("f_bD_SP.O.map.MR_mean","r_N2O_SP.O.map.MR_mean")]
SPOmap.RM.mean[,c("bD_sd","red_sd")] = data[,c("f_bD_SP.O.map.RM_sd","r_N2O_SP.O.map.RM_sd")]
SPOmap.MR.mean[,c("bD_sd","red_sd")] = data[,c("f_bD_SP.O.map.MR_sd","r_N2O_SP.O.map.MR_sd")]
gasflux15N.mean = data_temp; 
gasflux15N.mean[,c("bD_mean","red_mean")] = data[,c("f_N2O.NO3_15Ngasflux_mean","r_N2O_15Ngasflux_mean")]
gasflux15N.mean[,c("bD_sd","red_sd")] = data[,c("f_N2O.NO3_15Ngasflux_sd","r_N2O_15Ngasflux_sd")]
# For this, bD is actually all N2O from NO3 precursors...

# Compare to FRAME
pdf(file = "Comparison.pdf",width=9,height=5)
par(mfrow=c(2,3),mar=c(2,4,1,1))
for (f in pathways[c(1,5)]){ # For simplicity, only show bD and red
  for (mod in c("FRAME","SP-Omap","15Ntracer")){
    if (mod=="FRAME"){moddata = FRAME} # Choose what to plot against
    if (mod=="SP-Omap"){moddata = SPOmap.MR.mean}
    if (mod=="15Ntracer"){moddata = gasflux15N.mean} 
    # Plot results
    if (sum(!is.na(moddata[,paste0(f,"_mean")]))==0){ # If no data for this model/pathway
      plot(1,1); next
    } else {
      plot(moddata[,paste0(f,"_mean")],tf_ind[,paste0(f,"_mean")], type="n", xlim=c(0,1),ylim=c(0,1),
         xlab=paste(f, mod, sep="_"),ylab=paste(f, "TimeFRAME", sep="_")) 
    }
    # independent frame
    points(moddata[L1,paste0(f,"_mean")],tf_ind[L1,paste0(f,"_mean")], pch=18, col=alpha("darkorange1", 0.6))
    points(moddata[L2,paste0(f,"_mean")],tf_ind[L2,paste0(f,"_mean")], pch=20, col=alpha("darkorange1", 0.6))
    points(moddata[L1,paste0(f,"_mean")],tf_spline[L1,paste0(f,"_mean")], pch=18, col=alpha("blue4", 0.6))
    points(moddata[L2,paste0(f,"_mean")],tf_spline[L2,paste0(f,"_mean")], pch=20, col=alpha("blue4", 0.6))
    points(moddata[L1,paste0(f,"_mean")],tf_dgp[L1,paste0(f,"_mean")], pch=18, col=alpha("chartreuse3", 0.4))
    points(moddata[L2,paste0(f,"_mean")],tf_dgp[L2,paste0(f,"_mean")], pch=20, col=alpha("chartreuse3", 0.4))
    lines(c(0,1),c(0,1),col="grey")
    # Statistical comparison
    mod1 = summary(lm(tf_ind[,paste0(f,"_mean")]~moddata[,paste0(f,"_mean")]))
    mad1 = mean(abs(tf_ind[,paste0(f,"_mean")]-moddata[,paste0(f,"_mean")]),na.rm=TRUE)
    lines(c(0,1),c(0,1)*mod1$coefficients[2,1]+mod1$coefficients[1,1], col="darkorange1")
    mod2 = summary(lm(tf_spline[,paste0(f,"_mean")]~moddata[,paste0(f,"_mean")]))
    mad2 = mean(abs(tf_spline[,paste0(f,"_mean")]-moddata[,paste0(f,"_mean")]),na.rm=TRUE)
    lines(c(0,1),c(0,1)*mod2$coefficients[2,1]+mod2$coefficients[1,1], col="blue4")
    mod3 = summary(lm(tf_dgp[,paste0(f,"_mean")]~moddata[,paste0(f,"_mean")]))
    mad3 = mean(abs(tf_dgp[,paste0(f,"_mean")]-moddata[,paste0(f,"_mean")]),na.rm=TRUE)
    lines(c(0,1),c(0,1)*mod3$coefficients[2,1]+mod3$coefficients[1,1], col="chartreuse3")
    legend("topleft", legend=c(round(mad1,2),round(mad2,2),round(mad3,2)),col=c("darkorange1","blue4", "chartreuse3"), pch=1, cex=0.5) 
  }}
dev.off()

# Plot against WFPS
pdf(file = "WFPS.pdf",width=8,height=6)
par(mfrow=c(2,1),mar=c(2,4,1,1))
for (f in pathways[c(1,5)]){
  plot(data$WFPS,tf_spline[,paste0(f,"_mean")],col="blue4",ylab=f,ylim=c(0,1),pch=20)
  arrows(data$WFPS,tf_spline[,paste0(f,"_mean")]-tf_spline[,paste0(f,"_sd")],
         data$WFPS,tf_spline[,paste0(f,"_mean")]+tf_spline[,paste0(f,"_sd")],length=0.05,angle=90,code=3,col="blue4")
  points(data$WFPS,tf_dgp[,paste0(f,"_mean")],col="chartreuse3",pch=20)
  arrows(data$WFPS,tf_dgp[,paste0(f,"_mean")]-tf_dgp[,paste0(f,"_sd")],
         data$WFPS,tf_dgp[,paste0(f,"_mean")]+tf_dgp[,paste0(f,"_sd")],length=0.05,angle=90,code=3,col="chartreuse3")
  # Find weighted fits
  mod1 = summary(lm(tf_spline[,paste0(f,"_mean")]~data$WFPS),weights=1/tf_spline[,paste0(f,"_sd")])
  mod1str = paste0(f," = WFPS * ",round(mod1$coefficients[2,1],2)," + ",round(mod1$coefficients[1,1],2),"; R2 = ",round(mod1$r.squared,2))
  lines(c(55,90),c(55,90)*mod1$coefficients[2,1]+mod1$coefficients[1,1],col="blue4")
  mod2 = summary(lm(tf_dgp[,paste0(f,"_mean")]~data$WFPS),weights=1/tf_spline[,paste0(f,"_sd")])
  mod2str = paste0(f," = WFPS * ",round(mod2$coefficients[2,1],2)," + ",round(mod2$coefficients[1,1],2),"; R2 = ",round(mod2$r.squared,2))
  lines(c(55,90),c(55,90)*mod2$coefficients[2,1]+mod2$coefficients[1,1],col="chartreuse3")
  print(mod1); print(mod2)
  legend("topleft", legend=c(mod1str,mod2str),col=c("blue4", "chartreuse3"), pch=1, cex=0.6) 
}
dev.off()
  

