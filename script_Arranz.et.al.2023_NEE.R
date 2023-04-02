#=================================#

#'@Title Analysis of the temporal trends of fish size spectra in stream communities across France

#'@Description This R script is a brief overview of the main findings of the manuscript "Modulating 
#              effects of human pressures on climate warming-induced changes in size spectra of stream 
#             fish communities" published in Nature Ecology and Evolution

#'@Contact ignasi.arranz at urjc.es for further information on the details of the results

#'@Data_information: 
#'    @import_data: data frame with all variables
#'        @Code_ID: unique identification of the stream location
#'        @y.WGS84: latitudinal coordinates in WGS84
#'        @x.WGS84: longitudinal coordinates in WGS84
#'        @First_samp: first year of the fish survey
#'        @Num_samp: number of times surveyed te fish community
#'        @Trend_size_spectra: annual changes of size spectrum slopes
#'        @Updow_grad: upstream-downstream gradient
#'        @Clim_cond: climatic conditions in ºC
#'        @Clim_warm: climate warming in ºC/year
#'        @Hum_pres: human pressures

#=================================#

#Set the environment

rm(list = ls(all.names = TRUE))
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/tricomicet/Postdoc Toulouse/Manuscript temporal dynamics streams/")

#Load packages

packages <- c("nlme", "MASS", "wesanderson", "scales", "effects")
install.packages(setdiff(packages, rownames(installed.packages())))

library(nlme)
library(MASS)
library(wesanderson)
library(scales)
library(effects)

#Load data

import_data <- read.table(file = 'data_Arranz.et.al.2023_NEE', sep = ';',dec=".", header = TRUE,fill = TRUE)

#Transform data

import_data$Hum_pres_log <- log10(import_data$Hum_pres)
  
#Run Generalized Mixed Model

full.model <- lme(Trend_size_spectra ~
                    
                    scale(Updow_grad)*scale(Clim_warm) +
                    
                    scale(Clim_cond)*scale(Clim_warm) +
                    
                    scale(Hum_pres_log)*scale(Clim_warm),
                  
                    random= ~1|First_samp/Num_samp,
                  
                    method = "ML", data= import_data)

step.save <- stepAIC(full.model,direction=c("backward"))#

best.model <- lme(Trend_size_spectra ~ 
                    
                    scale(Updow_grad) + 
                    
                    scale(Clim_warm) + 
                    
                    scale(Hum_pres_log) + 
                  
                    scale(Clim_warm):scale(Hum_pres_log),
                  
                    random= ~1|First_samp/Num_samp,
                  
                    method = "ML", data= import_data) 

summary(best.model)#Outcome provided in Supplementary Table 

#Illustrate Figure 2
pdf(file = "~/Desktop/test.pdf",   # The directory you want to save the file in
    paper="a4r",width=18.27,height=21.69)#"a4r" for rotated (landscape)

par(mar=c(4,5,4,4),mfrow=c(1,2),pty="s")

#Figure 2A
pred <- as.data.frame(expand.grid(Clim_warm=seq(-0.03, 0.14, 0.001), 
                                  Hum_pres_log=seq(0.07, 1.699, 0.01)))#It is 10^1.699 because the max value that GHFI can take is 50

pred$Updow_grad <- mean(import_data$Updow_grad) 
pred$size_spectra_fit <- predict(best.model, pred,level=0)#check this on why it is important to add level= 0: https://stats.stackexchange.com/questions/58031/prediction-on-mixed-effect-models-what-to-do-with-random-effects
len <- length(unique(pred$Hum_pres_log))
pale <- wes_palette("Zissou1", len, type = "continuous")

plot(size_spectra_fit~Clim_warm,
     data=pred,type="n",
     xlab="Climate warming (°C/year)",
     ylab="Annual trends size spectrum slopes",
     ylim=c(-0.1,0.1),
     xlim=c(-0.03,0.14),
     cex.lab=1.3,
     cex.axis=1.1)

colfunc <- colorRampPalette(c(pale[len],pale[len/2], pale[1]))
legend_image <- as.raster(matrix(colfunc(len), ncol=1))
rasterImage(legend_image, -0.02, 0.09, -0.01, 0.05)     
mtext('Human \n pressures', 3, -1.6, cex=.6, adj=.1, font=2)
mtext('(A)', side=3, line=-1.5, at=0.135,cex=1.5)
pale2 <- wes_palette("Zissou1", 10, type = "continuous")
pal =  colorRampPalette(c(pale2[1], pale2[5],pale2[10]))
import_data$ord_hum_pres = findInterval(import_data$Hum_pres, sort(import_data$Hum_pres))
col_hum_pres <- alpha(pal(nrow(import_data))[import_data$ord_hum_pres],0.8)
points(import_data$Trend_size_spectra ~ import_data$Clim_warm,pch=21, bg= alpha(col_hum_pres,0.3),col= alpha("black", 0.0),cex=2)

for (i in 1:len){
  subset <- pred[which(pred$Hum_pres_log == unique(pred$Hum_pres_log)[i],),]
  points(subset$size_spectra_fit~subset$Clim_warm,type="l",col=pale[i])
}

abline(h=0,lty=2)

#Figure 2B
effx1 <- effect("scale(Updow_grad)", best.model,residuals=T)
df <- data.frame(updow_grad=effx1$x.all,residuals=effx1$residuals)

plot(df$Updow_grad,df$residuals,
     xlab="PCA Axis \nUpstream-downstream gradient",
     ylab="Annual trends size spectrum slopes",
     ylim=c(-0.1,0.1),
     cex.lab=1.3,
     cex.axis=1.1,
     cex=2, 
     pch=21,
     col="black",
     bg= alpha("grey",0.12))

mtext('(B)', side=3, line=-1.5, at=3.7,cex=1.5)
abline(h=0,lty=2)
abline(a=summary(best.model)$tTable[1,1],b=summary(best.model)$tTable[2,1],lty=1,lwd=4,col="black")

dev.off()

#=================================#
  