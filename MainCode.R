setwd("/home/joelsmith/Projects/brazil/data")

library(raster)
library(rgdal)
library(foreign)
library(RNetCDF)
library(ncdf4)
library(compare)
library(corrplot)
library(lmtest)
library(nlme)
library(forecast)

#-------------------------------------------------------------------------------
#------------------------Deforestation and IBIS Data-------------------------------------
#------------------------------------------------------------------------------

setwd("/home/joelsmith/Projects/brazil/data")

hist.tucurui = read.csv("formatted-historic-tucurui.csv")
hist.tapajos = read.csv("formatted-historic-tapajos.csv")
future.tapajos = read.csv("formatted-future-tapajos.csv")

hist.tucurui = hist.tucurui[complete.cases(hist.tucurui[,c("X","Y")]),]
hist.tapajos = hist.tapajos[complete.cases(hist.tapajos[,c("X","Y")]),]
future.tapajos = future.tapajos[complete.cases(future.tapajos[,c("X","Y")]),]

dim(hist.tucurui)
dim(hist.tapajos)
dim(future.tapajos)

tapajos.grass = read.csv("tapajos-grass.csv")
tapajos.forest = read.csv("tapajos-forest.csv")

tucurui.grass = read.csv("tucurui-grass.csv")
tucurui.forest = read.csv("tucurui-forest.csv")

dim(tapajos.grass)
dim(tapajos.forest)

dim(tucurui.grass)
dim(tucurui.forest)

tucurui.xy = c(-49.646667,-3.831667)
tapajos.xy = c(-56.2845,-4.603722)


#----------------------Future data Tapajos-------------------------

tapajos2030 = subset(future.tapajos,select = c("X","Y","Initial","Forest_2030","Nonfor_2030"))
tapajos2030$Year = rep(2030,nrow(tapajos2030))
tapajos2030$Month = rep(1,nrow(tapajos2030))
colnames(tapajos2030)[3:5] = c("State","Forest","Nonforest")
tapajos2030[,4:5] = tapajos2030[,4:5]/(10^6)
dim(tapajos2030)
head(tapajos2030)

tapajos2050 = subset(future.tapajos,select = c("X","Y","Initial","Forest_2050","Nonfor_2050"))
tapajos2050$Year = rep(2050,nrow(tapajos2050))
tapajos2050$Month = rep(1,nrow(tapajos2050))
colnames(tapajos2050)[3:5] = c("State","Forest","Nonforest")
tapajos2050[,4:5] = tapajos2050[,4:5]/(10^6)
dim(tapajos2050)
head(tapajos2050)

#-----------------Modeling for runoff---------------------------------

tapajos.runoff = read.csv("tapajos-runoff.csv")
tucurui.runoff = read.csv("tucurui-runoff.csv")

#tapajos.runoff = tapajos.runoff[-which(tapajos.runoff$Month%in%c(1,2,3,10,11,12)),]
#tapajos.runoff = tapajos.runoff[-which(tapajos.runoff$Month%in%c(4:9)),]


dim(tapajos.runoff)
dim(tucurui.runoff)

tapajos.runoff[,6:7] = tapajos.runoff[,6:7]*100
tucurui.runoff[,6:7] = tucurui.runoff[,6:7]*100

tucurui.time = (tucurui.runoff$Year-1985)*12 + tucurui.runoff$Month + 1
dam.dist.mat = t(tucurui.xy-t(data.matrix(tucurui.runoff[,c("X","Y")])))
dam.dist = 1/(sqrt(rowSums(dam.dist.mat^2))+0.01)
tucurui.runoff$Dam.dist = dam.dist
tucurui.runoff$Time.diff = tucurui.time

results.runoff <- data.analysis(tapajos.runoff,tucurui.runoff,"Qt","Nonforest")
hist.tapajos.model <- results.runoff$model1
hist.tucurui.model <- results.runoff$model2
runoff.summary = summary(hist.tapajos.model)
summary(hist.tucurui.model)
runoff.line = as.numeric(hist.tapajos.model$coefficients[1:2])
x.runoff = seq(0,100,by=0.1)
y.runoff = runoff.line[1] + runoff.line[2]*x
plot(x.runoff,y.runoff,type="l",xlab="Deforested (%)",ylab="Runoff")

predict.runoff <- results.runoff$predict.tapajos(tapajos2030)
tapajos2030$Qt = predict.runoff

proj =  crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

tapajos30runoff = tapajos2030[,c("X","Y","Qt")]
coordinates(tapajos30runoff) <- ~X+Y
projection(tapajos30runoff) = proj
r = raster(res = 1)
extent(r) = extent(tapajos30runoff)
r.runoff <- rasterize(tapajos30runoff, r, field = "Qt", fun = mean)

plot(r.runoff)


#---------------Modeling for ET---------------------------

tapajos.et = read.csv("tapajos-ET.csv")
tucurui.et = read.csv("tucurui-ET.csv")

#tapajos.et = tapajos.et[-which(tapajos.et$Month%in%c(1,2,3,10,11,12)),]
#tapajos.et = tapajos.et[-which(tapajos.et$Month%in%c(4:9)),]


dim(tapajos.et)
dim(tucurui.et)

tapajos.et[,6:7] = tapajos.et[,6:7]*100
tucurui.et[,6:7] = tucurui.et[,6:7]*100

tucurui.time = (tucurui.et$Year-1985)*12 + tucurui.et$Month + 1
dam.dist.mat = t(tucurui.xy-t(data.matrix(tucurui.et[,c("X","Y")])))
dam.dist = 1/(sqrt(rowSums(dam.dist.mat^2))+0.01)
tucurui.et$Dam.dist = dam.dist
tucurui.et$Time.diff = tucurui.time

head(tapajos.et)
head(tucurui.et)

results.et <- data.analysis(tapajos.et,tucurui.et,"Evap","Nonforest")
et.summary = hist.tapajos.model <- results.et$model1
hist.tucurui.model <- results.et$model2
et.summary = summary(hist.tapajos.model)
summary(hist.tucurui.model)
et.line = as.numeric(hist.tapajos.model$coefficients[1:2])
x.et = seq(0,100,by=0.1)
y.et = et.line[1] + et.line[2]*x
plot(x.et,y.et,type="l",xlab="Deforested (%)",ylab="Evaporation")

predict.et <- results.et$predict.tapajos(tapajos2030)
tapajos2030$Evap = predict.et

proj =  crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
tapajos30et = tapajos2030[,c("X","Y","Evap")]
coordinates(tapajos30et) <- ~X+Y
projection(tapajos30et) = proj
r = raster(res = 1)
extent(r) = extent(tapajos30et)
r.et <- rasterize(tapajos30et, r, field = "Evap", fun = mean)

plot(r.et)

#---------------Modeling for Temperature---------------------------

tapajos.temp = read.csv("tapajos-temp.csv")
tucurui.temp = read.csv("tucurui-temp.csv")

dim(tapajos.temp)
dim(tucurui.temp)

tapajos.temp[,6:7] = tapajos.temp[,6:7]*100
tucurui.temp[,6:7] = tucurui.temp[,6:7]*100

tucurui.time = (tucurui.temp$Year-1985)*12 + tucurui.temp$Month + 1
dam.dist.mat = t(tucurui.xy-t(data.matrix(tucurui.temp[,c("X","Y")])))
dam.dist = 1/(sqrt(rowSums(dam.dist.mat^2))+0.01)
tucurui.temp$Dam.dist = dam.dist
tucurui.temp$Time.diff = tucurui.time

head(tapajos.temp)
head(tucurui.temp)

results.temp <- data.analysis(tapajos.temp,tucurui.temp,"SoilTemp","Nonforest")
hist.tapajos.model <- results.temp$model1
hist.tucurui.model <- results.temp$model2
temp.summary = summary(hist.tapajos.model)
summary(hist.tucurui.model)
temp.line = as.numeric(hist.tapajos.model$coefficients[1:2])
x = seq(0,100,by=0.1)
y = temp.line[1] + temp.line[2]*x
plot(x,y,type="l",xlab="Deforested (%)",ylab="Temperature")

predict.temp <- results.temp$predict.tapajos(tapajos2030)
tapajos2030$SoilTemp = predict.temp

proj =  crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
tapajos30temp = tapajos2030[,c("X","Y","SoilTemp")]
coordinates(tapajos30temp) <- ~X+Y
projection(tapajos30temp) = proj
r = raster(res = 1)
extent(r) = extent(tapajos30temp)
r.temp <- rasterize(tapajos30temp, r, field = "SoilTemp", fun = mean)

plot(r.temp)

#-----------------Modeling for Rnet----------------------------------

tapajos.rnet = read.csv("tapajos-Rnet.csv")
tucurui.rnet = read.csv("tucurui-Rnet.csv")

dim(tapajos.rnet)
dim(tucurui.rnet)

tucurui.time = (tucurui.rnet$Year-1985)*12 + tucurui.rnet$Month + 1
dam.dist.mat = t(tucurui.xy-t(data.matrix(tucurui.rnet[,c("X","Y")])))
dam.dist = 1/(sqrt(rowSums(dam.dist.mat^2))+0.01)
tucurui.rnet$Dam.dist = dam.dist
tucurui.rnet$Time.diff = tucurui.time

results.rnet <- data.analysis(tapajos.rnet,tucurui.rnet,"Rn","Nonforest")
hist.tapajos.model <- results.rnet$model1
hist.tucurui.model <- results.rnet$model2
rnet.summary = summary(hist.tapajos.model)
summary(hist.tucurui.model)
rnet.line = as.numeric(hist.tapajos.model$coefficients[1:2])
x.rnet = seq(0,100,by=0.1)
y.rnet = rnet.line[1] + rnet.line[2]*x/100
plot(x.rnet,y.rnet,type="l",xlab="Deforested (%)",ylab="Rnet")

predict.rnet <- results.rnet$predict.tapajos(tapajos2030)
tapajos2030$Rn = predict.rnet

proj =  crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

tapajos30rnet = tapajos2030[,c("X","Y","Rn")]
coordinates(tapajos30rnet) <- ~X+Y
projection(tapajos30rnet) = proj
r = raster(res = 1)
extent(r) = extent(tapajos30rnet)
r.rnet <- rasterize(tapajos30rnet, r, field = "Rn", fun = mean)

plot(r.rnet)

par(mfrow=c(1,2))
plot(r.et)
plot(r.rnet)
#------------------------Plots---------------------------------------
months = c("Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
rnet.effects = matrix(nrow=11,ncol=3)
rnet.effects[,1] = rnet.summary$coefficients[c(3:13),1]
rnet.effects[,2] = rnet.summary$coefficients[c(3:13),1]+rnet.summary$coefficients[2,1]


et.effects = matrix(nrow=11,ncol=3)
et.effects[,1] = et.summary$coefficients[c(3:13),1]+50
et.effects[,2] = et.summary$coefficients[c(3:13),1]+et.summary$coefficients[2,1]+50


runoff.effects = matrix(nrow=11,ncol=3)
runoff.effects[,1] = runoff.summary$coefficients[c(3:13),1]+50
runoff.effects[,2] = runoff.summary$coefficients[c(3:13),1]+runoff.summary$coefficients[2,1]+50



pdf("Rnet_deforestation.pdf",width=10,height=9)
par(mar=c(6,7,5,5))
plot(rnet.effects[,1],type='l',lty=1,col='dodgerblue2',xlab='',lwd=5,ylab="Net Energy Effect",xaxt='n',ylim=c(-30,-2),main="",cex.axis=2,cex.lab=3)
axis(side = 1, at=c(1:11),labels = months,cex.axis=3,cex=2,tick=F,line=1.5)
lines(rnet.effects[,2],type='l',lty=2,col='dodgerblue2',lwd=5,ylab="Effect",xaxt='n')
legend(1,-25,c("With Deforestation","Without Deforestation"))
dev.off()


png("SeasonalEffects.png",width=800,height=800)
options(digits=7)
options(scipen=0)
par(mar=c(6,7,5,7))
plot(runoff.effects[,2]/runoff.effects[,1],type='l',lty=1,col='slateblue',xlab='',lwd=7,ylab="Discharge Effect Ratio",cex.axis=2,cex.lab=3,xaxt='n',main="")
legend(1,.995,c("Discharge","Evaporation"),col=c('slateblue','tomato3'),pch=19,cex=2.5)
par(new=TRUE)
plot(et.effects[,2]/et.effects[,1],type='l',lty=1,col='tomato3',xlab='',lwd=7,ylab="",xaxt='n',yaxt='n',main="")
axis(side = 1, at=c(1:11),labels = months,cex.axis=3,cex=2,tick=F,line=1.5)
axis(side=4,at=seq(0.9990882, 0.9995669,length.out=4),cex.axis=2)
mtext("Evaporation Effect Ratio",side=4,line=3.3,cex=3)
abline(h=1,lwd=5,lty=2)
dev.off()

pdf("Rnet_rec.pdf",width=10,height=9)
par(mar=c(6,7,5,6))
plot(x.rnet,y.rnet,type="l",xlab="Deforested (%)",ylab="Net Energy",lwd=7,cex.axis=2,cex.lab=3)
lines(x=c(20,20),y=c(0,130.9),lwd=7,col='slateblue')
lines(x=c(-5,20),y=c(130.9,130.9),lwd=7,col='slateblue')
points(x=20,y=130.9,col='slateblue',lwd=8,cex=5,pch=19)
lines(x=c(31.34,31.34),y=c(0,129.7),lwd=7,col='tomato3')
lines(x=c(-5,31.34),y=c(129.7,129.7),lwd=7,col='tomato3')
points(x=31.34,y=129.7,col='tomato3',lwd=8,cex=5,pch=19)
legend(35,133,c("% Deforested","% Legal Limit"),col=c("tomato3","slateblue"),cex=2.5,pch=20,pt.cex=4)
dev.off()



png("Rnet_presentation.png",width=800,height=800)
par(mar=c(6,7,5,6))
plot(x.rnet,y.rnet,type="l",xlab="Deforested (%)",ylab="Net Energy",lwd=7,cex.axis=2,cex.lab=3)
lines(x=c(20,20),y=c(0,130.9),lwd=7,col='slateblue')
lines(x=c(-5,20),y=c(130.9,130.9),lwd=7,col='slateblue')
points(x=20,y=130.9,col='slateblue',lwd=1,cex=5,pch=19)
lines(x=c(31.34,31.34),y=c(0,129.7),lwd=7,col='tomato3')
lines(x=c(-5,31.34),y=c(129.7,129.7),lwd=7,col='tomato3')
points(x=31.34,y=129.7,col='tomato3',lwd=1,cex=5,pch=19)
legend(35,133,c("% Deforested","% Legal Limit"),col=c("tomato3","slateblue"),cex=2.5,pch=20,pt.cex=4)
dev.off()

pdf("ET_presentation.pdf",width=10,height=9)
par(mar=c(6,7,5,6))
plot(x.et,y.et,type="l",xlab="Deforested (%)",ylab="Evapo-transpiration",lwd=7,cex.axis=2,cex.lab=3)
lines(x=c(20,20),y=c(0,129.0829),lwd=7,col='slateblue')
lines(x=c(-5,20),y=c(129.0829,129.0829),lwd=7,col='slateblue')
points(x=20,y=129.0829,col='slateblue',lwd=1,cex=5,pch=19)
lines(x=c(31.34,31.34),y=c(0,128.8641),lwd=7,col='tomato3')
lines(x=c(-5,31.34),y=c(128.8641,128.8641),lwd=7,col='tomato3')
points(x=31.34,y=128.8641,col='tomato3',lwd=1,cex=5,pch=19)
legend(35,133,c("% Deforested","% Legal Limit"),col=c("tomato3","slateblue"),cex=2.5,pch=20,pt.cex=4)
dev.off()

#to get mean deforestation in 2014 for the line plot
new.tap = tapajos.runoff[which(tapajos.runoff$Year==2014),]
new.tap = new.tap[which(new.tap$Month==12),]
mean(new.tap$Nonforest)

fut.tap = tapajos2030[which(tapajos2030$Year==2030),]
fut.tap = fut.tap[which(fut.tap$Month==1),]
mean(fut.tap$Nonforest)

#-------------Other calculations - rough-----------------------------

arima.res = auto.arima(tucurui$Rn)
ar.p = arima.res$arma[1]
ma.q = arima.res$arma[2]

model.rnet = gls(Rn ~ nonforest+time.diff+dam.dist+factor(month),
                 data = subset(tucurui,year<2014),correlation = corAR1())
testset = subset(tucurui,year==2014)
predictions <- predict(model.rnet,testset)
pred.mat = cbind(testset[,c("X","Y","year","month","Rn")],predictions)
pred.mat$error = predictions - pred.mat$Rn
sum(pred.mat$error^2)



trainset = subset(tucurui,year<2014)
testset = subset(tucurui,year==2014)
modeltrain <- data.analysis(tapajos,trainset,"Evap")
train.tucurui = modeltrain$model2
