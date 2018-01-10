

#-----------------------------------------------
#------------------Make Maps--------------------
#-----------------------------------------------

library(reshape2)
library(ggplot2)
library(ggmap)
library(RgoogleMaps)
library(rgdal)
library(broom)
library(gridExtra)

#--------------Deforestation in 2000, 2015, 2030-----------------

tapajos.munic = readOGR(".","tapajos_all")
tapajos.munic = spTransform(tapajos.munic, CRS("+datum=WGS84 +proj=longlat"))
tapajos.munic =  fortify(tapajos.munic)

mean.lon=mean(hist.tapajos$X)
mean.lat=mean(hist.tapajos$Y)
map <- get_googlemap(center=c(lon=mean.lon,lat=mean.lat),zoom=5, maptype='hybrid',path = "&style=feature:all|element:labels|visibility:off")

hist.data = hist.tapajos[which(hist.tapajos$Def15!=0),]
hist.data$Def15 = hist.data$Def15/1000000

plot2 = ggmap(map, extent = "normal",base_layer=ggplot(aes(x=hist.data$X,y=hist.data$Y),data=hist.data),maprange=FALSE) +
geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.2,color='grey', data=tapajos.munic, alpha=.2) +
theme(legend.position = "none",
plot.title = element_text(size = rel(3)),
axis.title.x=element_blank(),
axis.text.x=element_blank(),
axis.ticks.x=element_blank(),
axis.title.y=element_blank(),
axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.line = element_line(colour = "black"))

future.data = future.tapajos[which(future.tapajos$Def_2030!=0),]
future.data$Def_2030 = future.data$Def_2030/1000000

plot3 = ggmap(map, extent = "normal",base_layer=ggplot(aes(x=future.data$X,y=future.data$Y),data=future.data),maprange=FALSE) +
labs(title = "2030") +
geom_point(aes(x = future.data$X, y = future.data$Y, colour = Def_2030), data=future.data, alpha = 1,shape=15,size=.1) +
scale_color_gradient(low='orange',high='red') +
geom_polygon(aes(x=long, y=lat, group=group), fill='grey', size=.2,color='grey', data=tapajos.munic, alpha=.2) +
theme(legend.position = "none",
plot.title = element_text(size = rel(3)),
axis.title.x=element_blank(),
axis.text.x=element_blank(),
axis.ticks.x=element_blank(),
axis.title.y=element_blank(),
axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.line = element_line(colour = "black"))

png("deforest_compare4.png",width = 900, height = 600)
grid.arrange(plot2,plot3,ncol=2)
dev.off()


png("tapajos.region.png",width = 900, height = 900)
plot2
dev.off()

library(RNetCDF)
#Pick which variable
var.name = "SoilTemp"
#Open data file
fname = paste("~/Projects/brazil/data/",var.name,".nc",sep='')
fid = open.nc(fname)
#Print dataset summary information
print.nc(fid)
#Express time in years / months
time.index = make.time.index(1960,2015)
#Choose a time (year,month)
my.time = c(1960,3)
#Read data
var.values = var.get.nc(fid,var.name,start=c(1,1,which(time.index==paste(my.time[1],my.time[2]))),count=c(NA,NA,1))
time = var.get.nc(fid,"time",start=1,count=NA)
ylat = var.get.nc(fid,"latitude",start=1,count=NA)
xlon = var.get.nc(fid,"longitude",start=1,count=NA)
#Close file
close.nc(fid)









