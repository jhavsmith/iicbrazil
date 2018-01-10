
library(raster)
library(rgdal)
library(foreign)
library(RNetCDF)
library(ncdf4)
library(compare)

setwd("/home/joelsmith/Projects/brazil/data")

hist.tucurui = read.csv("formatted-historic-tucurui.csv")
hist.tapajos = read.csv("formatted-historic-tapajos.csv")
future.tapajos = read.csv("formatted-future-tapajos.csv")

hist.tucurui = hist.tucurui[complete.cases(hist.tucurui[,c("X","Y")]),]
hist.tapajos = hist.tapajos[complete.cases(hist.tapajos[,c("X","Y")]),]
future.tapajos = future.tapajos[complete.cases(future.tapajos[,c("X","Y")]),]

tapajos.grass = read.csv("tapajos-grass.csv")
tapajos.forest = read.csv("tapajos-forest.csv")

tucurui.grass = read.csv("tucurui-grass.csv")
tucurui.forest = read.csv("tucurui-forest.csv")

tucurui.xy = c(-49.646667,-3.831667)
tapajos.xy = c(-56.2845,-4.603722)


#-------------------------------------------------------------------------------
#---------------Prepare the historic data for new grids---------------------------
#------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#---------------Prepare the historic data for the Tapajos basin---------------------------
#------------------------------------------------------------------------------

new.grid.data = read.csv("newgrid-tapajos-10km.csv")
dim(new.grid.data)  
head(new.grid.data)

all.coords = unique(tapajos.forest[,c("X","Y")])
tempdata = new.grid.data
tempdata[,8:23] = t(apply(t(tempdata[,8:23]),2,cumsum))

unique.xy = unique(tempdata[,c("X","Y")])
no.of.timepoints = dim(unique(tapajos.forest[,c("year","month")]))[1]
year.month = unique(tapajos.forest[,c("year","month")])
year.col = year.month$year-2000+6
model1data = matrix(0,nrow = (nrow(tempdata)*no.of.timepoints),ncol = 23)


library(raster)
library(rgdal)
library(foreign)
library(RNetCDF)
library(ncdf4)
library(compare)

setwd("/home/joelsmith/Projects/brazil/data")

hist.tucurui = read.csv("formatted-historic-tucurui.csv")
hist.tapajos = read.csv("formatted-historic-tapajos.csv")
future.tapajos = read.csv("formatted-future-tapajos.csv")

hist.tucurui = hist.tucurui[complete.cases(hist.tucurui[,c("X","Y")]),]
hist.tapajos = hist.tapajos[complete.cases(hist.tapajos[,c("X","Y")]),]
future.tapajos = future.tapajos[complete.cases(future.tapajos[,c("X","Y")]),]

tapajos.grass = read.csv("tapajos-grass.csv")
tapajos.forest = read.csv("tapajos-forest.csv")

tucurui.grass = read.csv("tucurui-grass.csv")
tucurui.forest = read.csv("tucurui-forest.csv")

tucurui.xy = c(-49.646667,-3.831667)
tapajos.xy = c(-56.2845,-4.603722)


#-------------------------------------------------------------------------------
#---------------Prepare the historic data for new grids---------------------------
#------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#---------------Prepare the historic data for the Tapajos basin---------------------------
#------------------------------------------------------------------------------

new.grid.data = read.csv("newgrid-tapajos-10km.csv")
dim(new.grid.data)  
head(new.grid.data)

all.coords = unique(tapajos.forest[,c("X","Y")])
tempdata = new.grid.data
tempdata[,8:23] = t(apply(t(tempdata[,8:23]),2,cumsum))

unique.xy = unique(tempdata[,c("X","Y")])
no.of.timepoints = dim(unique(tapajos.forest[,c("year","month")]))[1]
year.month = unique(tapajos.forest[,c("year","month")])
year.col = year.month$year-2000+6
model1data = matrix(0,nrow = (nrow(tempdata)*no.of.timepoints),ncol = 23)

pb = txtProgressBar(1,nrow(tempdata),1,style=3)	    
for (i in 1:nrow(tempdata)){
  setTxtProgressBar(pb, i)
  xcoord = tempdata$X[i]
  ycoord = tempdata$Y[i]
  X = rep(xcoord,no.of.timepoints)
  Y = rep(ycoord,no.of.timepoints)
  state = rep(tempdata$State[i],no.of.timepoints)
  cell.data = as.numeric(tempdata[i,3:23])
  deforest = cell.data[year.col]
  nodata = rep(cell.data[1],no.of.timepoints)
  water = rep(cell.data[2],no.of.timepoints)
  nonforest = rep(cell.data[3],no.of.timepoints)+water+deforest
  cloud = rep(cell.data[4],no.of.timepoints)
  forest = rep(cell.data[5],no.of.timepoints)+cloud
  grids = small.grid(all.coords,xcoord,ycoord,0.09,0.09)
  energy = matrix(0,nrow = no.of.timepoints,ncol = 12)
  if (is.null(grids)==F){
    for (j in 1:nrow(grids)){
      idx = which((tapajos.forest$X==grids$X[j])&(tapajos.forest$Y==grids$Y[j]))
      grass.energy = as.matrix(tapajos.grass[idx,5:16])
      for.energy = as.matrix(tapajos.forest[idx,5:16])
      energy.mat = grass.energy*nonforest + for.energy*(1-nonforest-nodata)
      energy = energy + energy.mat*grids$cell.prop[j]
    }
  }
  temp = as.matrix(cbind(X,Y,state,year.month,nodata,water,cloud,deforest,forest,nonforest,energy))
  model1data[((i-1)*no.of.timepoints+1):(i*no.of.timepoints),] = temp 
}
close(pb)        
model1data = data.frame(model1data)
colnames(model1data) = c("X","Y","State","Year","Month","Nodata","Water","Cloud",
                         "Deforest","Forest","Nonforest","Evap","Precip","Qh",
                         "Qle","Rn","Qt","Qs","Qsb","wsoi","wisoi","SoilMoist","PAW")
model1data = model1data[complete.cases(model1data),]
row.names(model1data) = 1:nrow(model1data)
dim(model1data)
head(model1data)

tapajos = model1data
tapajos.et = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Evap"))
write.csv(tapajos.et,"tapajos-ET.csv",row.names = F)
tapajos.rnet = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Rn"))
write.csv(tapajos.rnet,"tapajos-Rnet.csv",row.names = F)
tapajos.runoff = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qt"))
write.csv(tapajos.runoff,"tapajos-runoff.csv",row.names = F)
tapajos.Qh = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qh"))
write.csv(tapajos.Qh,"tapajos-Qh.csv",row.names = F)


#-------------------------------------------------------------------------------
#---------------Prepare the historic data for the Tucurui basin-----------------
#-------------------------------------------------------------------------------

new.grid.data = read.csv("newgrid-tucurui-10km.csv")
dim(new.grid.data)  
head(new.grid.data)

all.coords = unique(tucurui.forest[,c("X","Y")])
tempdata = new.grid.data
tempdata[,8:23] = t(apply(t(tempdata[,8:23]),2,cumsum))

unique.xy = unique(tempdata[,c("X","Y")])
no.of.timepoints = dim(unique(tucurui.forest[,c("year","month")]))[1]
year.month = unique(tucurui.forest[,c("year","month")])
year.col = year.month$year-2000+6
model1data = matrix(0,nrow = (nrow(tempdata)*no.of.timepoints),ncol = 23)

pb = txtProgressBar(1,nrow(tempdata),1,style=3)	    
for (i in 1:nrow(tempdata)){
  setTxtProgressBar(pb, i)
  xcoord = tempdata$X[i]
  ycoord = tempdata$Y[i]
  X = rep(xcoord,no.of.timepoints)
  Y = rep(ycoord,no.of.timepoints)
  state = rep(tempdata$State[i],no.of.timepoints)
  cell.data = as.numeric(tempdata[i,3:23])
  deforest = cell.data[year.col]
  nodata = rep(cell.data[1],no.of.timepoints)
  water = rep(cell.data[2],no.of.timepoints)
  nonforest = rep(cell.data[3],no.of.timepoints)+water+deforest
  cloud = rep(cell.data[4],no.of.timepoints)
  forest = rep(cell.data[5],no.of.timepoints)+cloud
  grids = small.grid(all.coords,xcoord,ycoord,0.09,0.09)
  energy = matrix(0,nrow = no.of.timepoints,ncol = 12)
  if (is.null(grids)==F){
    for (j in 1:nrow(grids)){
      idx = which((tucurui.forest$X==grids$X[j])&(tucurui.forest$Y==grids$Y[j]))
      grass.energy = as.matrix(tucurui.grass[idx,5:16])
      for.energy = as.matrix(tucurui.forest[idx,5:16])
      energy.mat = grass.energy*nonforest + for.energy*(1-nonforest-nodata)
      energy = energy + energy.mat*grids$cell.prop[j]
    }
  }
  temp = as.matrix(cbind(X,Y,state,year.month,nodata,water,cloud,deforest,forest,nonforest,energy))
  model1data[((i-1)*no.of.timepoints+1):(i*no.of.timepoints),] = temp 
}
close(pb)
model1data = data.frame(model1data)
colnames(model1data) = c("X","Y","State","Year","Month","Nodata","Water","Cloud",
                         "Deforest","Forest","Nonforest","Evap","Precip","Qh",
                         "Qle","Rn","Qt","Qs","Qsb","wsoi","wisoi","SoilMoist","PAW")
model1data = model1data[complete.cases(model1data),]
row.names(model1data) = 1:nrow(model1data)
dim(model1data)
head(model1data)

tucurui = model1data
tucurui.et = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Evap"))
write.csv(tucurui.et,"tucurui-ET.csv",row.names = F)
tucurui.rnet = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Rn"))
write.csv(tucurui.rnet,"tucurui-Rnet.csv",row.names = F)
tucurui.runoff = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qt"))
write.csv(tucurui.runoff,"tucurui-runoff.csv",row.names = F)
tucurui.Qh = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qh"))
write.csv(tucurui.Qh,"tucurui-Qh.csv",row.names = F)


library(raster)
library(rgdal)
library(foreign)
library(RNetCDF)
library(ncdf4)
library(compare)

setwd("/home/joelsmith/Projects/brazil/data")

hist.tucurui = read.csv("formatted-historic-tucurui.csv")
hist.tapajos = read.csv("formatted-historic-tapajos.csv")
future.tapajos = read.csv("formatted-future-tapajos.csv")

hist.tucurui = hist.tucurui[complete.cases(hist.tucurui[,c("X","Y")]),]
hist.tapajos = hist.tapajos[complete.cases(hist.tapajos[,c("X","Y")]),]
future.tapajos = future.tapajos[complete.cases(future.tapajos[,c("X","Y")]),]

tapajos.grass = read.csv("tapajos-grass.csv")
tapajos.forest = read.csv("tapajos-forest.csv")

tucurui.grass = read.csv("tucurui-grass.csv")
tucurui.forest = read.csv("tucurui-forest.csv")

tucurui.xy = c(-49.646667,-3.831667)
tapajos.xy = c(-56.2845,-4.603722)


#-------------------------------------------------------------------------------
#---------------Prepare the historic data for new grids---------------------------
#------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#---------------Prepare the historic data for the Tapajos basin---------------------------
#------------------------------------------------------------------------------

new.grid.data = read.csv("newgrid-tapajos-10km.csv")
dim(new.grid.data)  
head(new.grid.data)

all.coords = unique(tapajos.forest[,c("X","Y")])
tempdata = new.grid.data
tempdata[,8:23] = t(apply(t(tempdata[,8:23]),2,cumsum))

unique.xy = unique(tempdata[,c("X","Y")])
no.of.timepoints = dim(unique(tapajos.forest[,c("year","month")]))[1]
year.month = unique(tapajos.forest[,c("year","month")])
year.col = year.month$year-2000+6
model1data = matrix(0,nrow = (nrow(tempdata)*no.of.timepoints),ncol = 23)

for (i in 1:nrow(tempdata)){
  xcoord = tempdata$X[i]
  ycoord = tempdata$Y[i]
  X = rep(xcoord,no.of.timepoints)
  Y = rep(ycoord,no.of.timepoints)
  state = rep(tempdata$State[i],no.of.timepoints)
  cell.data = as.numeric(tempdata[i,3:23])
  deforest = cell.data[year.col]
  nodata = rep(cell.data[1],no.of.timepoints)
  water = rep(cell.data[2],no.of.timepoints)
  nonforest = rep(cell.data[3],no.of.timepoints)+water+deforest
  cloud = rep(cell.data[4],no.of.timepoints)
  forest = rep(cell.data[5],no.of.timepoints)+cloud
  grids = small.grid(all.coords,xcoord,ycoord,0.09,0.09)
  energy = matrix(0,nrow = no.of.timepoints,ncol = 12)
  if (is.null(grids)==F){
    for (j in 1:nrow(grids)){
      idx = which((tapajos.forest$X==grids$X[j])&(tapajos.forest$Y==grids$Y[j]))
      grass.energy = as.matrix(tapajos.grass[idx,5:16])
      for.energy = as.matrix(tapajos.forest[idx,5:16])
      energy.mat = grass.energy*nonforest + for.energy*(1-nonforest-nodata)
      energy = energy + energy.mat*grids$cell.prop[j]
    }
  }
  temp = as.matrix(cbind(X,Y,state,year.month,nodata,water,cloud,deforest,forest,nonforest,energy))
  model1data[((i-1)*no.of.timepoints+1):(i*no.of.timepoints),] = temp 
}

model1data = data.frame(model1data)
colnames(model1data) = c("X","Y","State","Year","Month","Nodata","Water","Cloud",
                         "Deforest","Forest","Nonforest","Evap","Precip","Qh",
                         "Qle","Rn","Qt","Qs","Qsb","wsoi","wisoi","SoilMoist","PAW")
model1data = model1data[complete.cases(model1data),]
row.names(model1data) = 1:nrow(model1data)
dim(model1data)
head(model1data)

tapajos = model1data
tapajos.et = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Evap"))
write.csv(tapajos.et,"tapajos-ET.csv",row.names = F)
tapajos.rnet = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Rn"))
write.csv(tapajos.rnet,"tapajos-Rnet.csv",row.names = F)
tapajos.runoff = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qt"))
write.csv(tapajos.runoff,"tapajos-runoff.csv",row.names = F)
tapajos.Qh = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qh"))
write.csv(tapajos.Qh,"tapajos-Qh.csv",row.names = F)


#-------------------------------------------------------------------------------
#---------------Prepare the historic data for the Tucurui basin---------------------------
#------------------------------------------------------------------------------

new.grid.data = read.csv("newgrid-tucurui-10km.csv")
dim(new.grid.data)  
head(new.grid.data)

all.coords = unique(tucurui.forest[,c("X","Y")])
tempdata = new.grid.data
tempdata[,8:23] = t(apply(t(tempdata[,8:23]),2,cumsum))

unique.xy = unique(tempdata[,c("X","Y")])
no.of.timepoints = dim(unique(tucurui.forest[,c("year","month")]))[1]
year.month = unique(tucurui.forest[,c("year","month")])
year.col = year.month$year-2000+6
model1data = matrix(0,nrow = (nrow(tempdata)*no.of.timepoints),ncol = 23)

for (i in 1:nrow(tempdata)){
  xcoord = tempdata$X[i]
  ycoord = tempdata$Y[i]
  X = rep(xcoord,no.of.timepoints)
  Y = rep(ycoord,no.of.timepoints)
  state = rep(tempdata$State[i],no.of.timepoints)
  cell.data = as.numeric(tempdata[i,3:23])
  deforest = cell.data[year.col]
  nodata = rep(cell.data[1],no.of.timepoints)
  water = rep(cell.data[2],no.of.timepoints)
  nonforest = rep(cell.data[3],no.of.timepoints)+water+deforest
  cloud = rep(cell.data[4],no.of.timepoints)
  forest = rep(cell.data[5],no.of.timepoints)+cloud
  grids = small.grid(all.coords,xcoord,ycoord,0.09,0.09)
  energy = matrix(0,nrow = no.of.timepoints,ncol = 12)
  if (is.null(grids)==F){
    for (j in 1:nrow(grids)){
      idx = which((tucurui.forest$X==grids$X[j])&(tucurui.forest$Y==grids$Y[j]))
      grass.energy = as.matrix(tucurui.grass[idx,5:16])
      for.energy = as.matrix(tucurui.forest[idx,5:16])
      energy.mat = grass.energy*nonforest + for.energy*(1-nonforest-nodata)
      energy = energy + energy.mat*grids$cell.prop[j]
    }
  }
  temp = as.matrix(cbind(X,Y,state,year.month,nodata,water,cloud,deforest,forest,nonforest,energy))
  model1data[((i-1)*no.of.timepoints+1):(i*no.of.timepoints),] = temp 
}

model1data = data.frame(model1data)
colnames(model1data) = c("X","Y","State","Year","Month","Nodata","Water","Cloud",
                         "Deforest","Forest","Nonforest","Evap","Precip","Qh",
                         "Qle","Rn","Qt","Qs","Qsb","wsoi","wisoi","SoilMoist","PAW")
model1data = model1data[complete.cases(model1data),]
row.names(model1data) = 1:nrow(model1data)
dim(model1data)
head(model1data)

tucurui = model1data
tucurui.et = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Evap"))
write.csv(tucurui.et,"tucurui-ET.csv",row.names = F)
tucurui.rnet = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Rn"))
write.csv(tucurui.rnet,"tucurui-Rnet.csv",row.names = F)
tucurui.runoff = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qt"))
write.csv(tucurui.runoff,"tucurui-runoff.csv",row.names = F)
tucurui.Qh = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qh"))
write.csv(tucurui.Qh,"tucurui-Qh.csv",row.names = F)


for (i in 1:nrow(tempdata)){
  xcoord = tempdata$X[i]
  ycoord = tempdata$Y[i]
  X = rep(xcoord,no.of.timepoints)
  Y = rep(ycoord,no.of.timepoints)
  state = rep(tempdata$State[i],no.of.timepoints)
  cell.data = as.numeric(tempdata[i,3:23])
  deforest = cell.data[year.col]
  nodata = rep(cell.data[1],no.of.timepoints)
  water = rep(cell.data[2],no.of.timepoints)
  nonforest = rep(cell.data[3],no.of.timepoints)+water+deforest
  cloud = rep(cell.data[4],no.of.timepoints)
  forest = rep(cell.data[5],no.of.timepoints)+cloud
  grids = small.grid(all.coords,xcoord,ycoord,0.09,0.09)
  energy = matrix(0,nrow = no.of.timepoints,ncol = 12)
  if (is.null(grids)==F){
    for (j in 1:nrow(grids)){
      idx = which((tapajos.forest$X==grids$X[j])&(tapajos.forest$Y==grids$Y[j]))
      grass.energy = as.matrix(tapajos.grass[idx,5:16])
      for.energy = as.matrix(tapajos.forest[idx,5:16])
      energy.mat = grass.energy*nonforest + for.energy*(1-nonforest-nodata)
      energy = energy + energy.mat*grids$cell.prop[j]
    }
  }
  temp = as.matrix(cbind(X,Y,state,year.month,nodata,water,cloud,deforest,forest,nonforest,energy))
  model1data[((i-1)*no.of.timepoints+1):(i*no.of.timepoints),] = temp 
}

model1data = data.frame(model1data)
colnames(model1data) = c("X","Y","State","Year","Month","Nodata","Water","Cloud",
                         "Deforest","Forest","Nonforest","Evap","Precip","Qh",
                         "Qle","Rn","Qt","Qs","Qsb","wsoi","wisoi","SoilMoist","PAW")
model1data = model1data[complete.cases(model1data),]
row.names(model1data) = 1:nrow(model1data)
dim(model1data)
head(model1data)

tapajos = model1data
tapajos.et = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Evap"))
write.csv(tapajos.et,"tapajos-ET.csv",row.names = F)
tapajos.rnet = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Rn"))
write.csv(tapajos.rnet,"tapajos-Rnet.csv",row.names = F)
tapajos.runoff = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qt"))
write.csv(tapajos.runoff,"tapajos-runoff.csv",row.names = F)
tapajos.Qh = subset(tapajos,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qh"))
write.csv(tapajos.Qh,"tapajos-Qh.csv",row.names = F)


#-------------------------------------------------------------------------------
#---------------Prepare the historic data for the Tucurui basin---------------------------
#------------------------------------------------------------------------------

new.grid.data = read.csv("newgrid-tucurui-10km.csv")
dim(new.grid.data)  
head(new.grid.data)

all.coords = unique(tucurui.forest[,c("X","Y")])
tempdata = new.grid.data
tempdata[,8:23] = t(apply(t(tempdata[,8:23]),2,cumsum))

unique.xy = unique(tempdata[,c("X","Y")])
no.of.timepoints = dim(unique(tucurui.forest[,c("year","month")]))[1]
year.month = unique(tucurui.forest[,c("year","month")])
year.col = year.month$year-2000+6
model1data = matrix(0,nrow = (nrow(tempdata)*no.of.timepoints),ncol = 23)

for (i in 1:nrow(tempdata)){
  xcoord = tempdata$X[i]
  ycoord = tempdata$Y[i]
  X = rep(xcoord,no.of.timepoints)
  Y = rep(ycoord,no.of.timepoints)
  state = rep(tempdata$State[i],no.of.timepoints)
  cell.data = as.numeric(tempdata[i,3:23])
  deforest = cell.data[year.col]
  nodata = rep(cell.data[1],no.of.timepoints)
  water = rep(cell.data[2],no.of.timepoints)
  nonforest = rep(cell.data[3],no.of.timepoints)+water+deforest
  cloud = rep(cell.data[4],no.of.timepoints)
  forest = rep(cell.data[5],no.of.timepoints)+cloud
  grids = small.grid(all.coords,xcoord,ycoord,0.09,0.09)
  energy = matrix(0,nrow = no.of.timepoints,ncol = 12)
  if (is.null(grids)==F){
    for (j in 1:nrow(grids)){
      idx = which((tucurui.forest$X==grids$X[j])&(tucurui.forest$Y==grids$Y[j]))
      grass.energy = as.matrix(tucurui.grass[idx,5:16])
      for.energy = as.matrix(tucurui.forest[idx,5:16])
      energy.mat = grass.energy*nonforest + for.energy*(1-nonforest-nodata)
      energy = energy + energy.mat*grids$cell.prop[j]
    }
  }
  temp = as.matrix(cbind(X,Y,state,year.month,nodata,water,cloud,deforest,forest,nonforest,energy))
  model1data[((i-1)*no.of.timepoints+1):(i*no.of.timepoints),] = temp 
}

model1data = data.frame(model1data)
colnames(model1data) = c("X","Y","State","Year","Month","Nodata","Water","Cloud",
                         "Deforest","Forest","Nonforest","Evap","Precip","Qh",
                         "Qle","Rn","Qt","Qs","Qsb","wsoi","wisoi","SoilMoist","PAW")
model1data = model1data[complete.cases(model1data),]
row.names(model1data) = 1:nrow(model1data)
dim(model1data)
head(model1data)

tucurui = model1data
tucurui.et = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Evap"))
write.csv(tucurui.et,"tucurui-ET.csv",row.names = F)
tucurui.rnet = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Rn"))
write.csv(tucurui.rnet,"tucurui-Rnet.csv",row.names = F)
tucurui.runoff = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qt"))
write.csv(tucurui.runoff,"tucurui-runoff.csv",row.names = F)
tucurui.Qh = subset(tucurui,select=c("X","Y","State","Year","Month","Forest","Nonforest","Qh"))
write.csv(tucurui.Qh,"tucurui-Qh.csv",row.names = F)

#----------------------------------------------------------------
#--------------------Prepare the new grids-----------------------
#----------------------------------------------------------------

xlimits = range(hist.tapajos$X)
ylimits = range(hist.tapajos$Y)

new.grid.x = seq(xlimits[1],xlimits[2],by=0.09)
new.grid.y = seq(ylimits[1],ylimits[2],by=0.09)

new.grid.xy = expand.grid(new.grid.x,new.grid.y)
colnames(new.grid.xy) = c("X","Y")

tempdata = hist.tapajos[,6:29]
tempdata[,9] = tempdata[,8] + tempdata[,9]
colnames(tempdata)[9] = "Def00"
tempdata = tempdata[,-8]
dim(tempdata)
head(tempdata)

namevector = colnames(tempdata)[3:23]
new.grid.data = new.grid.xy
new.grid.data[,namevector] = 0
N = nrow(new.grid.data)
new.grid.data$State = NA

for (i in 1:N){
  xcoord = new.grid.data$X[i]
  xid = which((hist.tapajos$X>(xcoord-0.045)) & (hist.tapajos$X<(xcoord+0.045)))
  ycoord = new.grid.data$Y[i]
  yid = which((hist.tapajos$Y>(ycoord-0.045)) & (hist.tapajos$Y<(ycoord+0.045)))
  idx = intersect(xid,yid)
  if (length(idx)>0){
    new.grid.data[i,3:23] = colSums(tempdata[idx,3:23])
    new.grid.data$State[i] = as.character(hist.tapajos$Initial[idx[1]])
  }
}

dim(new.grid.data)
idx = which(is.na(new.grid.data$State)==T)
newgrid = new.grid.data[-idx,]
rownames(newgrid) = 1:nrow(newgrid)
newgrid[,3:23] = newgrid[,3:23]/(10^8)
idx = which(rowSums(newgrid[,3:23])>1)
newgrid[idx,3:23] = newgrid[idx,3:23]/rowSums(newgrid[idx,3:23])
head(newgrid)
write.csv(newgrid,"newgrid-tapajos-10km.csv",row.names = F)



#------------------Tucurui data - new grid--------------------------------

xlimits = range(hist.tucurui$X)
ylimits = range(hist.tucurui$Y)

new.grid.x = seq(xlimits[1],xlimits[2],by=0.09)
new.grid.y = seq(ylimits[1],ylimits[2],by=0.09)

new.grid.xy = expand.grid(new.grid.x,new.grid.y)
colnames(new.grid.xy) = c("X","Y")

tempdata = hist.tucurui[,6:29]
tempdata[,9] = tempdata[,8] + tempdata[,9]
colnames(tempdata)[9] = "Def00"
tempdata = tempdata[,-8]
dim(tempdata)
head(tempdata)

namevector = colnames(tempdata)[3:23]
new.grid.data = new.grid.xy
new.grid.data[,namevector] = 0
N = nrow(new.grid.data)
new.grid.data$State = NA

for (i in 1:N){
  xcoord = new.grid.data$X[i]
  xid = which((hist.tucurui$X>(xcoord-0.045)) & (hist.tucurui$X<(xcoord+0.045)))
  ycoord = new.grid.data$Y[i]
  yid = which((hist.tucurui$Y>(ycoord-0.045)) & (hist.tucurui$Y<(ycoord+0.045)))
  idx = intersect(xid,yid)
  if (length(idx)>0){
    new.grid.data[i,3:23] = colSums(tempdata[idx,3:23])
    new.grid.data$State[i] = as.character(hist.tucurui$Initial[idx[1]])
  }
}

dim(new.grid.data)
idx = which(is.na(new.grid.data$State)==T)
newgrid = new.grid.data[-idx,]
rownames(newgrid) = 1:nrow(newgrid)
newgrid[,3:23] = newgrid[,3:23]/(10^8)
idx = which(rowSums(newgrid[,3:23])>1)
newgrid[idx,3:23] = newgrid[idx,3:23]/rowSums(newgrid[idx,3:23])
head(newgrid)
write.csv(newgrid,"newgrid-tucurui-10km.csv",row.names = F)

#-------------------------------------------------------------------------------
#------------------------------IBIS Model data-------------------------------------
#------------------------------------------------------------------------------

setwd("C:/Academics/Chicago/ThirdYear/IIC/Dam-Project/Data/IBIS_data/IBIS_Newdata/IBIS_PotVeg/")

xlimits = range(hist.tapajos$X)
ylimits = range(hist.tapajos$Y)

var.name = c("Evap","Precip","Qh","Qle","Rn","runoff","SoilMoist")
general.info = c("latitude","longitude","time","time_weights","date")
fname = paste(var.name[1],".nc",sep='')
fid = open.nc(fname)
time.index = make.time.index(1960,2015)
required.time = make.time.index(2000,2015)

time = var.get.nc(fid,"time",start=1,count=NA)
ylat = var.get.nc(fid,"latitude",start=1,count=NA)
xlon = var.get.nc(fid,"longitude",start=1,count=NA)
coords = expand.grid(xlon,ylat)

no.of.timepoints = length(required.time)
no.of.obs = dim(coords)[1]
fulldata.forest = matrix(nrow = 0,ncol = 0)

for (i in 1:no.of.timepoints){
  time.year = as.numeric(substr(required.time[i],1,4))
  time.month = as.numeric(substr(required.time[i],6,7))
  temp = c()
  tempname = c()
  for (j in 1:length(var.name)){
    fname = paste(var.name[j],".nc",sep='')
    fid = open.nc(fname)
    nc = nc_open(fname)
    variables = names(nc[['var']])
    idx = variables[!(variables %in% general.info)]
    for (k in 1:length(idx)){
      var.values = var.get.nc(fid,idx[k],start=c(1,1,which(time.index==required.time[i])),count=c(NA,NA,1))
      temp = cbind(temp,as.vector(var.values))
    }
    tempname = c(tempname,idx)
    close.nc(fid)
    nc_close(nc)
  }
  year = rep(time.year,no.of.obs)
  month = rep(time.month,no.of.obs)
  newdata = cbind(year,month,coords,temp)
  fulldata.forest = rbind(fulldata.forest,newdata)
}
colnames(fulldata.forest) <- c("year","month","X","Y",tempname)

temp1 = which(fulldata.forest$X>=(xlimits[1]-0.25) & fulldata.forest$X<=(xlimits[2]+0.25))
temp2 = which(fulldata.forest$Y>=(ylimits[1]-0.25) & fulldata.forest$Y<=(ylimits[2]+0.25))
temp = intersect(temp1,temp2)

tapajos.forest = estimate.missing.data(fulldata.forest[temp,])


setwd("C:/Academics/Chicago/ThirdYear/IIC/Dam-Project/Data/IBIS_data/IBIS_Newdata/IBIS_Grass/")


fulldata.grass = matrix(nrow = 0,ncol = 0)

for (i in 1:no.of.timepoints){
  time.year = as.numeric(substr(required.time[i],1,4))
  time.month = as.numeric(substr(required.time[i],6,7))
  temp = c()
  tempname = c()
  for (j in 1:length(var.name)){
    fname = paste(var.name[j],".nc",sep='')
    fid = open.nc(fname)
    nc = nc_open(fname)
    variables = names(nc[['var']])
    idx = variables[!(variables %in% general.info)]
    for (k in 1:length(idx)){
      var.values = var.get.nc(fid,idx[k],start=c(1,1,which(time.index==required.time[i])),count=c(NA,NA,1))
      temp = cbind(temp,as.vector(var.values))
    }
    tempname = c(tempname,idx)
    close.nc(fid)
    nc_close(nc)
  }
  year = rep(time.year,no.of.obs)
  month = rep(time.month,no.of.obs)
  newdata = cbind(year,month,coords,temp)
  fulldata.grass = rbind(fulldata.grass,newdata)
}
colnames(fulldata.grass) <- c("year","month","X","Y",tempname)

temp1 = which(fulldata.grass$X>=(xlimits[1]-0.25) & fulldata.grass$X<=(xlimits[2]+0.25))
temp2 = which(fulldata.grass$Y>=(ylimits[1]-0.25) & fulldata.grass$Y<=(ylimits[2]+0.25))
temp = intersect(temp1,temp2)

tapajos.grass = estimate.missing.data(fulldata.grass[temp,])

setwd("C:/Academics/Chicago/ThirdYear/IIC/Dam-Project/Data/Code/")

write.csv(tapajos.grass,file = "tapajos-grass.csv",row.names = F)
write.csv(tapajos.forest,file = "tapajos-forest.csv",row.names = F)


#-------------------------------------------------------------------------------
#------------------------------IBIS Model data-------------------------------------
#------------------------------------------------------------------------------

xlimits = range(hist.tucurui$X)
ylimits = range(hist.tucurui$Y)
  
temp1 = which(fulldata.forest$X>=(xlimits[1]-0.25) & fulldata.forest$X<=(xlimits[2]+0.25))
temp2 = which(fulldata.forest$Y>=(ylimits[1]-0.25) & fulldata.forest$Y<=(ylimits[2]+0.25))
temp = intersect(temp1,temp2)

tucurui.forest = estimate.missing.data(fulldata.forest[temp,])


temp1 = which(fulldata.grass$X>=(xlimits[1]-0.25) & fulldata.grass$X<=(xlimits[2]+0.25))
temp2 = which(fulldata.grass$Y>=(ylimits[1]-0.25) & fulldata.grass$Y<=(ylimits[2]+0.25))
temp = intersect(temp1,temp2)

tucurui.grass = estimate.missing.data(fulldata.grass[temp,])

dim(tucurui.grass)
dim(tucurui.forest)

setwd("C:/Academics/Chicago/ThirdYear/IIC/Dam-Project/Data/Code/")

write.csv(tucurui.grass,file = "tucurui-grass.csv",row.names = F)
write.csv(tucurui.forest,file = "tucurui-forest.csv",row.names = F)

#---------------------------------------------------------------------
#-----------------------------Rough-------------------------------------
#-------------------------------------------------------------------------

tempdata = hist.tapajos[,c(6:12,14:29)]
names(tempdata)[8] = "Def00"
tempdata$Def00 = hist.tapajos$Def97 + hist.tapajos$Def98_00

proj =  crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

tapajos = tempdata[,c("X","Y","Forest")]
coordinates(tapajos) <- ~X+Y
#raster.tapajos <- rasterFromXYZ(tapajos)
projection(tapajos) = proj

r = raster(res = 1)
extent(r) = extent(tapajos)
r.00 <- rasterize(tapajos, r, field = 1, fun = sum)
r.00[is.na(r.00)] = 0

plot(r.00)


hist(r.00)
View(getValues(r.00))
#rasterVis::levelplot(r.00)

for (i in 1:nrow(unique.xy)){
  xcoord = unique.xy$X[i]
  ycoord = unique.xy$Y[i]
  idx = which((tapajos.forest$X==xcoord) & (tapajos.forest$Y==ycoord))
  grass.energy = as.matrix(tapajos.grass[idx,5:16])
  for.energy = as.matrix(tapajos.forest[idx,5:16])
  grids = big.grids(all.coords,xcoord,ycoord)
  X = rep(xcoord,no.of.timepoints)
  Y = rep(ycoord,no.of.timepoints)
  year = tapajos.forest$year[idx]
  month = tapajos.forest$month[idx]
  year.def.col = 6+year-2000
  deforest = numeric(no.of.timepoints)
  nonforest = numeric(no.of.timepoints)
  forest = numeric(no.of.timepoints)
  cloud = numeric(no.of.timepoints)
  water = numeric(no.of.timepoints)
  nodata = numeric(no.of.timepoints)
  energy = matrix(nrow = no.of.timepoints,ncol = 12)
  if (is.null(grids)==FALSE){
    cell.id = which(((tempdata$X %in% grids$X)==1) & ((tempdata$Y %in% grids$Y)==1))
    cell.data = colMeans(tempdata[cell.id,3:23])/1000000
    deforest = cell.data[year.def.col]
    nodata = rep(cell.data[1],no.of.timepoints)
    water = rep(cell.data[2],no.of.timepoints)
    nonforest = rep(cell.data[3],no.of.timepoints)+water+deforest
    cloud = rep(cell.data[4],no.of.timepoints)
    forest = rep(cell.data[5],no.of.timepoints)
    energy = grass.energy*nonforest + for.energy*(1-nonforest-nodata)
  }
  temp = cbind(X,Y,year,month,nodata,water,cloud,deforest,forest,nonforest,energy)
  model1data = rbind(model1data,temp)
}

row.names(model1data) = 1:nrow(model1data)
model1data = data.frame(model1data)
colnames(model1data)[11:22] = c("Evap","Precip","Qh","Qle","Rn","Qt","Qs","Qsb",
                                "wsoi","wisoi","SoilMoist","PAW")
model1data = model1data[complete.cases(model1data),]
row.names(model1data) = 1:nrow(model1data)
dim(model1data)
head(model1data)

tapajos = model1data

#-------------------------------------------------------------------------------
#---------------Prepare the historic data for the Tucurui basin---------------------------
#------------------------------------------------------------------------------

N = dim(hist.tucurui)[1]

ibisxy = tucurui.grass[,c("X","Y")]

all.coords = unique(hist.tucurui[,c("X","Y")])
tempdata = hist.tucurui[,c(1:2,7:10,12:27)]
names(tempdata)[7] = "def00"
tempdata$def00 = hist.tucurui$def97 + hist.tucurui$def98_00

nonmissing = which(rowSums(tempdata[,3:22])>0)
tempdata = tempdata[nonmissing,]

unique.xy = unique(tucurui.forest[,c("X","Y")])
no.of.timepoints = dim(unique(tucurui.forest[,c("year","month")]))[1]
model1data = c()

for (i in 1:dim(unique.xy)[1]){
  xcoord = unique.xy$X[i]
  ycoord = unique.xy$Y[i]
  idx = which((tucurui.forest$X==xcoord) & (tucurui.forest$Y==ycoord))
  grass.energy = as.matrix(tucurui.grass[idx,5:16])
  for.energy = as.matrix(tucurui.forest[idx,5:16])
  grids = big.grids(all.coords,xcoord,ycoord)
  X = rep(xcoord,no.of.timepoints)
  Y = rep(ycoord,no.of.timepoints)
  year = tucurui.forest$year[idx]
  month = tucurui.forest$month[idx]
  year.def.col = year-1999
  deforest = numeric(no.of.timepoints)
  nonforest = numeric(no.of.timepoints)
  forest = numeric(no.of.timepoints)
  cloud = numeric(no.of.timepoints)
  water = numeric(no.of.timepoints)
  energy = matrix(nrow = no.of.timepoints,ncol = 12)
  if (is.null(grids)==FALSE){
    cell.id = which(((tempdata$X %in% grids$X)==1) & ((tempdata$Y %in% grids$Y)==1))
    cell.data = colMeans(tempdata[cell.id,3:22])/sum(colMeans(tempdata[cell.id,3:22]))
    deforest = cumsum(cell.data[5:20])[year.def.col]
    water = rep(cell.data[1],no.of.timepoints)
    nonforest = rep(cell.data[2],no.of.timepoints)+water+deforest
    cloud = rep(cell.data[3],no.of.timepoints)
    forest = rep(cell.data[4],no.of.timepoints)
    energy = grass.energy*nonforest + for.energy*(1-nonforest-nodata)
  }
  temp = cbind(X,Y,year,month,water,cloud,deforest,forest,nonforest,energy)
  model1data = rbind(model1data,temp)
}

row.names(model1data) = 1:nrow(model1data)
model1data = data.frame(model1data)
colnames(model1data)[10:21] = c("Evap","Precip","Qh","Qle","Rn","Qt","Qs","Qsb",
                                "wsoi","wisoi","SoilMoist","PAW")
model1data = model1data[complete.cases(model1data),]
row.names(model1data) = 1:nrow(model1data)
dim(model1data)
head(model1data)

tucurui = model1data


