#####################################################################################
###############################         READ ME         #############################
#####################################################################################
# The following plotting code uses the Google Maps Platform. Google recently 
# tightened restrictions for access to the database so it is now necessary to register
# with Google to obtain an API key in order ggmap() to work. This is a free service, 
# although does require the input of a credit card (supposedly no charges will be made 
# without authorization from the user). To obtain an API key, register with Google here:
# https://cloud.google.com/maps-platform/#get-started
#
# Then, install ggmap() from github
# devtools::install_github("dkahle/ggmap")
#
# Then, save your API key info within R
# ggmap::register_google(key = "copy your API key here")

#####################################################################################
###############################  Code to Make Figures   #############################
#####################################################################################

require(ggplot2)
require(ggmap)
require(RANN)
require(chron)
require(fields)
theme_set(theme_bw(base_size=18))

load("pred_block.Rda")
load("data_blockmed.Rda")
df$Y_block_med <- df$Y
df$Car = ifelse(df$Car=="Car_B",2,1)

p0 = get_map(location = c(lon = -122.2505, lat = mean(pred_block$Latitude)),
             zoom = 12, maptype = 'toner-lite', source="stamen")

# load data ---------------------------------------------------------------
covs_name <- "Data_Covariates.csv" #GIS covariates
goog_name <- "oakland_data_simulated.csv"

covs      <- read.csv(covs_name,nrows=Inf)
goog      <- read.csv(goog_name,nrows=Inf,as.is=TRUE)
rm(segs_name,covs_name,goog_name)

goog      <- goog[!is.na(goog[,6]),] # work with NO2 data. rm NAs
s         <- goog[,4:3]
nn        <- as.vector(nn2(covs[,2:3],s,k=1)$nn.idx)
spID <- covs$ID[nn]
Y         <- goog[,6]
Y         <- log(Y)

pb.txt <- goog[,1]
pb.date <- as.POSIXct(pb.txt, origin = "2015-01-01 00:00:00",tz="America/Los_Angeles")
attributes(pb.date)$tzone <- "America/Los_Angeles"
t <- pb.date

car       <- goog[,2]
s         <- goog[,4:3]
year      <- as.numeric(substr(t,1,4))
month     <- as.numeric(substr(t,6,7))
day       <- as.numeric(substr(t,9,10))
hour      <- as.numeric(substr(t,12,13))
min       <- as.numeric(substr(t,15,16))
sec       <- as.numeric(substr(t,18,19))

dd <- julian(month, day, year, origin=c(1,1,2015))
dt <- hour/24+min/(60*24)+sec/(60*60*24)

############################### Figure 1 Drive Maps ############################

sdf <- data.frame(x=covs[,2],y=covs[,3]) # covs [,2:3] are the coordinates corresponding to locID
scardf <- data.frame(lon=s[,1],lat=s[,2],car=car,time=t)
Ainit = min(which(car=="Car_A")) #day that car A starts
idx = which(day==day[Ainit]&month==month[Ainit]&year==year[Ainit])

df_217 <- scardf[dd==217,]
ggmap(p0) +
  geom_point(data=df_217,aes(x=lon,y=lat,colour=car),size=0.5) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  ggtitle("Roads Driven on August 6, 2015") + scale_colour_discrete(name="") +
  guides(colour = guide_legend(override.aes = list(size=2)))

df_490 <- scardf[dd==490,]
ggmap(p0) +
  geom_point(data=df_490,aes(x=lon,y=lat,colour=car),size=0.5) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  ggtitle("Roads Driven on May 5, 2016") + scale_colour_discrete(name="") +
  guides(colour = guide_legend(override.aes = list(size=2)))

############################### Figure 2 Temporal Aggregation ##############################

day <- 454
hr <- 15
kp <- (dd == day & hour == hr)

rawY <- Y[kp]
tY <- dt[kp] + dd[kp]
car_kp <- car[kp]
dY <- t[kp]
sp_kp <- s[kp,]

sub_1min <- subset(df_block_tseg1min, (floor(t) == 454) & (hour == 15))
sub_15sec <- subset(df_block_tseg15sec, floor(t) == 454 & hour == 15)

a <- c(1, 357, 779, 1286,  1773, 2224)
b <- c("15:00", "15:10", "15:20", "15:30", "15:40", "15:50")

qplot(x=tY, y=rawY, xlab="Time", ylab="log(NO2)",main = "March 3, 2016, 15:00", colour = "Raw Data", size=I(1.5), alpha=I(0.5)) + 
  geom_point(aes(x=sub_15sec$tseg15sec, y=sub_15sec$Y_block_med, col="15sec Median"), size=2.5, alpha=1) + 
  geom_point(aes(x=sub_1min$tseg1min, y=sub_1min$Y_block_med, col="1min Median"), size=2.5, alpha=1) + 
  scale_x_continuous(breaks=tY[a], labels=b) + theme(plot.title = element_text(hjust = 0.5), legend.position="bottom") +
  scale_colour_manual("", 
                      breaks = c("Raw Data", "15sec Median", "1min Median"),
                      values = c("red" ,"green2", "black")) + 
  scale_shape_manual("", 
                     breaks = c("Raw Data", "15sec Median", "1min Median"),
                     values = c("21" ,"15", "18"))

############################## Figure 3 PCs ##########################

m1 <- ggmap(p0) +
  geom_point(data=pred_block,aes(x=Longitude,y=Latitude,colour=PC1),size=0.2) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=tim.colors(256)) + coord_map() + ggtitle("PC 1")

m2 <- ggmap(p0) +
  geom_point(data=pred_block,aes(x=Longitude,y=Latitude,colour=PC3),size=0.2) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=tim.colors(256), limits=c(-3.5,5.5)) + coord_map() + ggtitle("PC 3")

m3 <- ggmap(p0) +
  geom_point(data=pred_block,aes(x=Longitude,y=Latitude,colour=PC4),size=0.2) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=tim.colors(256), limits=c(-5,5)) + coord_map() + ggtitle("PC 4")

m4 <- ggmap(p0) +
  geom_point(data=pred_block,aes(x=Longitude,y=Latitude,colour=PC5),size=0.2) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=tim.colors(256), limits=c(-6.5,6)) + coord_map() + ggtitle("PC 5")

gridExtra::grid.arrange(m1,m2,m3,m4,ncol=2)

########################## Figure 4 Hourly Regression Means ###########################

df_block = df_block_tseg15sec
keep = c("Y_block_med","locID","Longitude","Latitude","Car","speed","year",
         "month","day","hour","min","sec","wday","t",paste0("PC",1:7),paste0("H",1:4))
df_block[,colnames(df_block)%in%paste0("PC",1:7)] = scale(df_block[,colnames(df_block)%in%paste0("PC",1:7)])
df_block=df_block[order(df_block$t),keep]

load("s_map_preds253.Rda")
load("st_map_preds253.Rda")

hr <- 14

np <- floor(nrow(sp_dat)/nrow(pred_block))

big_dat <- do.call(rbind, rep(list(pred_block), np))
                   
sp_alldat <- data.frame(sp_dat, big_dat)
st_alldat <- data.frame(st_dat, big_dat)
cov_names <- c(paste0("PC",1:7),paste0("H",1:4))

load(file=paste0("linmod_type2.Rda"))
lm_dat <- do.call(rbind,rep(list(pred_block), 4))
lm_dat$hour <- rep(c(9,12,15,18), each=nrow(pred_block))
lm_dat = lm_dat[order(lm_dat$hour),]

lm_dat$H1         <- sin(2*pi*lm_dat$hour/24)
lm_dat$H2         <- cos(2*pi*lm_dat$hour/24)
lm_dat$H3         <- sin(4*pi*lm_dat$hour/24)
lm_dat$H4         <- cos(4*pi*lm_dat$hour/24)

lm_dat$fit <- predict(linmod, newdata = lm_dat)

m1 <- ggmap(p0) +
  geom_point(data=subset(lm_dat, hour==9),aes(x=Longitude,y=Latitude,colour=fit),size=0.2) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=c(0.15,6.8)) + coord_map() + ggtitle("9:00")

m2 <- ggmap(p0) +
  geom_point(data=subset(lm_dat, hour==12),aes(x=Longitude,y=Latitude,colour=fit),size=0.2) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=c(2.1, 3.65)) + coord_map() + ggtitle("12:00")

m3 <- ggmap(p0) +
  geom_point(data=subset(lm_dat, hour==15),aes(x=Longitude,y=Latitude,colour=fit),size=0.2) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=c(2,3.5)) + coord_map() + ggtitle("15:00")

m4 <- ggmap(p0) +
  geom_point(data=subset(lm_dat, hour==18),aes(x=Longitude,y=Latitude,colour=fit),size=0.2) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=c(-0.5, 8.3)) + coord_map() + ggtitle("18:00")

gridExtra::grid.arrange(m1,m2,m3,m4,ncol=2)


################### Figure 6 15min Forecast Comparison (S vs ST model) #####################

kp <- which(floor(df_block$t)==253 & df_block$hour %in% 9:(hr-1))
df_block$resid <- df_block$Y_block_med - predict(linmod, newdata = df_block)

cond_dat <- df_block[kp,c("Longitude", "Latitude", "t", cov_names, "Y_block_med","resid")]

lim1 <- range(quantile(subset(cond_dat, Longitude > -122.22 & Latitude < 37.79)$Y_block_med, c(0.01, 0.99)), quantile(subset(sp_alldat, hour==hr & Longitude > -122.22 & Latitude < 37.79)$Pred_S, c(0.001, 0.999)), quantile(subset(st_alldat, hour==hr & Longitude > -122.22 & Latitude < 37.79)$Pred_St, c(0.001, 0.999)))
lim2 <- range(quantile(subset(cond_dat, Longitude > -122.22 & Latitude < 37.79)$resid), range(subset(sp_alldat, hour==hr & Longitude > -122.22 & Latitude < 37.79)$Pred_S - subset(sp_alldat, hour==hr & Longitude > -122.22 & Latitude < 37.79)$Pred_Xb), range(subset(st_alldat, hour==hr & Longitude > -122.22 & Latitude < 37.79)$Pred_St - subset(st_alldat, hour==hr & Longitude > -122.22 & Latitude < 37.79)$Pred_Xb))

theme_set(theme_bw(base_size=18))
ggmap(p0) +
  geom_point(data=cond_dat,aes(x=Longitude,y=Latitude,colour=Y_block_med),size=1) +
  xlim(c(-122.22, max(pred_block$Longitude))) + ylim(c(37.73, 37.79)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(2,"cm"), legend.key.height=unit(.5,"cm"), legend.text=element_text(size=16), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=lim1) + coord_map() + ggtitle("Observed Log(NO2)")

ggmap(p0) +
  geom_point(data=subset(sp_alldat, hour==hr), aes(x=Longitude,y=Latitude,colour=Pred_S),size=0.5) +
  xlim(c(-122.22, max(pred_block$Longitude))) + ylim(c(37.73, 37.79)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=lim1) + coord_map() + ggtitle("Spatial Forecast")

ggmap(p0) +
  geom_point(data=subset(st_alldat, hour==hr), aes(x=Longitude,y=Latitude,colour=Pred_St),size=0.5) +
  xlim(c(-122.22, max(pred_block$Longitude))) + ylim(c(37.73, 37.79)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(2,"cm"), legend.key.height=unit(.5,"cm"), legend.text=element_text(size=16), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=lim1) + coord_map() + ggtitle("Space-Time Forecast")

ggmap(p0) +
  geom_point(data=cond_dat,aes(x=Longitude,y=Latitude,colour=resid),size=1) +
  xlim(c(-122.22, max(pred_block$Longitude))) + ylim(c(37.73, 37.79)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradient2(name="", limits=lim2) + coord_map() + ggtitle("Observed Log(NO2) - Xbeta")

ggmap(p0) +
  geom_point(data=subset(sp_alldat, hour==hr), aes(x=Longitude,y=Latitude,colour=Pred_S - Pred_Xb),size=0.5) +
  xlim(c(-122.22, max(pred_block$Longitude))) + ylim(c(37.73, 37.79)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradient2(name="", limits=lim2) + coord_map() + ggtitle("Spatial Forecast - Xbeta")

ggmap(p0) +
  geom_point(data=subset(st_alldat, hour==hr), aes(x=Longitude,y=Latitude,colour=Pred_St - Pred_Xb),size=0.5) +
  xlim(c(-122.22, max(pred_block$Longitude))) + ylim(c(37.73, 37.79)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(1.5,"cm"), legend.key.height=unit(.4,"cm"), legend.position="bottom")  + 
  scale_colour_gradient2(name="", limits=lim2) + coord_map() + ggtitle("Space-Time Forecast - Xbeta")

####################### Figure 5 Example 15min Forecasts (ST model) #######################

load("s_map_preds490.Rda")
load("st_map_preds490.Rda")

hr <- 15

np <- floor(nrow(sp_dat)/nrow(pred_block))

big_dat <- do.call(rbind, rep(list(pred_block), np))

sp_alldat <- data.frame(sp_dat, big_dat)
st_alldat <- data.frame(st_dat, big_dat)

kp <- which(floor(df_block$t)==490 & df_block$hour %in% 9:(hr-1))
cond_dat <- df_block[kp,c("Longitude", "Latitude", "t", cov_names, "Y_block_med","resid")]

lim1 <- range(quantile(subset(cond_dat, Longitude > -122.22 & Latitude < 37.79)$Y_block_med, c(0.0001, 0.9999)), quantile(subset(st_alldat, hour==hr & Longitude > -122.22 & Latitude < 37.79)$Pred_St, c(0.001, 0.999)))

theme_set(theme_bw(base_size=18))
ggmap(p0) +
  geom_point(data=cond_dat,aes(x=Longitude,y=Latitude,colour=Y_block_med),size=1) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(2,"cm"), legend.key.height=unit(.5,"cm"), legend.text=element_text(size=16), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=lim1) + coord_map() + ggtitle("Observed Data, 9:00 - 14:45")

ggmap(p0) +
  geom_point(data=subset(st_alldat, hour==hr), aes(x=Longitude,y=Latitude,colour=Pred_St),size=0.3) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(2,"cm"), legend.key.height=unit(.5,"cm"), legend.text=element_text(size=16), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=rev(heat.colors(256)), limits=lim1) + coord_map() + ggtitle("15:00 Forecast")

ggmap(p0) +
  geom_point(data=subset(st_alldat, hour==hr),aes(x=Longitude,y=Latitude,colour=sqrt(Var_St)),size=0.3) +
  xlim(range(pred_block$Longitude)) + ylim(range(pred_block$Latitude)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.key.width=unit(2,"cm"), legend.key.height=unit(.5,"cm"), legend.text=element_text(size=16), legend.position="bottom")  + 
  scale_colour_gradientn(name="",colours=tim.colors(256)) + coord_map() + ggtitle("15:00 Standard Error")


#######################   Figure 8 Example Simulated Paths   ####################

load(file="design_output_all.Rda")
load(file="design_deployment.Rda")

p1 = get_map(location = c(lon =  -122.2323, lat = 37.78477),
             zoom = 12, maptype = 'toner-lite', source="stamen")

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

datfull <- savedatfull[[1]]

ggmap(p1) +
  geom_point(data=subset(datfull,car!='grid'),aes(x=Longitude,y=Latitude,group=car,colour=car,size=car,shape=car)) +
  xlim(range(datfull$Longitude)) + ylim(range(datfull$Latitude)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5), panel.grid = element_blank(),
        legend.position="none")  + 
  scale_color_manual(name="",values=c(ggplotColours(15),"black")) + 
  scale_size_manual(name="",values = c(rep(0.9,15),2)) + 
  scale_shape_manual(name="",values = c(rep(16,15),17)) + coord_map()

#######################   Figure 9 Simulation MSPE   ####################

nstation = 15
ncar = 15
predcar <- colMeans(predcarall)
predst <-  colMeans(predstall)
interpcar <- colMeans(interpcarall)
interpst <- colMeans(interpstall)
predcarCI <- apply(predcarall,2,quantile,prob = c(0.025,0.975))
interpcarCI <- apply(interpcarall,2,quantile,prob = c(0.025,0.975))

plot(NA, xlim = c(1,15), ylim = range(0.32,0.46),ylab = "MSPE",xlab="Number of Monitors")
lines(1:nstation, predst,col="red",lwd=2,lty=1)
lines(1:nstation, interpst,col="red",lwd=2,lty=2)
polygon(c(1:ncar, ncar:1),c(predcarCI[1,], rev(predcarCI[2,])),col=alpha("lightblue",0.5),border =NA,xpd=TRUE)
lines(1:ncar, predcar,col="lightblue4",lwd=2,lty=1)
lines(1:ncar, predcarCI[1,],col="lightblue",lwd=2,lty=1)
lines(1:ncar, predcarCI[2,],col="lightblue",lwd=2,lty=1)
polygon(c(1:ncar, ncar:1),c(interpcarCI[1,], rev(interpcarCI[2,])),col=alpha("skyblue",0.5),border = NA,xpd=TRUE)
lines(1:ncar, interpcar,col="skyblue4",lwd=2,lty=2)
lines(1:ncar, interpcarCI[1,],col="skyblue",lwd=2,lty=1)
lines(1:ncar, interpcarCI[2,],col="skyblue",lwd=2,lty=1)
legend("topright", legend=c("Station Forecast","Station Interpolation","Mobile Forecast","Mobile Interpolation"),
       col = c("red","red","lightblue4","skyblue4"),lwd=2,lty=c(1,2,1,2))




