# Install sparklyr
# install.packages("devtools")
# devtools::install_github("rstudio/sparklyr")
# library(sparklyr)
# spark_install()

################
# 1 - Set-up   #
################
# Load required libraries
library(sparklyr)
library(dplyr)
library(flows)
library(ggmap)
library(rgeos)
library(circlize)
library(RColorBrewer)
library(reshape2)

# Connect and load the dataset to local Spark
sc <- spark_connect(master = "local")
setwd("~/R3PI/")
dataset <- spark_read_parquet(sc, "file1", path = "zubie_trips_anonymous/")

# Write the dataset back to R and close local Spark
car_data <- dataset %>% collect
spark_disconnect(sc)


#####################
# 2 - Data clean-up #
#####################

# Count/display unique devices   
devices <- distinct(car_data,device_key)
summary(devices)
View(devices)

# Remove duplicate rows
car_data <- distinct(car_data[,-50])

# Remove observations with NAs in duration (they also have no end point) 
car_data <- filter(car_data,!is.na(duration_seconds))

# Remove observations with start/end point latitude 0
car_data <- filter(car_data,end_point_latitude!=0 & start_point_latitude!=0)

# Remove observations with top speed = 0
car_data <- filter(car_data,top_speed!=0)

# Remove observations with gps_distance = 0
car_data <- filter(car_data,!is.na(gps_distance))

# Convert times from strings to timestamps
car_data$end_point_timestamp <- as.POSIXct(strptime(car_data$end_point_timestamp, "%Y-%m-%d %H:%M:%S", tz = "Europe/London")) 
car_data$start_point_timestamp <- as.POSIXct(strptime(car_data$start_point_timestamp, "%Y-%m-%d %H:%M:%S", tz = "Europe/London")) 
range(car_data$start_point_timestamp) # check period covered in the dataset

# Add a field with a weekend flag
days_set <- data.frame('day'=weekdays(car_data$start_point_timestamp))
lookup_days <- data.frame('day'=c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'),
                          'flag'=c(0,0,0,0,0,1,1))
car_data$weekend_flag <- left_join(days_set,lookup_days,by='day')[,2]

# Add field day number
car_data$day <- as.POSIXlt(car_data$start_point_timestamp)$yday

# Fill in the missing START point address cities using Google Maps API 
# ATTENTION: it may take 7-8 min
missing <- as.data.frame(car_data[is.na(car_data$start_point_address_city),c(44,43)])
res <- lapply(with(missing,paste(start_point_latitude,start_point_longitude, sep = ",")), 
                geocode, output = "more")
cities <- transform(missing, city = as.matrix(lapply(sapply(res,"[[","postal_town"), function(x){as.character(x)[1]}),ncol=1))
car_data[is.na(car_data$start_point_address_city),]$start_point_address_city <- cities$city
# we delete the few observations where Google Maps failed to return a town name
car_data <- filter(car_data,!is.na(start_point_address_city))

# Fill in the missing END point address cities using Google Maps API
# ATTENTION: it may take 7-8 min
missing2 <- as.data.frame(car_data[is.na(car_data$end_point_address_city),c(10,9)])
res2 <- lapply(with(missing2,paste(end_point_latitude,end_point_longitude, sep = ",")), 
              geocode, output = "more")
cities2 <- transform(missing2, city = as.matrix(lapply(sapply(res2,"[[","postal_town"), function(x){as.character(x)[1]}),ncol=1),
                    zip = as.matrix(lapply(sapply(res2,"[[","postal_code"), function(x){as.character(x)[1]}),ncol=1))
car_data[is.na(car_data$end_point_address_city),]$end_point_address_city <- cities2$city
# we delete the few observations where Google Maps failed to return a town name
car_data <- filter(car_data,!is.na(end_point_address_city))

# Fill in the missing zipcodes from a lookup list 
# Load lookup list (obtained from Free Map Tools) 
# Attention: unzip and place in the same directory the file ukpostcodes.csv
postcodes <- read.csv("ukpostcodes.csv",stringsAsFactors = F)
names(postcodes) <- c('id','zip','lat','lon')
postcodes <- filter(postcodes,!is.na(lat) & !is.na(lon))
sp_lookup <- SpatialPoints(postcodes[,c(4,3)])
# Fill in gaps START zip code
missing <- as.data.frame(car_data[is.na(car_data$start_point_address_zipcode),c(44,43)])
names(missing) <- c('lon','lat')
sp_missing <- SpatialPoints(missing)
car_data[is.na(car_data$start_point_address_zipcode),]$start_point_address_zipcode <- 
  postcodes[apply(gDistance(sp_missing, sp_lookup, byid=TRUE),2,which.min),2]
# Fill in gaps END zip code
missing2 <- as.data.frame(car_data[is.na(car_data$end_point_address_zipcode),c(10,9)])
names(missing2) <- c('lon','lat')
sp_missing2 <- SpatialPoints(missing2)
car_data[is.na(car_data$end_point_address_zipcode),]$end_point_address_zipcode <- 
  postcodes[apply(gDistance(sp_missing2, sp_lookup, byid=TRUE),2,which.min),2]


#################
# 3 - Analysis  #
#################

# Select fields that will be used for the analysis
analysis <- as.data.frame(select(car_data,duration_seconds,end_point_address_city,end_point_address_state,
                   end_point_address_zipcode,end_point_latitude,end_point_longitude,
                   end_point_timestamp,gps_distance,vehicle_nickname,weekend_flag,
                   start_point_address_city,start_point_address_state,start_point_address_zipcode,
                   start_point_latitude,start_point_longitude,start_point_timestamp,day))

# Separate weekdays and weekends
week <- filter(analysis,weekend_flag==0)[,-10]
weekend <- filter(analysis,weekend_flag==1)[,-10]

# Create full origin-destination matrix
od_week <- week %>% group_by(start_point_address_city,end_point_address_city) %>% 
  summarize(count=n()) %>% arrange(desc(count))
od_weekend <- weekend %>% group_by(start_point_address_city,end_point_address_city) %>% 
  summarize(count=n()) %>% arrange(desc(count))
# Transform it into an adjecency matrix
adj_week <- prepflows(mat = od_week, i = "start_point_address_city", 
                      j = "end_point_address_city", fij = "count")
adj_weekend <- prepflows(mat = od_weekend, i = "start_point_address_city", 
                      j = "end_point_address_city", fij = "count")

# Create reduced origin-destination matrix with the M towns with the biggest flows
# weekdays
M <- 29
ordered <- sort(rowSums(adj_week)+colSums(adj_week),decreasing = T)
lookup <- data.frame('orig'= names(ordered[1:M]), 'new'=names(ordered[1:M]),stringsAsFactors = F)
lookup <- rbind(lookup, data.frame('orig'= names(ordered[M+1:length(ordered)]), 'new'='Other', stringsAsFactors = F))
od_week <- week %>% left_join(lookup,by=c("start_point_address_city"="orig")) %>%
  left_join(lookup,by=c("end_point_address_city"="orig"),suffix=c('s','e')) %>%
  group_by(news,newe) %>% summarize(count=n()) %>% arrange(desc(count)) %>% as.data.frame(stringsAsFactors = F)
adjR_week <- prepflows(mat = as.data.frame(od_week, stringsAsFactors = F),i = "news",j = "newe", fij = "count")
# Output circular chart weekdays
# format data
city_names <- unlist(distinct(od_week,news),use.names = F)
city_short_names <- city_names
city_short_names[13] <- paste(strwrap(city_names[13], width=12), collapse="\n") # Shorten long name
city_clk <- c("",city_short_names[2:(M+1)])
city_obl <- c(city_short_names[1],rep("",M))
# Create chord chart
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(-0.2, 0.2), points.overflow.warning = FALSE)
par(mar = rep(1, 4), oma=rep(0,4))
chordDiagram(x = adjR_week, grid.col = c(brewer.pal(10,"PRGn"),brewer.pal(10,"RdYlGn"),brewer.pal(10,"BrBG")),
             transparency = 0.45, order = city_names, directional = 1, 
             direction.type = c("arrows", "diffHeight"), diffHeight  = -0.04,
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.1),
             link.arr.type = "big.arrow", link.sort = TRUE, link.largest.ontop = TRUE)
# Add labels and axes
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    reg1 = city_clk[city_names == sector.index]
    reg2 = city_obl[city_names == sector.index]
    circos.text(x = mean(xlim), y = 5, adj = c(1,0.5),
                labels = reg1, facing = "reverse.clockwise", cex = ifelse(nchar(reg1) > 10,1.2,1.4))
    circos.text(x = mean(xlim), y = 10, adj = c(0,0.5), 
                labels = reg2, facing = "bending.inside", cex = 2)
    circos.axis(h = "top", 
                major.at = seq(from = 0, to = xlim[2], by = 100),
                labels.cex = 0.6,
                minor.ticks = 1, major.tick.percentage = 0.5,
                labels.niceFacing = T, labels.facing = "clockwise")
  }
)
text(x =0.6, y = 0.8, pos = 4, cex = 4, col = 'darkblue', labels = "Week")
dev.copy2pdf(file = "Circular_chart_week2.pdf", height=15, width=15)

# Create reduced origin-destination matrix with the M towns with the biggest flows
# weekends
M <- 20
ordered <- sort(rowSums(adj_weekend)+colSums(adj_weekend),decreasing = T)
lookup <- data.frame('orig'= names(ordered[1:M]), 'new'=names(ordered[1:M]),stringsAsFactors = F)
lookup <- rbind(lookup, data.frame('orig'= names(ordered[M+1:length(ordered)]), 'new'='Other', stringsAsFactors = F))
od_weekend <- weekend %>% left_join(lookup,by=c("start_point_address_city"="orig")) %>%
  left_join(lookup,by=c("end_point_address_city"="orig"),suffix=c('s','e')) %>%
  group_by(news,newe) %>% summarize(count=n()) %>% arrange(desc(count)) %>% as.data.frame(stringsAsFactors = F)
adjR_weekend <- prepflows(mat = as.data.frame(od_weekend, stringsAsFactors = F),i = "news",j = "newe", fij = "count")
# Output chart weekends
# Format data
city_names <- unlist(distinct(od_weekend,news),use.names = F)
city_short_names <- city_names
city_clk <- c("",city_short_names[2:(M+1)])
city_obl <- c(city_short_names[1],rep("",M))
# Create chord chart
circos.clear()
circos.par(start.degree = 90, gap.degree = 2, track.margin = c(-0.3, 0.3), points.overflow.warning = FALSE)
par(mar = rep(1, 4), oma=rep(0,4))
chordDiagram(x = adjR_weekend, grid.col = c(brewer.pal(10,"PRGn"),brewer.pal(10,"RdYlGn"),brewer.pal(10,"BrBG"))[1:(M+1)],
             transparency = 0.45, order = city_names, directional = 1, 
             direction.type = c("arrows", "diffHeight"), diffHeight  = -0.04,
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.1),
             link.arr.type = "big.arrow", link.sort = TRUE, link.largest.ontop = TRUE)
# Add labels and axes
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    reg1 = city_clk[city_names == sector.index]
    reg2 = city_obl[city_names == sector.index]
    circos.text(x = mean(xlim), y = 5, adj = c(1,0.5),labels.niceFacing = T,
                labels = reg1, facing = "reverse.clockwise", cex = 1.4)
    circos.text(x = mean(xlim), y = 10, adj = c(0,0.5), 
                labels = reg2, facing = "bending.inside", cex = 2)
    circos.axis(h = "top", 
                major.at = seq(from = 0, to = xlim[2], by = 100),
                labels.cex = 0.8,labels.away.percentage = 0.1,
                minor.ticks = 1,
                labels.niceFacing = T, labels.facing = "clockwise")
  }
)
text(x =0.6, y = 0.8, pos = 4, cex = 4, col = 'darkblue', labels = "Weekend")
dev.copy2pdf(file = "Circular_chart_weekend2.pdf", height=15, width=15)


# Hourly profile main London-London connections
# weekdays
London_week <- week %>% filter(start_point_address_city=='London' & end_point_address_city=='London')
cars <- distinct(London_week,vehicle_nickname)[,1] 
n <- length(cars)
hcount <- data.frame(matrix(ncol = n, nrow = 24),stringsAsFactors = F)
names(hcount) <- cars
# count per car
for(i in 0:23){
  for(j in 1:n){
    tmp <- filter(London_week, 
                   as.POSIXlt(start_point_timestamp)$hour < (i+1) & 
                     as.POSIXlt(end_point_timestamp)$hour >= i &
                     vehicle_nickname == cars[j])
    hcount[i+1,j] <- dim(tmp)[1]
  }
}
# plot
# format & flatten data
hcount['Hour'] <- row.names(hcount)
bardata <- melt(hcount, id="Hour") 
bardata$Hour <- factor(bardata$Hour, levels = as.character(1:24))
bardata$Hour <- as.integer(bardata$Hour)
ax <- c("",0); for(i in 1:24){ax <- append(ax,c("",i))} #create manual axis
# output barchart
ggplot(data=bardata, aes(x=Hour, y=value, fill=variable)) +
  geom_bar(stat="identity")+
  theme_classic()+
  ylab('Trips in the period (weekdays)')+
  scale_fill_manual(values=brewer.pal(n,"Paired"),name="")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,13),breaks=seq(0,20,1))+
  scale_x_continuous(breaks=c(seq(0,24.5,0.5)), labels=ax, limits=c(7,24))+
  theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))+
  theme(legend.text=element_text(size=14),legend.position=c(0.2,0.7))+
  theme(axis.title.y = element_text(size = rel(1.8)))+
  theme(axis.title.x = element_text(size = rel(1.8)))+
  theme(axis.text.x = element_text(size = rel(1.8)))+      
  theme(axis.text.y = element_text(size = rel(1.8)))+
  guides(fill=guide_legend(keywidth=1,keyheight=1,default.unit="cm"))+
  theme(axis.ticks.x = element_line(color = rep(c('black', NA),24)))
ggsave("Profile_week.pdf",width=898/72,height=610/72)


# Hourly profile main connection
# weekends
London_weekend <- weekend %>% filter(start_point_address_city=='London' & end_point_address_city=='London')
cars <- distinct(London_weekend,vehicle_nickname)[,1] 
n <- length(cars)
hcount <- data.frame(matrix(ncol = n, nrow = 24),stringsAsFactors = F)
names(hcount) <- cars
# count per car
for(i in 0:23){
  for(j in 1:n){
    tmp <- filter(London_weekend, 
                  as.POSIXlt(start_point_timestamp)$hour < (i+1) & 
                    as.POSIXlt(end_point_timestamp)$hour >= i &
                    vehicle_nickname == cars[j])
    hcount[i+1,j] <- dim(tmp)[1]
  }
}
# plot
# format & flatten the data
hcount['Hour'] <- row.names(hcount)
bardata <- melt(hcount, id="Hour") 
bardata$Hour <- factor(bardata$Hour, levels = as.character(1:24))
bardata$Hour <- as.integer(bardata$Hour)
ax <- c("",0); for(i in 1:24){ax <- append(ax,c("",i))} #create manual axis
# output barchart
ggplot(data=bardata, aes(x=Hour, y=value, fill=variable)) +
  geom_bar(stat="identity")+
  theme_classic()+
  ylab('Trips in the period (weekend)')+
  scale_fill_manual(values=brewer.pal(n,"Paired"),name="")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,13),breaks=seq(0,20,1))+
  scale_x_continuous(breaks=c(seq(0,24.5,0.5)), labels=ax, limits=c(7,24))+
  theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))+
  theme(legend.text=element_text(size=14),legend.position=c(0.7,0.7))+
  theme(axis.title.y = element_text(size = rel(1.8)))+
  theme(axis.title.x = element_text(size = rel(1.8)))+
  theme(axis.text.x = element_text(size = rel(1.8)))+      
  theme(axis.text.y = element_text(size = rel(1.8)))+
  guides(fill=guide_legend(keywidth=1,keyheight=1,default.unit="cm"))+
  theme(axis.ticks.x = element_line(color = rep(c('black', NA),24)))
ggsave("Profile_weekend.pdf",width=898/72,height=610/72)

# Check car sharing possibilities
# weekdays
T <- dim(week)[1]
delta_t <- 60 # acceptable time difference (in min) for car sharing
delta_k <- 2000 # acceptable distance (in m) for car sharing
# Specify the projections to calculate planar distances in the UK
epsg.27700 <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
wgs.84    <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# We loop to record the matches for a range of values of delta_k
trend <- data.frame(matrix(ncol = 2, nrow = 10),stringsAsFactors = F) # data frame to record the results
for(delta_k in seq(from=2000, to=20000, by=2000)){
# Initialize matches dataframe
matches <- data.frame(matrix(ncol = 11, nrow = 0),stringsAsFactors = F)
names(matches) <- c('car1','car2','start_city_car1','start_city_car2','end_city',
                    'start_time_car1','start_time_car2','start_lon_car1','start_lat_car1',
                    'start_lon_car2','start_lat_car2')
for(i in 1:(T-1)){
  current <- week[i,]
  cmp <- week[i+1:T,]
  # Pre-filter by vehicle name, start points and day 
  cmp <- cmp %>% filter(vehicle_nickname != current$vehicle_nickname & 
                        end_point_address_city == current$end_point_address_city &
                        abs(day-current$day)< 2)
  if(dim(cmp)[1]!= 0){
    # within the same day, filter by time
    cmp <- cmp %>% 
      filter(abs(as.numeric(cmp$start_point_timestamp-
                          current$start_point_timestamp,units="mins")) < delta_t)
    if(dim(cmp)[1]!= 0){
      # Read current lon/lat coordinates, project into the UK planar CRS & filter according to distance 
      cmp <- cmp %>% 
        filter(gDistance(spTransform(SpatialPoints(data.frame(lon=current$start_point_longitude,
                                                              lat=current$start_point_latitude,
                                                              stringsAsFactors = F),
                                                   proj4string = CRS(wgs.84)),CRS(epsg.27700)),
                         spTransform(SpatialPoints(data.frame(lon=cmp$start_point_longitude,
                                                              lat=cmp$start_point_latitude,
                                                              stringsAsFactors = F),
                                                   proj4string = CRS(wgs.84)),CRS(epsg.27700))) <delta_k)
      if(dim(cmp)[1]!= 0){
        # if there are some matches
        for(j in dim(cmp)[1]){
          # record each match
          matches <- rbind(matches, data.frame(car1=current$vehicle_nickname, 
                                               car2=cmp[j,]$vehicle_nickname,
                                               start_city_car1=current$start_point_address_city,
                                               start_city_car2=cmp[j,]$start_point_address_city,
                                               end_city=current$start_point_address_city,
                                               start_time_car1=current$start_point_timestamp,
                                               start_time_car2=cmp[j,]$start_point_timestamp,
                                               start_lon_car1=current$start_point_longitude,
                                               start_lat_car1=current$start_point_latitude,
                                               start_lon_car2=cmp[j,]$start_point_longitude,
                                               start_lat_car2=cmp[j,]$start_point_latitude,
                                               stringsAsFactors = F),stringsAsFactors = F)
        }
      }
    }
  }
}
# record number of matches
trend[delta_k/2000,1] <- dim(matches)[1]
# save csv file with the record of the matches
write.csv(matches,paste0("matches_",delta_k/1000,"km_weekdays.csv"))
}

# Same analysis but for weekends
# weekend
T <- dim(weekend)[1]
delta_t <- 60 # acceptable time difference (in min) for car sharing
delta_k <- 2000 # acceptable distance (in m) for car sharing
# Specify the projections to calculate planar distances in the UK
epsg.27700 <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
wgs.84    <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# We loop to record the matches for a range of values of delta_k
for(delta_k in seq(from=2000, to=20000, by=2000)){
  # Initialize matches dataframe
  matches <- data.frame(matrix(ncol = 11, nrow = 0),stringsAsFactors = F)
  names(matches) <- c('car1','car2','start_city_car1','start_city_car2','end_city',
                      'start_time_car1','start_time_car2','start_lon_car1','start_lat_car1',
                      'start_lon_car2','start_lat_car2')
  for(i in 1:(T-1)){
    current <- weekend[i,]
    cmp <- weekend[i+1:T,]
    # Pre filter by vehicle name, start & end points and day 
    cmp <- cmp %>% filter(vehicle_nickname != current$vehicle_nickname & 
                            end_point_address_city == current$end_point_address_city &
                            abs(day-current$day)< 2)
    if(dim(cmp)[1]!= 0){
      # within the same day, filter by time
      cmp <- cmp %>% 
        filter(abs(as.numeric(cmp$start_point_timestamp-
                                current$start_point_timestamp,units="mins")) < delta_t)
      if(dim(cmp)[1]!= 0){
        # Read current lon/lat coordinates, project into the UK planar CRS & filter according to distance 
        cmp <- cmp %>% 
          filter(gDistance(spTransform(SpatialPoints(data.frame(lon=current$start_point_longitude,
                                                                lat=current$start_point_latitude,
                                                                stringsAsFactors = F),
                                                     proj4string = CRS(wgs.84)),CRS(epsg.27700)),
                           spTransform(SpatialPoints(data.frame(lon=cmp$start_point_longitude,
                                                                lat=cmp$start_point_latitude,
                                                                stringsAsFactors = F),
                                                     proj4string = CRS(wgs.84)),CRS(epsg.27700))) <delta_k)
        if(dim(cmp)[1]!= 0){
          # if there are some matches
          for(j in dim(cmp)[1]){
            # record each match
            matches <- rbind(matches, data.frame(car1=current$vehicle_nickname, 
                                                 car2=cmp[j,]$vehicle_nickname,
                                                 start_city_car1=current$start_point_address_city,
                                                 start_city_car2=cmp[j,]$start_point_address_city,
                                                 end_city=current$start_point_address_city,
                                                 start_time_car1=current$start_point_timestamp,
                                                 start_time_car2=cmp[j,]$start_point_timestamp,
                                                 start_lon_car1=current$start_point_longitude,
                                                 start_lat_car1=current$start_point_latitude,
                                                 start_lon_car2=cmp[j,]$start_point_longitude,
                                                 start_lat_car2=cmp[j,]$start_point_latitude,
                                                 stringsAsFactors = F),stringsAsFactors = F)
          }
        }
      }
    }
  }
  # record result
  trend[delta_k/2000,2] <- dim(matches)[1]
  # save csv file with the matches
  write.csv(matches,paste0("matches_",delta_k/1000,"km_weekend.csv"))
}

# Print results
# format & flatten dataframe
names(trend) <- c('weekdays','weekend')
trend$distance <- seq(2,20,by=2)
trend <- melt(trend, id="distance") 
# plot line chart
ggplot(data=trend, aes(x=distance, y=value, group=variable)) +
  geom_line(aes(colour = variable),size=1.5) +
  theme_classic()+ 
  theme(axis.title.y=element_text(margin=margin(0,20,0,0)),axis.title.x=element_text(margin=margin(20,0,0,0)))+
  ylab("Number of potential car sharings")+
  xlab("Maximum start distance allowed between car sharers (km)")+
  scale_x_continuous(breaks=seq(0,20,2))+
  scale_y_continuous(breaks=seq(0,35,5),
                     labels=function(x) format(x, big.mark = " ", scientific = FALSE))+
  scale_color_manual(values=brewer.pal(3,'BuGn')[c(2,3)],name="")+
  theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))+
  theme(panel.grid.major = element_line(colour="grey80", linetype="dashed",size=0.2),
        panel.grid.major.x = element_line(colour="white", linetype="dashed",size=0))+
  theme(legend.text=element_text(size=rel(2)),legend.direction="vertical",
        legend.position=c(0.3,0.85),
        legend.key.height=unit(3,"line"),legend.key.size=unit(1,"cm"),
        legend.background = element_rect(fill="transparent"))+
  theme(axis.title.y = element_text(size = rel(2)))+
  theme(axis.title.x = element_text(size = rel(2)))+
  theme(axis.text.x = element_text(size = rel(2)))+      
  theme(axis.text.y = element_text(size = rel(2)))      
# save file
ggsave("Matches.pdf",width=898/72,height=610/72)
