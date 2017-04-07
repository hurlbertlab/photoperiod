# Preliminary script for reading in geolocator data and extracting daylengths
# based on date and latitude

# Load libraries
library(lubridate)
library(geosphere)
library(dplyr)

# Functions

# preFormat takes a raw Movebank dataset and adds date, year, julian day, 
# and daylength columns, and eliminates unneeded fields.
preFormat = function(data, format = "%m/%d/%y %H:%M") {
  dataOut = select(data, individual.taxon.canonical.name, individual.local.identifier,
                   location.long, location.lat, timestamp)
  dataOut$date = as.POSIXct(strptime(dataOut$timestamp, format = format))
  dataOut$jd = yday(dataOut$date)
  dataOut$daylength = daylength(dataOut$location.lat, dataOut$jd)
  dataOut$year = as.numeric(format(dataOut$date, '%Y'))
  names(dataOut)[1:4] = c('species', 'individual', 'lon', 'lat')
  return(dataOut[, c('species', 'individual', 'lon', 'lat',
                     'date', 'jd', 'year', 'daylength')])
}

# plotDayLength takes a Movebank dataset that has been formatted with preFormat,
# and plots daylength experienced by an individual as a function of julian day
plotDaylength = function(data,              # preFormatted movement data
                         individualID = NULL,      # unique individual to plot
                         years = NULL,             # year or years of data to plot
                         new = TRUE,        # create new plot or add to existing
                         color = 'purple',             # color of experienced daylength line
                         refLines = FALSE,  # include daylength at reference latitudes
                         refLatitudes = NULL,   # vector of min, max ref latitudes
                         ref12hr = FALSE,     # ref line for 12 hr photoperiod
                         movingAvg = 1,     # number of days over which to calc a moving average
                         ...) {
  if (!is.null(individualID) & !is.null(years)) {
    temp = filter(data, individual == individualID, year %in% years)
  } else if (!is.null(individualID) & is.null(years)) {
    temp = filter(data, individual == individualID)
  } else if (is.null(individualID) & !is.null(years)) {
    temp = filter(data, year %in% years)
  } else {
    temp = data
  }
  temp = temp[order(temp$jd, decreasing = F),]
  maxDaylength = max(temp$daylength, na.rm = T)

  temp$daylength = stats::filter(temp$daylength, rep(1/movingAvg, movingAvg), sides = 2)

    if(new) {
    plot(temp$jd, temp$daylength, type = 'l', xlim = c(1, 365), ylim = c(9, maxDaylength + 2), 
         xlab = "Julian day", ylab = "Day length (h)", col = color, las = 1, ...)
  } else {
    points(temp$jd, temp$daylength, type = 'l', col = color, ...)
  }
  if (refLines) {
    if (is.null(refLatitudes)) {
      maxLatitude = max(temp$lat, na.rm = T)
      minLatitude = min(temp$lat, na.rm = T)
    } else {
      maxLatitude = refLatitudes[2]
      minLatitude = refLatitudes[1]
    }
    
    points(1:365, daylength(maxLatitude, 1:365), type = 'l', col = 'gray50', lty='dashed', lwd = 2)
    points(1:365, daylength(minLatitude, 1:365), type = 'l', col = 'gray80', lty='dashed', lwd = 2)
    legend("topleft", legend = c(round(maxLatitude, 1), round(minLatitude, 1)), lwd = c(2, 2),
           col = c('gray50', 'gray80'), lty = c('dashed', 'dashed'))
  }
  if (ref12hr) {
    abline(h = 12, col = 'red')
  }
}


# ----------------------------------------------------------------------------------

# GPS Data processing

# Read in blackpoll data
bp = read.csv('data/file-178761964.csv', stringsAsFactors = F, header = T)
bp2 = preFormat(bp)

os = read.csv('data/file-108623259.csv', stringsAsFactors = F, header = T)
os2 = preFormat(os)


# Blackpoll
plotDaylength(bp2, "A", years = c(2013, 2014), new = T, color = 'purple', lwd = 3, 
              refLines = TRUE, refLatitudes = c(11.5, 43.9))

# Osprey
os_summary = os2 %>% group_by(individual, year) %>%
  summarize(latRange = max(lat, na.rm = T) - min(lat, na.rm = T), 
            nRecs = n(), firstJD = min())

record_threshold = 2000
osp_inds = unique(os2$individual)

pdf('figs/osprey_individuals.pdf', height = 8, width = 10)
par(mfrow = c(3,4), mar = c(2.5, 2.5, 2.5, 1), oma = c(4, 4, 0, 0)) 
panel = 0
for (o in osp_inds) {
  tmp = filter(os2, individual == o)
  recs_by_year = tmp %>% group_by(year) %>% 
    summarize(latRange = max(lat, na.rm = T) - min(lat, na.rm = T), 
              nRecs = n())
  recs_by_year$latRange[is.na(recs_by_year$latRange)] = 0
  year = recs_by_year$year[recs_by_year$latRange == max(recs_by_year$latRange) &
                             recs_by_year$nRecs >= record_threshold]
  if (length(year) > 0) {
    plotDaylength(os2, o, years = year, new = T, col = 'purple', lwd = 4, 
                  refLines = TRUE, ref12hr = TRUE)
    mtext(paste(o, year, sep = ', '), 3)
    panel = panel+1
    if (panel%%12 == 0) {
      mtext("Julian day", 1, outer = TRUE, line = 1, cex = 2)
      mtext("Day length (h)", 2, outer = TRUE, line = 1, cex = 2)
    }
  }
}
dev.off()


# Swainson's hawk
sw = read.csv("data/Swainson's Hawks.csv")
swha = preFormat(dat, format = "%Y-%m-%d %H:%M:%S")

rec_threshold = 50

swct = count(swha, individual, year) %>%
        group_by(individual) %>% 
        summarize(ct = sum(n >= rec_threshold)) %>% 
        filter(ct==2) %>% 
        data.frame()

sw_inds = swct$individual

pdf('figs/swainsonshawk_individuals.pdf', height = 8, width = 10)
par(mfrow = c(3,4), mar = c(2.5, 2.5, 2.5, 1), oma = c(4, 4, 0, 0)) 
panel = 0
for (s in sw_inds) {
  tmp = filter(swha, individual == s)
  plotDaylength(swha, s, years = unique(tmp$year), new = T, col = 'purple', lwd = 4, 
                refLines = TRUE, ref12hr = TRUE)
  mtext(paste(s, unique(tmp$year)[1], sep = ', '), 3)
  panel = panel+1
  if (panel%%12 == 0) {
    mtext("Julian day", 1, outer = TRUE, line = 1, cex = 2)
    mtext("Day length (h)", 2, outer = TRUE, line = 1, cex = 2)
  }
}
dev.off()




# ---------------------------------------------------------------------------------

# eBird data processing

# Read in data
filelist = list.files('data/ebird')

ebird = c()
for (f in filelist) {
  tmp = read.csv(paste('data/ebird/', f, sep = ''), header = T)
  tmp$species = substr(f, 1, nchar(f)-4)
  ebird = rbind(ebird, tmp)
}
names(ebird)[1] = 'jd'
ebird$daylength = daylength(ebird$lat, ebird$jd)

splist = unique(ebird$species)

pdf('figs/centroid_photoperiods.pdf', height = 8, width = 10)
par(mfrow = c(3, 4), mar = c(2.5, 2.5, 2.5, 1), oma = c(4, 4, 0, 0))
panel = 0
for (s in splist) {
  centroids = filter(ebird, species == s)
  plotDaylength(centroids, new = T, col = 'purple', lwd = 4, refLines = TRUE, ref12hr = TRUE)
  mtext(s, 3)
  panel = panel+1
  if (panel%%12 == 0) {
    mtext("Julian day", 1, outer = TRUE, line = 1, cex = 2)
    mtext("Day length (h)", 2, outer = TRUE, line = 1, cex = 2)
  }
}
dev.off()

#--------------------------------------------------------------------------------
# Comparing Blackpoll warbler eBird centroid to geolocated individuals

blpw.cent = read.csv('data/ebird/Blackpoll_Warbler.csv', header = T)
names(blpw.cent)[1] = 'jd'
blpw.cent$daylength = daylength(blpw.cent$lat, blpw.cent$jd)


pdf('figs/ebird_v_geoloc_blackpoll.pdf', height = 5, width = 6)
par(mar = c(4,4,1,1), oma = c(0,0,0,0), mgp = c(2.5, 1, 0))
# Centroid
plotDaylength(blpw.cent, new = T, col = 'purple', lwd = 6, refLines = TRUE, ref12hr = TRUE)

# Geolocated individuals
plotDaylength(bp2, "A", years = c(2013, 2014), new = F, col = 'blue', lwd = 2, 
              refLines = F, movingAvg = 10)
plotDaylength(bp2, "B", years = c(2013, 2014), new = F, col = 'skyblue', lwd = 2, 
              refLines = F, movingAvg = 10)
plotDaylength(bp2, "E", years = c(2013, 2014), new = F, col = 'darkblue', lwd = 2, 
              refLines = F, movingAvg = 20)
# Reference lines for geolocated individuals
points(1:365, daylength(12, 1:365), type = 'l', lty = 'dotted', lwd = 2, col = 'black')
points(1:365, daylength(45, 1:365), type = 'l', lty = 'dotted', lwd = 2, col = 'gray50')
text(330, 21.5, "Geolocation")
text(30, 19.5, "Centroid")
legend(300, 21.3, legend = c(45, 12), lty = 'dotted', lwd = 2, col = c('gray50', 'black'), bty = 'n')
dev.off()