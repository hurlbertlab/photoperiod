# Preliminary script for reading in geolocator data and extracting daylengths
# based on date and latitude

# Load libraries
library(lubridate)
library(geosphere)
library(dplyr)

# Functions

# preFormat takes a raw Movebank dataset and adds date, year, julian day, 
# and daylength columns, and eliminates unneeded fields.
preFormat = function(data) {
  dataOut = select(data, individual.taxon.canonical.name, individual.local.identifier,
                   location.long, location.lat, timestamp)
  dataOut$date = as.POSIXct(strptime(dataOut$timestamp, format = "%m/%d/%y %H:%M"))
  dataOut$jd = yday(dataOut$date)
  dataOut$daylength = daylength(dataOut$location.lat, dataOut$jd)
  dataOut$year = format(dataOut$date, '%Y')
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
                         color,             # color of experienced daylength line
                         refLines = FALSE,  # include daylength at reference latitudes
                         refLatitudes = NULL,   # vector of min, max ref latitudes
                         ref12hr = FALSE,     # ref line for 12 hr photoperiod
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
  if(new) {
    plot(temp$jd, temp$daylength, type = 'l', xlim = c(1, 365), ylim = c(9, maxDaylength + 2), 
         xlab = "Julian day", ylab = "Day length (h)", col = color, las = 1, ...)
  } else {
    points(temp$jd, temp$daylength, type = 'l', ...)
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
