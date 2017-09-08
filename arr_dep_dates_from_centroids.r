# Calculating "arrival" and "departure" dates on breeding and wintering grounds
# from eBird centroid data (from La Sorte et al. 2016)


# Based on the figure created below, this current method does not work.
# May instead need to use % of deviation from max or min latitude instead of
# differences/derivatives.


filelist = list.files('data/ebird')

pdf('figs/centroid_latitudes.pdf', height = 8, width = 10)
par(mfrow = c(3, 4), mar = c(2.5, 2.5, 2.5, 1), oma = c(4, 4, 0, 0))
panel = 0

mig_summary = data.frame(species = NULL, w_lat = NULL, b_lat = NULL, w_depdate = NULL,
                         b_arrdate = NULL, b_depdate = NULL, w_arrdate = NULL)
for (f in filelist) {
  tmp = read.csv(paste('data/ebird/', f, sep = ''), header = T)
  tmp$species = substr(f, 1, nchar(f)-4)

  tmp$diff = c(NA, diff(tmp$lat))
  tmp$diff2 = c(NA, diff(tmp$diff))
  
  #Spring migration
  s_maxd2 = max(tmp$diff2[tmp$day < 180], na.rm = T)
  s_mind2 = min(tmp$diff2[tmp$day < 180], na.rm = T)
  
  #Wintering grounds
  f_maxd2 = max(tmp$diff2[tmp$day > 180], na.rm = T)
  f_mind2 = min(tmp$diff2[tmp$day > 180], na.rm = T)
  
  df = data.frame(species = tmp$species[1],
                  w_lat = min(tmp$lat),
                  b_lat = max(tmp$lat),
                  w_depdate = tmp$day[tmp$diff2 == s_maxd2 & !is.na(tmp$diff2) & tmp$day < 180],
                  b_arrdate = tmp$day[tmp$diff2 == s_mind2 & !is.na(tmp$diff2) & tmp$day < 180],
                  b_depdate = tmp$day[tmp$diff2 == f_mind2 & !is.na(tmp$diff2) & tmp$day > 180],
                  w_arrdate = tmp$day[tmp$diff2 == f_maxd2 & !is.na(tmp$diff2) & tmp$day > 180])
  
  mig_summary = rbind(mig_summary, df)
  
  plot(tmp$day, tmp$lat, type = 'l', xlab = "", ylab = "")
  abline(v = c(df$w_depdate, df$b_arrdate, df$b_depdate, df$w_arrdate), lty = 'dotted', col = 'red')
  mtext(tmp$species[1], 3)
  panel = panel+1
  if (panel%%12 == 0) {
    mtext("Julian day", 1, outer = TRUE, line = 1, cex = 2)
    mtext("Latitude", 2, outer = TRUE, line = 1, cex = 2)
  }
  
}

dev.off()