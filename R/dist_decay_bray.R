##------------------------------------------------------------
## recreates distance-decay analyses for
## american gut project
## - ashkaan fahimipour, may 17 2017
##------------------------------------------------------------
## set working directory
setwd('/Users/Ashkaan/Dropbox/AGP/ag_distance_decay/') ## CHANGE ME to your wd
set.seed(777)

## load libraries for plotting
library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)

## for data munging
library(dplyr)
library(parallel)

## mantel test
library(vegan)

## maps
library(fields)

##----------------
## first, i'll import the data
## and do some munging
##----------------
##----
## import mapping file
##----
map <- read.csv('./data/map.txt',
                strip.white = TRUE,
                sep = '\t',
                header = TRUE,
                row.names = 1)

## subset mapping file by USA samples
map <- subset(map, country == 'USA')

## subset mapping file to contiguous US
map <- subset(map, state != 'AK' & state != 'HI')

## only consider folks who haven't had 
## antibiotics in the last year
map <- map[grepl('I have not taken antibiotics in the past year.', map$antibiotic_history), ]

## retrieve IDs for samples with coordinates
have.coords <- rownames(map)[which(!grepl('Unspecified', map$latitude) &
                                   !grepl('Missing: Not collected', map$latitude) &
                                   !grepl('Unspecified', map$longitude) &
                                   !grepl('Missing: Not collected', map$longitude)) %>%
                             unique()]

## subset map and S by dropping samples with no coordinates
map <- map[have.coords, ]

## extract coordinates and dump in new array. coerce into numeric
S <- data.frame('longitude' = as.numeric(as.character(map$longitude)),
                'latitude' = as.numeric(as.character(map$latitude)))
rownames(S) <- rownames(map) ## add sample IDs

##----
## import bray-curtis distance matrix
##----
bray <- read.csv('./data/bray_curtis.txt',
                 row.names = 1,
                 header = TRUE,
                 sep = '\t') %>%
  as.matrix()

## fix column names (character is automatically appended to numeric column names)
colnames(bray) <- gsub('X', '', colnames(bray))

## match samples in bray matrix to those retained in mapping file
bray <- bray[rownames(map), rownames(map)]

##----------------
## now I'll generate a spatial distance matrix
## and prep data for distance-decay
##----------------
## compute great circle distance between samples
gc <- rdist.earth(S, miles = FALSE, R = 6371) ## earth's radius km
rownames(gc) <- colnames(gc) <- rownames(S) ## add row/col names

##-----
## take upper triangle of each matrix
##-----
## bray matrix
tri.bray <- bray
diag(tri.bray) <- NA
tri.bray[lower.tri(tri.bray)] <- NA

## gc matrix
tri.gc <- gc
diag(tri.gc) <- NA
tri.gc[lower.tri(tri.gc)] <- NA

## merge into melted data array
gg.dat <- data.frame('bray' = na.omit(c(tri.bray)), 'gc' = na.omit(c(tri.gc)))

##----------------
## analysis starts here
##----------------
## calculate the number of cores for parallelization
no.cores <- detectCores()

## window radii for mantels
windows <- c(50, 100, 500, 1000, 2500, 4500)

## empty data frame to populate
mant.results <- data.frame('window' = windows, 'r' = NA, 'p' = NA)

##----
## record system time
##----
start.time <- Sys.time()

## loop through and fit for each window value
for(f in 1:length(windows)){

  ## subsample data by gc distance
  k <- windows[f]
  gg.dat.sub <- subset(gg.dat, gc <= k)
  
  ##----------------
  ## analyses start here
  ##----------------
  ## subset gc distance matrix by pairwise comps that that fall in window
  gc.sub <- gc
  drop.em <- which(gc.sub > k) ## index elements to exclude
  gc.sub[drop.em] <- NA
  gc.sub <- gc.sub %>% as.dist()
  
  ## subset bray distance matrix by index above
  bray.sub <- bray
  bray.sub[drop.em] <- NA
  bray.sub <- bray.sub %>% as.dist()
  
  ##------
  ## mantel test
  ##-----
  ## initiate cluster
  cl <- makeCluster(no.cores)
  
  ## mantel
  mant <- mantel(bray.sub, gc.sub, permutations = 999, na.rm = TRUE, parallel = cl)
  mant.results$r[f] <- mant$statistic
  mant.results$p[f] <- mant$signif
  
  ## stop cluster
  stopCluster(cl)
 
  ## drop a line
  cat(f / length(windows) %>% round() * 100, '% done\n')
  
}

## adjust p values with fdr
mant.results$p.adj <- mant.results$p %>% p.adjust('fdr')

## save output
write.csv(mant.results, './data/mantel-results-bray.csv')

## print end time
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

##----------------
## plotting starts here
##----------------
## add binary variable to mantel results if p is significant
mant.results$p.sig <- rep(0)
mant.results$p.sig[mant.results$p < 0.05] <- 1 %>% factor()
mant.results$p.sig <- mant.results$p.sig %>% as.factor()

## dummy group variable
mant.results$group <- rep(1)

##----
## correlogram
##----
gg.corr <- ggplot(data = mant.results, aes(x = window, y = r, fill = p.sig, group = group)) +
  theme_classic() +
  scale_y_continuous(limits = c(-0.002, 0.04)) +
  ylab('Mantel r') +
  xlab('Distance (km)') +
  geom_hline(yintercept = 0, size = 0.5, linetype = 2) +
  geom_line(size = 0.5, linetype = 3) +
  geom_point(shape = 21, size = 4.5, stroke = 0.5) +
  scale_fill_manual(values = c('#2b8cbe', '#cb181d')) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9.5))
gg.corr

## save figure
pdf('./figs/corr-plot-bray.pdf', width = 7, height = 6)
gg.corr

## null device
dev.off()

##----
## distance-decay figures
##----
## subsample data by gc distance
k <- 100
gg.dat.sub <- subset(gg.dat, gc <= k)

## color palette
get.palette <- colorRampPalette(brewer.pal(9, 'YlGnBu'))

## 100 km window ggplot
dd.100 <- ggplot(gg.dat.sub, aes(x = gc, y = bray, fill = ..count..)) +
  theme_classic() +
  ylab('Bray-Curtis Distance') +
  xlab('Great Circle Distance (km)') +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  stat_bin_hex(bins = 18) +
  # stat_bin2d(bins = 20, drop = TRUE) +
  stat_smooth(data = gg.dat.sub, aes(x = gc, y = bray),
              size = 0.5, linetype = 2, colour = 'black',
              method = 'glm', method.args = list(family = 'quasibinomial'),
              level = 0, inherit.aes = FALSE) +
  scale_fill_gradientn(colors = get.palette(11),
                       na.value = 'white',
                       trans = 'sqrt',
                       name = 'Frequency',
                       breaks = pretty_breaks(n = 4)) +
  theme(legend.text = element_text(size = 7),
        legend.title = element_text(face = 'bold', size = 8),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9.5))
dd.100

## create rectangle for inset plot
rect <- data.frame(xmin = 0, xmax = 100, ymin = -Inf, ymax = Inf)

## 4500 km window ggplot
dd.4500 <- ggplot(gg.dat, aes(x = gc, y = bray, fill = ..count..)) +
  theme_classic() +
  ylab('') +
  xlab('') +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  stat_bin_hex(bins = 30) +
  # stat_bin2d(bins = 20, drop = TRUE) +
  stat_smooth(data = gg.dat, aes(x = gc, y = bray),
              size = 0.5, linetype = 2, colour = 'black',
              method = 'glm', method.args = list(family = 'quasibinomial'),
              level = 0, inherit.aes = FALSE) +
  geom_rect(data = rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = 'grey20',
            size = 0,
            alpha = 0.35,
            inherit.aes = FALSE) +
  scale_fill_gradientn(colors = get.palette(11),
                       na.value = 'white',
                       trans = 'sqrt',
                       name = 'Frequency',
                       breaks = pretty_breaks(n = 6)) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 7),
        axis.title = element_blank())
dd.4500

##--------
## inset plot
##--------
grid.newpage()

## save figure
pdf('./figs/inset-plot-bray.pdf', width = 7, height = 6)
v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2 <- viewport(width = 0.3, height = 0.3, x = 0.72, y = 0.28) #plot area for the inset map
print(dd.100, vp = v1)
print(dd.4500, vp = v2)

## null device
dev.off()




