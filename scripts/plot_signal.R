#!/usr/bin/env Rscript
# -*- mode: R -*-

# plot_signal: produces a plot of ATAC-seq signal around features,
# using the output of the measure_signal script

options(stringsAsFactors = F)
library(entropy)
library(ggplot2)
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggplot2"))

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

argv <- commandArgs(TRUE)

signal_file = argv[1]
title = argv[2]
output_filename <- argv[3]
plot_filename <- argv[4]

## Step 1: Read and pre-process data
# Read data table
d <- read.delim(signal_file)
# Trim data to our region of interest
d <- d[d$FragmentSize >= 36 & d$FragmentSize <= 400,]
# Transform frag sizes in factors to avoid 0-count sizes not being included
d2 <- d
d2$FragmentSize <- factor(d2$FragmentSize, levels = 36:400)
# Make counts matrix
d2 <- table(d2$FragmentMidpointDistanceToFeatureMidpoint, d2$FragmentSize)
d2 <- as.data.frame.matrix(d2)


## Step 2: Make compressed matrices
# Compress matrix in 5 and 11 bp windows to boost signal
comp_mat <- matrix(nrow = 91, ncol = 365)
for(i in 0:90){
    j <- 11*i
    comp_mat[i+1,] <- as.numeric(apply(d2[(j+1):(j+11),], 2, sum))
}; rm(i,j)
comp_mat2 <- matrix(nrow = 91, ncol = 73)
for(i in 0:72){
    j <- 5*i
    comp_mat2[,i+1] <- as.numeric(apply(comp_mat[,(j+1):(j+5)], 1, sum))
}; rm(i,j)

# Compress matrix in 11 bp windows and 3 fragment sizes bins
sizes <- list(s = c(36,149), m = c(150, 325), l = c(325,400))
comp_mat_bins <- matrix(nrow = 91, ncol = 3)
for(i in 1:3){
    start = sizes[[i]][1]-35
    end =   sizes[[i]][2]-35
    comp_mat_bins[,i] <- as.numeric(apply(comp_mat[,(start):(end)], 1, sum))/(end-start)
}; rm(i, sizes, start, end)


## Step 3: Calculate entropies
# Normalized entropy function
normentropy <- function(x){
    maxentro <- entropy(rep(1,length(x)))
    normentro <- ( maxentro - entropy(x) ) / maxentro
    return(normentro)
}
# Calculate entropy
entro_bins100 <- apply(comp_mat_bins[36:55,],1,normentropy)/20
entro_bins500 <- apply(comp_mat_bins[c(1:35,56:91),],1,normentropy)/71
entro_full <- sum(apply(comp_mat_bins,1,normentropy))
entro100 <- sum(entro_bins100)
entro500 <- sum(entro_bins500)
entro_ratio <- entro100/entro500

## Step 4: Output data
write.table(matrix(c(entro100, entro500, entro_ratio, entro_full),nrow = 1), row.names = F,
            col.names = F, file = output_filename, sep = '\t')
cat(paste("Entropies written to", output_filename, "\n"))

## Step 5: Plot information content per region

posinfo  <- data.frame(info = apply(comp_mat2, 1, normentropy), pos = -501+c(1:91)*11)
fraginfo <- data.frame(info = apply(comp_mat2, 2, normentropy), size = 35+c(1:73)*5)

p1 <- ggplot(posinfo, aes(x = pos , y= info)) + geom_line() +
    theme_bw() + theme(
                     strip.background = element_rect(fill="gray90", colour=FALSE),
                     panel.border = element_rect(colour="gray90"),
                     plot.title=element_text(vjust= 1, face = 'bold'),
                     axis.title.x=element_text(vjust = -0.5, size = 7),
                     axis.title.y=element_text(vjust = 1, size = 7),
                     axis.text.x = element_text(size = 5),
                     axis.text.y = element_text(size = 5)
                 ) +
    labs(x = 'Fragment midpoint position relative to feature midpoint',
         y = 'Normalized information') +
    ylim(0, 0.28)

p2 <- ggplot(fraginfo, aes(x = size, y = info)) + geom_line() + coord_flip() +
    theme_bw() + theme(
                     strip.background = element_rect(fill="gray90", colour=FALSE),
                     panel.border = element_rect(colour="gray90"),
                     plot.title=element_text(vjust= 1, face = 'bold'),
                     axis.title.x=element_text(vjust = -0.5, size = 7),
                     axis.title.y=element_text(vjust = 1, size = 7),
                     axis.text.x = element_text(size = 5),
                     axis.text.y = element_text(size = 5)
                 ) +
    labs(x = 'Fragment size (bp)',
         y = 'Normalized information',
         title = '') +
    ylim(0, 0.2)

if(nrow(d) >= 250000){
    d_to_plot <- d[sample(nrow(d), 250000, replace = F),]
} else {d_to_plot <- d}

p <- ggplot(d_to_plot, aes(x=FragmentMidpointDistanceToFeatureMidpoint, y=FragmentSize,
                           colour=FragmentPositionRelativeToFeatureMidpoint))
p <- p + geom_point(size=0.025, alpha=0.05)
p <- p + ggtitle(title)
p <- p + theme_bw(base_size=10)
p <- p + theme(
             strip.background = element_rect(fill="gray90", colour=FALSE),
             panel.border = element_rect(colour="gray90"),
             panel.margin = unit(1, "lines"),
             plot.title=element_text(vjust=1, hjust = 0.5),
             axis.title.y=element_text(vjust=1),
             axis.title.x=element_blank()
         ) +
    labs(y = '')
p <- p + scale_colour_manual(values = c("blue", "red", "black"), guide = F)


p0 <- ggplot() + theme_bw() + theme(strip.background = element_rect(fill="white", colour=FALSE),
                                    panel.border = element_rect(colour="white"))


layout <- matrix(c(1,2,2,
                   1,2,2,
                   4,3,3), ncol=3, byrow = T)
pdf(file = plot_filename, width = 6, height = 5)
multiplot(layout = layout, p2, p, p1, p0)
dev.off()
cat(paste("Plot written to", plot_filename, "\n"))
