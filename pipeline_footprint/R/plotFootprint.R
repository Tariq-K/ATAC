### packages

stopifnot(
	require(optparse),
	require(ggplot2),
	require(wesanderson)
	)

### command line args
option_list <- list(
	    make_option(c("--infiles", "-i"), default=NULL, help="BED file of ATAC cutsites, or comma seperated list of files"),
   	    make_option(c("--outfile", "-o"), default=NULL, help="Filename for plot (*.png)"),
	    make_option(c("--bandwidth", "-b"), default=3, help="Bandwidth for smoothing"),
   	    make_option(c("--labels", "-l"), default="", help="Labels for plot"),
       	    make_option(c("--title", "-t"), default="", help="Title for plot"),
	    make_option(c("--xlims"), default=250, type="integer", help="Window limits for plot, e.g. -250,250"),
	    make_option(c("--motif_width"), default=0, help="Optional, width of motif for plot annotation"),
	    make_option(c("--unsmoothed"), default=FALSE, help="TRUE/FALSE, show unsmoothed data"),
   	    make_option(c("--wes_palette"), default="Royal1", help="Named Wes Anderson palette to use")
)

opts <- parse_args(OptionParser(option_list=option_list))


science_theme <- theme(
  panel.grid.major = element_line(size = 0.5, color = "grey"),
  panel.grid.minor = element_blank(),
  text = element_text(size = 14),
  axis.line = element_line(color="black", size = 0.7),
  axis.line.x = element_line(color="black", size = 0.7),
  axis.line.y = element_line(color="black", size = 0.7)
  # plot.margin = unit(c(0.7,0.7,0.7,0.7), "lines")
)


SmoothProfile <- function(profile, bandwidth=5){
    # Smooth a profile given the specified badnwidth and a normal kernel
    #
    # Args:
    #   profile: smoothed profile
    #   badnwidth: bandwidth for normal kernel smoothing (default=5)
    #
    # Returns:
    #   smoothed profile

    profile <- ksmooth(c(1:length(profile)), profile, kernel="normal", bandwidth=bandwidth)$y

    return(profile)

}

palette <- wes_palette(opts$wes_palette)



plotFootprint <- function(footprint, b=3, title="", labels="", xlims=c(-250,250), motif_width=FALSE, plot_unsmoothed=FALSE,  palette=palette){

    # set delims
    xmin = as.integer(-xlims)
    xmax = as.integer(xlims)

    # subset on window to plot
    footprint <- subset(footprint, position >=xmin & position <= xmax) # subset

    # normalise cut frequency to window
    footprint[["motif"]] <- sapply(footprint[["motif"]], function(x){(x/sum(footprint[["motif"]])*100)})

    # smooth signal
    footprint[["smooth"]] <- SmoothProfile(footprint[["motif"]], bandwidth=b)

    # plot
    p1 <- ggplot(footprint, aes(x=position, y=smooth)) + 
            geom_line(size=1, color=palette[2]) +
            labs(x="Relative position (b.p.)", y="Relative cut frequency", title=title) + 
            coord_cartesian(ylim=c(0, max(footprint[["motif"]])), xlim=c(xmin, xmax), expand=FALSE) + 
            theme_bw() + science_theme + 
            theme(panel.grid = element_blank(), text = element_text(size=20)) 
    
    if (motif_width==0){
        p1 <- p1 + geom_vline(xintercept=0, lty="dashed", col="black")
	} else {
        p1 <- p1 + geom_vline(xintercept = c(0-motif_width/2, 0+motif_width/2), lty="dashed", col="black")
    }
    
    if (plot_unsmoothed==TRUE){
        p1 <- p1 + geom_line(aes(x=position, y=motif), colour=palette[1]) + 
                   geom_line(aes(x=position, y=smooth), colour=palette[2], size=1) 
    }

    return(p1)
}

plot2Footprints <- function(footprints, labels="", title="", b=3, xlims=250, motif_width=FALSE, plot_unsmoothed=FALSE,  palette=palette){

    # first subset on window to plot
    xmin = as.integer(-xlims) # min delim
    xmax = as.integer(xlims) # max delim

    n = 0 
    for (f in footprints){
    	n <- n +1
    	df <- read.table(gzfile(f), header=TRUE, sep=",")
	df <- subset(df, position >=xmin & position <= xmax) # subset	
	df[["motif"]] <- sapply(df[["motif"]], function(x){(x/sum(df[["motif"]])*100)}) # then normalise signal to window
    	df[["smooth"]] <- SmoothProfile(df[["motif"]], bandwidth=b) # now smooth plots
    	df[["label"]] <- strsplit(labels, ",")[[1]][n] # add labels
	
	if (n==1){
	   footprint <- df
	} else {
	   footprint <- rbind(footprint, df) # cat dfs
	}
    }

    # plot
    p1 <- ggplot(footprint, aes(x=position, y=smooth, colour=label)) + 
            geom_line(size=1) +
            scale_colour_manual(values=palette, name ="") +
            labs(x="Relative position (b.p.)", y="Relative cut frequency", title=title) + 
            coord_cartesian(ylim=c(0, max(footprint[["motif"]])), expand=FALSE) + 
            theme_bw() + science_theme +
            theme(panel.grid = element_blank(), text = element_text(size=20),
                  legend.position="bottom", legend.direction="horizontal")
    
    if (motif_width==0){
        p1 <- p1 + geom_vline(xintercept=0, lty="dashed", col="black")
	} else {
        p1 <- p1 + geom_vline(xintercept = c(0-motif_width/2, 0+motif_width/2), lty="dashed", col="black")
    }
    
    if (plot_unsmoothed==TRUE){
        p1 <- p1 + geom_line(aes(x=position, y=motif),  lty="dashed") + 
                   geom_line(aes(x=position, y=smooth)) 
    }
    
    return(p1)
}

# read infiles
infiles <- strsplit(opts$infiles, ",")[[1]]

# plot 1 or 2 footprints
if (length(infiles) > 1){

   p <- plot2Footprints(infiles, labels=opts$labels, title=opts$title, xlims=opts$xlims, palette=palette,
     				  b=opts$bandwidth, motif_width=opts$motif_width, plot_unsmoothed=opts$unsmoothed)

  } else {

  infile1 <- read.table(gzfile(infiles[1]), header=TRUE, sep=",")

  p <- plotFootprint(infile1, labels=opts$labels, title=opts$title, b=opts$bandwidth, xlims=opts$xlims, palette=palette,
       			      motif_width=opts$motif_width, plot_unsmoothed=opts$unsmoothed)

  
}

# save plot
ggsave(p, file=opts$outfile)

