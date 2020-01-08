### packages

stopifnot(
	require(optparse),
	require(gplots),
	require(grid),
	require(gridExtra),
	require(reshape2),
	require(ggplot2),
	require(RSQLite),
	require(stats)
	)

### command line args
option_list <- list(
	    make_option(c("--sample", "-s"), default=NULL, type="character", action="store", help="Name of sample"),
	    make_option(c("--database", "-d"), default=NULL, type="character", action="store", help="Path to database")
)

opts <- parse_args(OptionParser(option_list=option_list))

theme_science <- function(base_size=18, base_family="helvetica"){
                           (theme_set(theme_bw(base_size=base_size)) +
		           theme(plot.title = element_text(size=base_size+2, hjust=0.5),
			   panel.grid.major = element_blank(),
			   panel.grid.minor = element_blank(),
			   panel.background = element_rect(colour='black', fill=NA),
			   strip.background = element_rect(fill=NA, colour=NA),
			   strip.text = element_text(size = base_size+1),
			   text = element_text(size = base_size),
			   axis.title.y = element_text(angle=90, vjust=2, size=base_size+2),
			   axis.title.x = element_text(vjust=-0.2, size=base_size+2),
			   plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines")
			   ))
}

analyse_counts <- function(sample, db){
	       # get data from db and normalise for input into SE analysis

	       table_a <- paste0(sample, "_norm_counts")
	       table_b <- paste0("ENHANCERSmerged_", sample)

	       con = dbConnect(RSQLite::SQLite(), dbname=db)
	       statement <- paste('select a.contig, a.start, a.end, a.peak_id, b.peak_no as no_atac_peaks,
	                    		  a.raw_counts, a.RPM, a.RPM_width_norm, a.width, a.feature from',
			   table_a,
			   'a inner join',
			   table_b,
			   'b on a.peak_id = b.peak_id and  a.feature = "mergedEnhancerPeaks" where feature = "mergedEnhancerPeaks" ')

		df = dbGetQuery(con,  statement)

		# normalise data
		df$ynorm <- df$RPM - sum(df$RPM) # signal - total RPM in input
		df$ynorm <- (df$ynorm - min(df$ynorm) ) / (max(df$ynorm) - min(df$ynorm) ) # min-max score
		df <- df[order(df$ynorm), ] # sort by increasing score
		rownames(df) <- 1:nrow(df) # reset index
		df$x <- rownames(df) # number enhancers in order of increasing ynorm

		return(df)
}

calculate_cutoff <- function(df, sample, drawPlot=TRUE){
		 ### Function from ROSE: https://bitbucket.org/young_computation/rose/src/master/ROSE_callSuper.R ###

		 #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)

		 inputVector <- sort(df$ynorm)

		 inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero

		 slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.

		 #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
		 xPt <- floor(optimize(numPts_below_line, lower=1, upper=length(inputVector), myVector=inputVector, slope=slope)$minimum)

		 y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.

		 if(drawPlot){  #if TRUE, draw the plot

		 		df$se <- "darkgray"
				df$se[df$ynorm > y_cutoff] <- "red3"

				p <- ggplot(df, aes(y=ynorm, x=1:length(ynorm))) +
				  geom_line() +
				  geom_point(colour=df$se) +
				  geom_vline(xintercept=xPt, linetype="dashed", colour="black")  +
				  labs(title="Super Enhancer Cutoff", y="Cumulative ATAC signal", x="Ranked Peaks") +
				  theme_science()

				  fname <- paste0("superenhancer.dir/", sample, ".se_cutoff.pdf")

				  ggsave(fname, device="pdf")
}

		return(y_cutoff)
}

numPts_below_line <- function(myVector,slope,x){
    yPt <- myVector[x]
    b <- yPt-(slope*x)
    xPts <- 1:length(myVector)
    return(sum(myVector<=(xPts*slope+b)))
}

get_enhancers <- function(df, se_cutoff, sample_name){
	      # apply se cutoff to ranked enhancers - yield enhancers and super enhancers

	      df$enhancer <- "enhancer"
	      df$enhancer[df$ynorm > se_cutoff] <- "super_enhancer"

	      se <- subset(df, enhancer=="super_enhancer")
	      e <- subset(df, enhancer=="enhancer")

	      se_name = paste0("superenhancer.dir/", sample_name, ".superenhancers.bed")
	      e_name = paste0("superenhancer.dir/", sample_name, ".enhancers.bed")

	      write.table(se, file=se_name, sep="\t", col.names=F, row.names=F, quote=F)
	      write.table(e, file=e_name, sep="\t", col.names=F, row.names=F, quote=F)

	      return(df)
}

qc_plot <- function(df, sample){
	# plot boxplots

	p1 <- ggplot(df, aes(y=width/1000, x=enhancer, fill=enhancer)) +
	   geom_boxplot() +
	   scale_fill_manual(values=c("darkgray", "red3"), guide=F) +
	   theme_science() +
	   labs(y="Enhancer Width (kb)", x="")

	p2 <- ggplot(df, aes(y=RPM, x=enhancer, fill=enhancer)) +
	   geom_boxplot() +
	   scale_fill_manual(values=c("darkgray", "red3"), guide=F) +
	   theme_science() +
	   labs(y="Enhancer signal (RPM)", x="")

	p3 <- ggplot(df, aes(y=no_atac_peaks, x=enhancer, fill=enhancer)) +
	   geom_boxplot() +
	   scale_fill_manual(values=c("darkgray", "red3"), guide=F) +
	   theme_science() +
	   labs(y="No ATAC peaks", x="")

	name = paste0("superenhancer.dir/", sample, ".QCboxplots.pdf")

	pdf(name,
	    width=12,
	    height=4)

	grid.arrange(p1, p2, p3, ncol=3, nrow=1)

	dev.off()
}

df <- analyse_counts(opts$sample, opts$database)
se_cutoff <- calculate_cutoff(df, opts$sample)
df <- get_enhancers(df, se_cutoff, opts$sample)
qc_plot(df, opts$sample)