#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
          args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      mergeR.R - Merge 2 files by column
 
      Arguments:
      --i=someValue            - input directory with bed files
      --from=someNumber        - \"from\" parameter of plotTracks()
      --to=someNumber          - \"to\" parameter of plotTracks()
      --b=someValue            - bedGraph
      --p=someValue            - protein Name
      --out=someValue          - output file
      --group=someValue        - group list
      --help                   - print this Help

      Example:
      ./covplotteR.R --i=\"inputdir/\" --f=\"0\" --to=\"1000\" --b=\"input.bedGraph\" --p=\"ProteinName\" --out=\"output.png\" --group=\"ATH,GMA\"
      
      Daniel Guariz Pinheiro
      FCAV/UNESP - Univ Estadual Paulista
      dgpinheiro@gmail.com
      \n\n")

  q(save="no")
}


## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(is.null(argsL[['i']])) {
        sink(stderr())
        cat("\nERROR: Missing input directory with .bed files !\n\n")
        sink()
        q(save="no")
}

if(is.null(argsL[['b']])) {
        sink(stderr())
        cat("\nERROR: Missing input .bedGraph file !\n\n")
        sink()
        q(save="no")
}

if(is.null(argsL[['from']])) {
argsL[['from']] = NULL
}

if(is.null(argsL[['to']])) {
argsL[['to']] = NULL
}

if(is.null(argsL[['p']])) {
argsL[['p']] = 'Some protein';
}

if(is.null(argsL[['out']])) {
        sink(stderr())
        cat("\nERROR: Missing output file !\n\n")
        sink()
        q(save="no")
}
if(is.null(argsL[['group']])) {
        sink(stderr())
        cat("\nERROR: Missing group list !\n\n")
        sink()
        q(save="no")
}

suppressPackageStartupMessages( library("Gviz") )
suppressPackageStartupMessages( library("rtracklayer") )
suppressPackageStartupMessages( library("RColorBrewer") )

options(ucscChromosomeNames=FALSE)

anntrack <- list()

annheight <- 0

anntrack[['gtrack']] <- GenomeAxisTrack()

displayPars(anntrack[['gtrack']]) <- list(size=100)
annheight <- annheight+100

anntrack[['ctrack']] <- DataTrack(range=import.bedGraph(argsL[['b']]),
                      type = "histogram",
                      name = "Coverage",
                      fill.histogram="gray", 
                      col.histogram="NA",
                      col.axis="black",
                      col="black",
                      window=-1)

displayPars(anntrack[['ctrack']]) <- list(size=200)
annheight <- annheight+200

grouplst <- unlist(strsplit(argsL[['group']], ','))
npal <- 5
if ( length(grouplst) > npal ) {
	npal <- length(grouplst)
}
colorlst <- brewer.pal(n = npal, name = "Set3")

for (i in seq(1,length(grouplst))) {
	grp <- grouplst[i]
	col <- colorlst[i]
	for (inbed in list.files(argsL[['i']], pattern=paste(grp,'.*\\.bed$',sep=""))) {
		id <- gsub('.bed$','',basename(inbed))
		print(paste("Loading", id , '...'))
		
		bedf <- import.bed(paste(argsL[['i']],inbed,sep="/"))
		
		anntrack[[id]] <- AlignmentsTrack(bedf, 
				name=id,
				fill=col,
				stacking = "squish",
				coverageHeight=0.7,
				col.axis="black",
				isPaired=FALSE
				)
		
		displayPars(anntrack[[id]]) <- list(size=350)
		annheight <- annheight+350
	}
}

png(filename=argsL[['out']], width=1024, height=annheight, type="cairo")
plotTracks(anntrack,
           from = as.numeric(as.character(argsL[['from']])),
           to = as.numeric(as.character(argsL[['to']])),
           shape="box"
           )
graphics.off()
