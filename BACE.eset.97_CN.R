#	BACE
#	Breakthrough Array CGH and Expression script
#	A suite of basic functions to process and integrate expression array data and arrayCGH
#	Written to accompany in-house spotted microarrays used at the Breakthrough Research Centre
#	Written and tested under OSX10.8.5
#	Current R version 3.0.2
#	
#	Alan Mackay
#
#	Comments, suggestions and requests for support are welcome
#	Please contact alan.mackay@icr.ac.uk 
#
#	Libraries required
#

library(Biobase)
library(affy)
library(annotate)
library(limma)
library(impute)
library(lattice)
library(grid)
library(marray)
library(RColorBrewer)

dietCGH <- function(design.file=NULL, subtractBG=FALSE, flag.values = c(-50,-75,-100), byDye=T, fdata.file = "Ann32K.txt", rhm=F, MAD=2)
#	Function to wrap a number of basic functions
#	for the default loading and pre-processing of aCGH data
#	fdata.file - featureData annotation file
#	rhm - repHybMean averages dye swap pairs if you have them
{
	require(Biobase)
	design <- readDesign(design.file)
	if(!is.element(fdata.file, list.files())) stop("fdata.file not found")
	fdata <- read.table(file=fdata.file, header=T, sep="\t", stringsAsFactors=F)
	if(is.null(fdata$probeID)) stop("fdata.file should contain \"probeID\"")
	if(is.null(fdata$chrom)|is.null(fdata$start)|is.null(fdata$end)) stop("fdata.file should contain \"chrom\" \"start\" and \"end\"")
	rm(fdata)
	raw <- getRawData(design.file=design.file, subtractBG=subtractBG, flag.values=flag.values, byDye=byDye)
	cgh.norm <- loessFitBlock(raw)
	if(identical(rhm,T)) cgh.norm <- repHybMean(cgh.norm, design$sampleNames)
	cgh.norm <- annotateCGH(cgh.norm, fdata.file=fdata.file)
	cgh.norm <- removeOutliersByMAD(cgh.norm, MAD=MAD)
	exprs(cgh.norm) <- imputeCGH(cgh.norm, span=2)
	cgh.norm
}


dietExpression <- function(design.file=NULL, subtractBG=F, flag.values = c(-50,-75,-100), byDye=T, fdata.file = "Hs20Kv3.txt", rhm=F)
#	Function to wrap a number of basic functions
#	for the default loading and pre-processing of expression data
#	ann.file - anotation file
#	rhm - repHybMean averages dye swap pairs if you have them
#	flag.limit - specifies the number of NAs allowed per probe
{
	require(Biobase)
	design <- readDesign(design.file)
	raw <- getRawData(design.file=design.file, subtractBG=subtractBG, flag.values=flag.values, byDye=byDye)
	cgh.raw <- annotateExpression(raw, fdata.file=fdata.file)
	
	norm <- loessFitMA(ma)
	if(identical(rhm,T)) norm <- repHybMean(norm, design$sampleNames)

	norm
}



readDesign <- function(design.file=NULL)
#	Reads design table
#	Design table is a tab delimited text file
#	file, sampleNames, dye and barcode as minimum information
{
	if(is.null(design.file)) {design.file <- file.choose()}

	design <- read.delim(design.file, header=T, sep="\t", stringsAsFactors=F)
	
	if(is.null(design$file)) stop("Please specify \"file\" in design file")
	if(is.null(design$sampleNames)) stop("Please specify \"sampleNames\" in design file")
	if(is.null(design$dye)) stop("Please specify \"dye\" in design file")

	design$file <- as.character(design$file)
	design$sampleNames <- as.character(design$sampleNames)
	design$dye <- as.integer(design$dye)
	if(!is.null(design$barcode)) design$barcode <- as.integer(design$barcode)

	design
}


getRawData <- function(design.file=NULL, subtractBG=FALSE, flag.values=c(-50,-75,-100), byDye=T)
{
	require(Biobase)
	
	commands <- data.frame(commands.history="Analysis script: BACE", stringsAsFactors=F)
	commands <- rbind(commands,date())
	commands <- rbind(commands,deparse(match.call()))
	
	if(is.null(design.file)){design.file <- file.choose()}
	cat("Design file chosen was \"",basename(design.file),"\"\n", sep="")
	
	design <- read.delim(design.file, header=T, sep="\t", stringsAsFactors=F)
	design$file <- as.character(design$file)
	design$sampleNames <- as.character(design$sampleNames)
	design$dye <- as.integer(design$dye)
#	if(sum(duplicated(design$sampleNames)) > 0) stop("sampleNames must be unique")
	
	commands <- rbind(commands,paste("Design file =", basename(design.file)))
	
	cat("Checking directory for results files...")
	directory.list <- list.files()
	gpr.list <- directory.list[grep(".*gpr.*",directory.list)]
	missing <- design$file[which(!is.element(design$file,gpr.list))]
	if(length(missing) != 0){
	cat(paste(missing, collapse="\n"),"\n")
	stop("These gpr files in your design table are not in your directory\n")
	}else{cat("OK\n")}	
	
	cat("Reading results files\n")
	gpr <- getGenepix(design$file[1])
	
	cat("Reading results file 1",design$file[1],nrow(gpr$raw),"probes\n")
	Description <- gpr$raw$Name
	probeID <- gpr$raw$ID
	block <- gpr$raw$Block
	column <- gpr$raw$Column
	row <- gpr$raw$Row
	cy3 <- cy5 <- flags <-matrix(0,nrow(gpr$raw),nrow(design))
	
	sampleNames.eset <- make.unique(design$sampleNames)
	featureNames.eset <- make.unique(probeID)

	if(subtractBG==F)
	{
		cy5[,1] <- gpr$raw$F5md
		cy3[,1] <- gpr$raw$F3md
	}
	if(subtractBG==T)
	{
		cy5[,1] <- gpr$raw$F5md.bg
		cy3[,1] <- gpr$raw$F3md.bg
	}
	X <- gpr$raw$X; Y <- gpr$raw$Y
	
	flags[,1] <- gpr$raw$Flags
	flag.50 <- rep("", nrow(design))
	flag.75 <- rep("", nrow(design))
	flag.100 <- rep("", nrow(design))
	flag.50[1] <- 100*length(gpr$raw$Flags[which(gpr$raw$Flags == -50)])/length(gpr$raw$Flags)
	flag.75[1] <- 100*length(gpr$raw$Flags[which(gpr$raw$Flags == -75)])/length(gpr$raw$Flags)
	flag.100[1] <- 100*length(gpr$raw$Flags[which(gpr$raw$Flags == -100)])/length(gpr$raw$Flags)
	
	B635.1SD <- rep("", nrow(design))
	B532.1SD <- rep("", nrow(design))
	B635.1SD[1] <- 100*length(gpr$raw$B635.1SD[which(gpr$raw$B635.1SD > 75)])/nrow(gpr$raw)
	B532.1SD[1] <- 100*length(gpr$raw$B532.1SD[which(gpr$raw$B532.1SD > 75)])/nrow(gpr$raw)
	B635.2SD <- rep("", nrow(design))
	B532.2SD <- rep("", nrow(design))
	B635.2SD[1] <- 100*length(gpr$raw$B635.2SD[which(gpr$raw$B635.2SD > 75)])/nrow(gpr$raw)
	B532.2SD[1] <- 100*length(gpr$raw$B532.2SD[which(gpr$raw$B532.2SD > 75)])/nrow(gpr$raw)
	gal.file <- rep("",nrow(design))
	gal.file[1] <- gpr$gal.file
	
	gpr.version <- rep(" ",nrow(design))
	gpr.version[1] <- gpr$gpr.version
	
	image <- rep("",nrow(design))
	if(is.null(gpr$image)){image[1] <- NA}else{image[1] <- as.character(gpr$image)}

	SNR532 <- rep("",nrow(design))
	SNR635 <- rep("",nrow(design))
	if(!is.null(gpr$raw$SNR532))
	{
	SNR532[1] <- median(as.numeric(gpr$raw$SNR532), na.rm=T)
	SNR635[1] <- median(as.numeric(gpr$raw$SNR635), na.rm=T)
	}
	
	B532 <- rep("",nrow(design))
	B635 <- rep("",nrow(design))
	B532[1] <- median(as.numeric(gpr$raw$B532md), na.rm=T)
	B635[1] <- median(as.numeric(gpr$raw$B635md), na.rm=T)
	
	spot.size <- rep("",nrow(design))
	spot.size[1] <- median(as.numeric(gpr$raw$Size), na.rm=T)
	
	cy3.pmt <- rep("",nrow(design))
	cy3.pmt[1] <- gpr$cy3.pmt	
	cy5.pmt <- rep("",nrow(design))
	cy5.pmt[1] <- gpr$cy5.pmt
	
	if(nrow(design) >= 2)
	{
		for(i in 2:nrow(design))
		{
	   		gpr <- getGenepix(design$file[i])
	   		probeID.i <- gpr$raw$ID
	   		if(!all(probeID.i == probeID)) warning("Results file ",i," has different probeIDs")           
	   		cat("Reading results file",i,design$file[i],nrow(gpr$raw),"probes\n")
	   		if(subtractBG==F)
			{
				cy5[,i] <- gpr$raw$F5md
				cy3[,i] <- gpr$raw$F3md
			}
			if(subtractBG==T)
			{
				cy5[,i] <- gpr$raw$F5md.bg
				cy3[,i] <- gpr$raw$F3md.bg
			}
	   		flags[,i] <- gpr$raw$Flags
	   		flag.50[i] <- 100*length(gpr$raw$Flags[which(gpr$raw$Flags == -50)])/length(gpr$raw$Flags)
			flag.75[i] <- 100*length(gpr$raw$Flags[which(gpr$raw$Flags == -75)])/length(gpr$raw$Flags)
			flag.100[i] <- 100*length(gpr$raw$Flags[which(gpr$raw$Flags == -100)])/length(gpr$raw$Flags)
			gal.file[i] <- gpr$gal.file
			gpr.version[i] <- gpr$gpr.version
			
			if(is.null(gpr$image)){image[i] <- NA}
			else{image[i] <- as.character(gpr$image)}
			
			if(!is.null(gpr$raw$SNR532))
			{
			SNR532[i] <- median(as.numeric(gpr$raw$SNR532[i]), na.rm=T)
			SNR635[i] <- median(as.numeric(gpr$raw$SNR635[i]), na.rm=T)
			}
			B532[i] <- median(as.numeric(gpr$raw$B532md), na.rm=T)
			B635[i] <- median(as.numeric(gpr$raw$B635md), na.rm=T)
			B635.1SD[i] <- 100*length(gpr$raw$B635.1SD[which(gpr$raw$B635.1SD > 75)])/nrow(gpr$raw)
			B532.1SD[i] <- 100*length(gpr$raw$B532.1SD[which(gpr$raw$B532.1SD > 75)])/nrow(gpr$raw)
			B635.2SD[i] <- 100*length(gpr$raw$B635.2SD[which(gpr$raw$B635.2SD > 75)])/nrow(gpr$raw)
			B532.2SD[i] <- 100*length(gpr$raw$B532.2SD[which(gpr$raw$B532.2SD > 75)])/nrow(gpr$raw)
			spot.size[i] <- median(as.numeric(gpr$raw$Size), na.rm=T)
			cy3.pmt[i] <- gpr$cy3.pmt
			cy5.pmt[i] <- gpr$cy5.pmt
		} 
	}
	
	cat("Assigning flags...")
	commands <- rbind(commands,paste("Flag values", paste(flag.values, collapse=" ")))
	flags[!is.element(flags,flag.values)] <- F
	flags[is.element(flags,flag.values)] <- T
	
	is.na(cy3[flags == 1]) <- T
	is.na(cy5[flags == 1]) <- T
	
	flags[is.na(cy3) | is.na(cy5)] <- T
	flags[(cy3<0) | (cy5<0)] <- T
	cat("Done\n")

	cat("Converting raw data to MA...")
	Lg.cy3 <- logb(cy3,2)
	Lg.cy5 <- logb(cy5,2)
		
	M <- Lg.cy5 - Lg.cy3
	A <- (Lg.cy5 + Lg.cy3)/2
	flags[is.infinite(M)|is.infinite(A)] <- T
	flags[is.na(M)|is.na(A)] <- T
	is.na(M) <- is.na(A) <- flags

	if(byDye == T) M[,which(design$dye == 3)] <- -1.0 * M[,which(design$dye == 3)]
	cat("Done\n")
	
	colnames(M) <- colnames(A) <-colnames(flags) <- sampleNames.eset
	rownames(M) <- rownames(A) <- rownames(flags) <- featureNames.eset
	
#	Write Eset	
#	Make QC data from gprs to write out	
	QC <- data.frame(sampleNames=sampleNames.eset)
	QC$file <- design$file
	QC$dye <- design$dye
	QC$barcode <- design$barcode
	QC$gal.file = gal.file
	QC$gpr.version = gpr.version
	if(!is.null(image))
	{
	QC$image <- image
	}
	QC$cy5.pmt = cy5.pmt
	QC$cy3.pmt = cy3.pmt
	QC$F635 <- apply(cy5,2,median,na.rm=T)
	QC$F532 <- apply(cy3,2,median,na.rm=T)
	QC$B635 <- B635
	QC$B532 <- B532
	QC$SNR532 <- SNR532
	QC$SNR635 <- SNR635
	QC$B635.1SD <- signif(as.numeric(B635.1SD), digits=3)
	QC$B532.1SD <- signif(as.numeric(B532.1SD), digits=3)
	QC$B635.2SD <- signif(as.numeric(B635.2SD), digits=3)
	QC$B532.2SD <- signif(as.numeric(B532.2SD), digits=3)
	QC$spot.size <- spot.size
	QC$flag.50 <- round(as.numeric(flag.50), 2)
	QC$flag.75 <- round(as.numeric(flag.75), 2)
	QC$flag.100 <- round(as.numeric(flag.100), 2)
	write.table(QC, "Genepix.QC.data.xls", row.names=F, sep="\t")
	
	fd <- data.frame(probeID=probeID, block=block, column=column, row=row, X=X, Y=Y, stringsAsFactors=F)
	rownames(fd) <- featureNames.eset
	
	eset.fdata <- new("AnnotatedDataFrame", data=fd)
	varMetadata(eset.fdata)$labelDescription[grep("probeID",names(fd))] <- "probeID"
	varMetadata(eset.fdata)$labelDescription[grep("block",names(fd))] <- "Block (printTip) of the array"
	varMetadata(eset.fdata)$labelDescription[grep("column",names(fd))] <- "Column within each block of the array"
	varMetadata(eset.fdata)$labelDescription[grep("row",names(fd))] <- "Row within each block of the array"
	varMetadata(eset.fdata)$labelDescription[grep("X",names(fd))] <- "X co-ordinate (microns)"
	varMetadata(eset.fdata)$labelDescription[grep("Y",names(fd))] <- "Y co-ordinate (microns)"
	
	eset <- new("ExpressionSet", exprs=M, featureData=eset.fdata)
	eset <- as(eset, "BACE.cgh")
	assayDataElement(eset, "A") <- A
	assayDataElement(eset, "flags") <- flags
	
	if(validObject(eset)) cat ("valid BACE.cgh Object\n")
	eset
}
	

getGenepix <- function(fname)
#	Read input array data from individual GenePix results file (gpr)
{  
	require(Biobase)
	
	header <- scan(file=fname, sep="\t", nlines=28,what=character(),quiet=T, strip.white=T)
	no.to.skip <- as.integer(header[3]) + 2
	gpr.version <- header[5]
	gpr.version <- strsplit(gpr.version, "=")
	gpr.version <- unlist(gpr.version)
	gpr.version <- gpr.version[2]

	gal <- header[grep("gal",header)]
	gal <- lapply(strsplit(gal,"\\\\"),rev)
	gal <- unlist(gal)
	gal <- gal[1]
	
	image <- header[grep("ImageFiles",header)]
	image <- lapply(strsplit(image,"\\\\"),rev)
	image <- unlist(image)
	if(!is.null(image))
	{	
	image <- image[1]
	image <- strsplit(image, " ")
	image <- unlist(image)
	image <- image[1]
	}
	
	pmt <- header[grep("PMT",header)]
	pmt <- lapply(strsplit(pmt,"="),rev)
	pmt <- unlist(pmt)
	pmt <- pmt[1]
	pmt <- strsplit(pmt,"\t")
	pmt <- unlist(pmt)
	
	cy5.pmt <- pmt[1]
	cy3.pmt <- pmt[2]

	header.list <- list(no.to.skip, gal, gpr.version, image, cy3.pmt, cy5.pmt)
	names(header.list) <- c("no.to.skip", "gal.file","gpr.version", "image", "cy3.pmt","cy5.pmt")
	header.list

	gpr <- read.table(file=fname,sep="\t",header=T, check.names=F ,skip=header.list$no.to.skip, stringsAsFactors=F)

	collabs <- NULL

	collabs <- c("Block","Column","Row","Name","ID","X","Y","Dia.","F635 Median","F532 Median","B635 Median", "B532 Median",
	"F635 Mean","F532 Mean","B635 Mean", "B532 Mean",
	"F635 Median - B635","F532 Median - B532","F635 Mean - B635","F532 Mean - B532",
	"% > B635+2SD","% > B532+2SD", "% > B635+1SD","% > B532+1SD","SNR 635","SNR 532","Flags")

	cols.to.read <- match(collabs,names(gpr),nomatch=NA)
	gpr <- gpr[,cols.to.read[which(!is.na(cols.to.read))]]

	if(length(grep("SNR",names(gpr))) ==0)
	{
	colnms <- c("Block","Column","Row","Name","ID","X","Y","Size","F5md","F3md","B635md","B532md",
	"F5mn","F3mn","B635mn","B532mn",
	"F5md.bg", "F3md.bg","F5mn.bg", "F3mn.bg",
	"B635.2SD", "B532.2SD","B635.1SD", "B532.1SD","Flags")	}
	if(length(grep("SNR",names(gpr))) >0)
	{
	colnms <- c("Block","Column","Row","Name","ID","X","Y","Size","F5md","F3md","B635md","B532md",
	"F5mn","F3mn","B635mn","B532mn",
	"F5md.bg", "F3md.bg","F5mn.bg", "F3mn.bg",
	"B635.2SD", "B532.2SD","B635.1SD", "B532.1SD","SNR635", "SNR532","Flags")	}
	names(gpr) <- colnms
	
	return(list(raw=gpr,gal.file=header.list$gal.file, gpr.version=header.list$gpr.version, image=header.list$image, cy3.pmt=header.list$cy3.pmt, cy5.pmt=header.list$cy5.pmt))

}

annotateCGH <- function(eset, fdata.file){

	require(Biobase)
	cat("annotating ExpressionSet CGH\n")
	fdata <- read.table(file=fdata.file, header=T, sep="\t", stringsAsFactors=F)
	if(is.null(fdata$probeID)) stop("fdata.file should contain \"probeID\"")
	if(is.null(fdata$chrom)|is.null(fdata$start)|is.null(fdata$end)) stop("fdata.file should contain \"chrom\" \"start\" and \"end\"")
	if(!is.null(fdata$chrom))
		{
		fdata$chrom[fdata$chrom == "X"] <- 23
		fdata$chrom[fdata$chrom == "Y"] <- 24
		}
	
	if(!is.null(fdata$start)) fdata$chrom <- as.integer(fdata$chrom)
	if(!is.null(fdata$start)) fdata$start <- as.integer(fdata$start)
	if(!is.null(fdata$end)) fdata$end <- as.integer(fdata$end)
	
	eset.fdata <- merge(fData(eset), fdata, by="probeID", all.x=T, all.y=F)
	eset.fdata <- eset.fdata[match(fData(eset)$probeID, eset.fdata$probeID, nomatch=0),]
	rownames(eset.fdata) <- featureNames(eset)
	eset.fdata.adf <- new("AnnotatedDataFrame", data=eset.fdata)
	
	cat(sum(is.element(fdata$probeID,fData(eset)$probeID)), "annotated probes\n")

    varMetadata(eset.fdata.adf)[which(is.element(varLabels(eset.fdata.adf), fvarLabels(eset))),] <- fvarMetadata(eset)
    varMetadata(eset.fdata.adf)$labelDescription[grep("chrom",names(eset.fdata))] <- "Chromosome"
	varMetadata(eset.fdata.adf)$labelDescription[grep("start",names(eset.fdata))] <- "Probe start position (bp)"
	varMetadata(eset.fdata.adf)$labelDescription[grep("end",names(eset.fdata))] <- "Probe end position (bp)"
	varMetadata(eset.fdata.adf)$labelDescription[grep("bac.id",names(eset.fdata))] <- "Probe BAC Clone ID"
	varMetadata(eset.fdata.adf)$labelDescription[grep("cytoband",names(eset.fdata))] <- "Cytoband"
	varMetadata(eset.fdata.adf)$labelDescription[grep("MB",names(eset.fdata))] <- "Mbp genome position"
	featureData(eset) <- eset.fdata.adf
	
	eset <- eset[which(is.element(fData(eset)$probeID, fdata$probeID)),]
	eset <- eset[order(fData(eset)$chrom, fData(eset)$start),]
	eset@annotation <- fdata.file
	cat("Done\n")
	
	cgh <- as(eset, "BACE.cgh")
	cgh
}



annotateExpression <- function(eset, fdata.file){
	
	require(Biobase)
	cat("annotating ExpressionSet\n")
	fdata <- read.table(file=fdata.file, header=T, sep="\t", stringsAsFactors=F, quote="\"")
	if(is.null(fdata$probeID)) stop("fdata.file should contain \"probeID\"")
	
	if(!is.null(fdata$chrom))
		{
		fdata$chrom[fdata$chrom == "X"] <- 23
		fdata$chrom[fdata$chrom == "Y"] <- 24
		}
	if(!is.null(fdata$start)) fdata$start <- as.integer(fdata$start)
	if(!is.null(fdata$end)) fdata$end <- as.integer(fdata$end)
	cat("Done\n")
	
	eset.fdata <- merge(fData(eset), fdata, by="probeID", all.x=T, all.y=F)
	eset.fdata <- eset.fdata[match(fData(eset)$probeID, eset.fdata$probeID, nomatch=0),]
	rownames(eset.fdata) <- featureNames(eset)
	eset.fdata.adf <- new("AnnotatedDataFrame", data=eset.fdata)
	
	cat(sum(is.element(fdata$probeID,fData(eset)$probeID)), "annotated probes\n")

    varMetadata(eset.fdata.adf)[which(is.element(varLabels(eset.fdata.adf), fvarLabels(eset))),] <- fvarMetadata(eset)
    if(is.element("chrom",names(eset.fdata))) varMetadata(eset.fdata.adf)$labelDescription[grep("chrom",names(eset.fdata))] <- "Chromosome"
	if(is.element("start",names(eset.fdata))) varMetadata(eset.fdata.adf)$labelDescription[grep("start",names(eset.fdata))] <- "Probe start position (bp)"
	if(is.element("end",names(eset.fdata))) varMetadata(eset.fdata.adf)$labelDescription[grep("end",names(eset.fdata))] <- "Probe end position (bp)"

	featureData(eset) <- eset.fdata.adf
	
	eset <- eset[which(is.element(fData(eset)$probeID, fdata$probeID)),]

	eset@annotation <- fdata.file
	cat("Done\n")
	
	expression <- new("BACE.exp", exprs=exprs(eset), phenoData=phenoData(eset), featureData=featureData(eset), annotation=fdata.file)
	
	for(i in 1:length(assayDataElementNames(eset))) assayDataElement(expression, assayDataElementNames(eset)[i]) <- assayDataElement(eset, assayDataElementNames(eset)[i])
	expression
}



singleMAplot <- function(eset, case, main=NULL)
{  
	require(Biobase)
	if(is.null(main)) main <- paste("MA plot", sampleNames(eset)[case])
	a <- assayData(eset)$A[,case]
	m <- exprs(eset)[,case]
	ok <- is.finite(m) & is.finite(a)
	aq <- quantile(a[ok],c(0.1,0.25,0.5,0.75,0.9))
	plot(a[ok],m[ok],type="n",xlab="A",ylab="M",main=main,las=1)
 	points(a[ok],m[ok],pch=16,col=1,cex=0.5)
	abline(v=aq,lty=8)
	abline(h=0,lty=1)
}


loessFitBlock <- function(eset)
#	loess normalisation based upon Block i.e. PrintTipLoess
{
	require(Biobase)
	require(limma)	
#	commands <- c(ma$commands,deparse(match.call()))
#	if(class(eset)!="ExpressionSet") stop("Function requires an ExpressionSet object")
	if(is.null(fData(eset)$block)) stop("Expression set feature data should be annotated with \"block\"")
	if(is.null(assayData(eset)$A)) stop("Expression set assay data should contain intensity matrix \"A\"")
	cat("Normalising by Block\n")
	return.order <- (1:nrow(eset))[order(fData(eset)$block)]
	eset <- eset[order(fData(eset)$block),]
		
	for(i in 1:ncol(eset))
	{
		M.blocks <- split(exprs(eset)[,i], fData(eset)$block)
		blocks <- unique(fData(eset)$block)
		fitted.M <- M.blocks

		for(j in 1:length(blocks))
		{
			fitted.M[[j]] <- loessFit(y=exprs(eset)[fData(eset)$block==blocks[j],i],x=assayData(eset)$A[fData(eset)$block==blocks[j],i])$residuals
		}
		exprs(eset)[,i] <- unlist(fitted.M)
		
		cat(i,sampleNames(eset)[i],"\n")
	}

	eset <- eset[order(return.order),]
	cat("Done\n")
	return(eset)	
}



loessFitMA <- function(eset,span=0.5)
# loess normalisation based upon loessFit{limma}
{
	require(Biobase)
	require(limma)
#	if(class(eset)!="ExpressionSet") stop("Function requires an ExpressionSet object")
	if(is.null(assayData(eset)$A)) stop("Expression set assay data should contain intensity matrix \"A\"")
	cat("Normalising with loessFit\n")
#	commands <- c(ma$commands,deparse(match.call()))
	
	for(i in 1:ncol(exprs(eset)))
	{
		exprs(eset)[,i] <- loessFit(y=exprs(eset)[,i],x=assayData(eset)$A[,i],span=span)$residuals
		cat(i,sampleNames(eset)[i],"\n")
	}
#	ma$commands <- commands
	cat("Done\n")
	return(eset)	
}


loessFitXY <- function(eset,span=0.5)
#	loess normalisation based upon XY position
{
	require(Biobase)
#	if(class(eset)!="ExpressionSet") stop("Function requires an ExpressionSet object")
	if(is.null(fData(eset)$X)|is.null(fData(eset)$Y)) stop("Expression set feature data should be annotated with \"X\" and \"Y\"")
	if(is.null(assayData(eset)$A)) stop("Expression set assay data should contain intensity matrix \"A\"")
	cat("Normalising Expression values based upon XY coordinates\n")
#	commands <- c(ma$commands,deparse(match.call()))

	for(i in 1:ncol(exprs(eset)))
	{
		ok <- is.finite(exprs(eset)[,i])
		exprs(eset)[ok,i] <- loessFit(exprs(eset)[ok,i]~fData(eset)$X[ok]*fData(eset)$Y[ok],degree=1,span=span)$residuals
		cat(i,sampleNames(eset)[i],"\n")
	}
#	ma$commands <- commands
	cat("Done\n")
	return(eset)
}

loessNormBlock <- function(eset,span=0.5)
#	loess normalisation based upon Block i.e. PrintTipLoess
{
	require(Biobase)
#	if(class(eset)!="ExpressionSet") stop("Function requires an ExpressionSet object")
	if(is.null(fData(eset)$block)) stop("Expression set feature data should be annotated with \"block\"")
	if(is.null(assayData(eset)$A)) stop("Expression set assay data should contain intensity matrix \"A\"")
	cat("Normalising by Block\n")
	
#	commands <- c(ma$commands,deparse(match.call()))

	for(i in 1:ncol(eset))
	{
#		i <- 1
		ok <- is.finite(exprs(eset)[,i]) & is.finite(assayData(eset)$A[,i])
		exprs(eset)[ok,i] <- loess(exprs(eset)[ok,i] ~ assayData(eset)$A[ok,i]*fData(eset)$block[ok],degree=1,span=span)$residuals
		cat(i,sampleNames(eset)[i],"\n")
	}
#	ma$commands <- commands
	cat("Done\n")
	return(eset)	
}

loessNormMA <- function(ma,span=0.5)
#	Redundant loessFitMA is very much faster and returns the same thing
{
	require(Biobase)
#	if(class(eset)!="ExpressionSet") stop("Function requires an ExpressionSet object")
	if(is.null(assayData(eset)$A)) stop("Expression set assay data should contain intensity matrix \"A\"")

	cat("Normalising with loess\n")
	
#   commands <- c(ma$commands,deparse(match.call()))


   for(i in 1:ncol(eset))
   {
	  ok <- is.finite(exprs(eset)[,i]) & is.finite(assayData(eset)$A[,i])
	  exprs(eset)[ok,i] <- loess(exprs(eset)[ok,i]~assayData(eset)$A[ok,i],degree=1,span=span, na.action=na.omit)$residuals
	  cat(i,sampleNames(eset)[i],"\n")
	}
#	ma$commands <- commands
	cat("Done\n")
	return(eset)
}

loessNormXY <- function(eset,span=0.5)
#	loess normalisation based upon XY position
{
	require(Biobase)
#	if(class(eset)!="ExpressionSet") stop("Function requires an ExpressionSet object")
	if(is.null(fData(eset)$X)|is.null(fData(eset)$Y)) stop("Expression set feature data should be annotated with \"X\" and \"Y\"")
	if(is.null(assayData(eset)$A)) stop("Expression set assay data should contain intensity matrix \"A\"")
	cat("Normalising Expression values based upon XY coordinates\n")
#	commands <- c(ma$commands,deparse(match.call()))

	for(i in 1:ncol(exprs(eset)))
	{
		ok <- is.finite(exprs(eset)[,i])
		exprs(eset)[ok,i] <- loess(exprs(eset)[ok,i]~fData(eset)$X[ok]*fData(eset)$Y[ok],degree=1,span=span)$residuals
		cat(i,sampleNames(eset)[i],"\n")
	}
#	ma$commands <- commands
	cat("Done\n")
	return(eset)
}



readPheno<- function(eset, pheno.file)
#	pheno file must carry column "sampleNames"
#	Reads pheno file into phenoData of ExpressionSet
#	overwrite pheno table in expression set
#	just keep sampleNames
{
	require(Biobase)
	cat("Reading pheno file\n")
	pheno <- read.delim(pheno.file, header=T, sep="\t", stringsAsFactors=T, na.strings=c("", " ", "NA", "N/A", "#N/A"))
	if(is.null(pheno$sampleNames)) stop("Pheno file must be indexed with \"sampleNames\"\n")
	
	cat(length(pheno$sampleNames),"Pheno sampleNames\n")
	pheno <- pheno[match(sampleNames(eset),pheno$sampleNames),]
	
	cat(ncol(eset),"Eset sampleNames\n")
	
	matching.ids <- as.character(pheno$sampleNames[!is.na(match(pheno$sampleNames,sampleNames(eset)))])
	cat(length(matching.ids),"Matching sampleNames\n",matching.ids,"\n")

	eset <- eset[,!is.na(pheno$sampleNames)]
	pheno <- pheno[!is.na(pheno$sampleNames),]
	rownames(pheno) <- as.character(pheno$sampleNames)

	sample.names <- pheno$sampleNames
	pheno <- pheno[,-which(names(pheno) == "sampleNames"), drop=F] 

	pData(eset) <- pheno
	sampleNames(eset) <- sample.names
	cat("Done\n")
	return(eset)
}

spitTables <- function(eset, output.directory, project=NULL, fdata.labels="all")
{
	require(Biobase)
	
	if(is.null(project)) project <- output.directory
	
	if(!is.element(output.directory, list.files())) dir.create(output.directory)
	
	if(identical(fdata.labels,"all")){
		fdata.cols <- fData(eset)
		}else{
			fdata.cols <- fData(eset)[,which(is.element(fvarLabels(eset), fdata.labels)),drop=F]
			names(fdata.cols) <- fvarLabels(eset)[which(is.element(fvarLabels(eset), fdata.labels))]
			}
	
	fdata.labels <- fvarLabels(eset)[which(is.element(fvarLabels(eset), fdata.labels))]
	
	cat("writing featureData table\n")
	probes.table <- data.frame(featureNames=featureNames(eset), fData(eset))
	write.table(probes.table, file=paste(paste(output.directory, project, sep="/"),"fdata.xls", sep="."), sep="\t",row.names=F, na="")
	
	cat("Writing phenoData table\n")
	pdata.table <- data.frame(sampleNames=sampleNames(eset), pData(eset))
	write.table(pdata.table, file=paste(paste(output.directory, project, sep="/"),"pdata.xls", sep="."), sep="\t",row.names=F, na="")
	
	for(i in 1:length(assayDataElementNames(eset))){
		cat("Writing",assayDataElementNames(eset)[i],"table\n")
		AD.table <- data.frame(featureNames=featureNames(eset), fdata.cols, assayData(eset)[[assayDataElementNames(eset)[i]]])
		write.table(AD.table, file=paste(paste(output.directory, project, sep="/"),".",assayDataElementNames(eset)[i],".xls", sep=""), sep="\t",row.names=F, na="")
			}
	
	if(is.element("GL", assayDataElementNames(eset))){
		cat("Writing GLAD table\n")
		gains <- apply(assayData(eset)$GL, 1, function(x) sum(x >= 1))
		losses <- apply(assayData(eset)$GL, 1, function(x) sum(x <= -1))
		amps <- apply(assayData(eset)$GL, 1, function(x) sum(x > 1))
		dels <- apply(assayData(eset)$GL, 1, function(x) sum(x < -1))
		GLAD.table <- data.frame(featureNames=featureNames(eset), fdata.cols, gains, losses, amps, dels)
		write.table(GLAD.table, file=paste(paste(output.directory, project, sep="/"),".GLAD.xls", sep=""), sep="\t",row.names=F, na="")
		}
	
		
	cat("Done\n")
	
}


removeProbesByMAD <- function(eset, mad.limit = NULL, n.probes)
{
	require(Biobase)
	featureData(eset)$MAD <- apply(exprs(eset),1,mad,na.rm=T)
	fvarMetadata(eset)$labelDescription[length(fvarMetadata(eset)$labelDescription)] <- "Median absolute deviation of expression"

	if(!is.null(mad.limit)){
		eset <- eset[featureData(eset)$MAD > mad.limit,]
		cat("",nrow(eset),"probes with a median absolute deviation greater than",mad.limit,"\n","Done\n") 
		}else{
		if(is.numeric(n.probes)){
			eset <- eset[rev(order(featureData(eset)$MAD)),]
			eset <- eset[1:n.probes,]
			min.mad <- min(featureData(eset)$MAD)
			cat(n.probes,"probes remaining\n","Done\n")
			}}

	eset
}


makeGTR <- function(gt,project)
#	gt is the result of hierarchical clustering of rows e.g.hclust(dist(ma$M), method="")
{
	height <- gt$height
 	height <- height + 1
  	height <- height[1]/height
	height <- signif(height, digits=5)
	node <- 1:length(height)
	node <- paste ('NODE',node,'X',sep='')
	merge1  <- gt$merge[,1]
	merge11 <- paste ('NODE',merge1,'X',sep='')
	merge12 <- paste ('GENE',-1-merge1,'X',sep='')
	merge1[gt$merge[,1]>0] <- merge11[gt$merge[,1]>0]
	merge1[gt$merge[,1]<0] <- merge12[gt$merge[,1]<0]
	merge2  <- gt$merge[,2]
	merge11 <- paste ('NODE',merge2,'X',sep='')
	merge12 <- paste ('GENE',-1-merge2,'X',sep='')
	merge2[gt$merge[,2]>0] <- merge11[gt$merge[,2]>0]
	merge2[gt$merge[,2]<0] <- merge12[gt$merge[,2]<0]
	data  <- data.frame(cbind(node,merge1,merge2))
	data  <- cbind(data,height)
	gtr.filename <- paste(project,"gtr", sep=".")
	write.table(data,file=gtr.filename,row.names=F, col.names=F,sep='\t')
}


makeATR <- function(at,project)
#	at is the result of hierarchical clustering of columns e.g.hclust(dist(t(ma$M)), method="")
{
  height <- at$height
  height <- height + 1
  height <- height[1]/height
  height  <- signif(height, digits = 5)
  node <- 1:length(height)
  node <- paste ('NODE',node,'X',sep='')
  merge1  <- at$merge[,1]
  merge11 <- paste ('NODE',merge1,'X',sep='')
  merge12 <- paste ('ARRY',-1-merge1,'X',sep='')
  merge1[at$merge[,1]>0] <- merge11[at$merge[,1]>0]
  merge1[at$merge[,1]<0] <- merge12[at$merge[,1]<0]
  merge2  <- at$merge[,2]
  merge11 <- paste ('NODE',merge2,'X',sep='')
  merge12 <- paste ('ARRY',-1-merge2,'X',sep='')
  merge2[at$merge[,2]>0] <- merge11[at$merge[,2]>0]
  merge2[at$merge[,2]<0] <- merge12[at$merge[,2]<0]
  data  <- data.frame(cbind(node,merge1,merge2))
  data  <- cbind(data,height)
  atr.filename <- paste(project,"atr", sep=".")
  write.table(data,file=atr.filename,row.names=F, col.names=F,sep='\t')

}


makeCDT <- function(gt, at, eset, Mtable=exprs(eset), project, pheno, sample.ids=NULL)
#	gt is the result of hierarchical clustering of rows e.g.hclust(dist(ma$M), method="")
#	at is the result of hierarchical clustering of columns e.g.hclust(dist(t(ma$M)), method="")
#	ma is the ma list object to be clustered
{
	require(Biobase)
	n <- nrow(Mtable)
	gweight <- matrix(1,ncol=1, nrow=nrow(Mtable))
	data <- data.frame(GID = paste("GENE", 0:(n - 1), "X", sep = ""),fData(eset), GWEIGHT=gweight, Mtable[, at$order], stringsAsFactors=F)
	data <- data[gt$order,  ]
	if(is.null(sample.ids)) sample.ids <- sampleNames(eset)
	names(data) <- c("GID",names(fData(eset)),"GWEIGHT", sample.ids[at$order])

	spacer <- rep(NA, length(names(fData(eset)))+1)
	nom <- c("EWEIGHT", spacer, rep("1", ncol(Mtable)))
	data <- rbind(nom, data)

	nom <- c("AID", spacer, paste("ARRY", (at$order - 1), "X", sep = ""))
	data <- rbind(nom, data)
	if(!is.null(pheno)){
		nom <- c("pheno", spacer, as.vector(pheno[at$order]))
		data <- rbind(nom, data)
		}
	
	cdt.filename <- paste(project,"cdt", sep=".")
	write.table(data, file = cdt.filename, row.names=F, col.names=T, sep = "\t", na="")
}


toTreeview <-function(eset, Mtable=exprs(eset), dist.method="correlation", clust.method="ward", cor.method="pearson", pheno=NULL, project)
#	hclust methods are
#	"ward", "single", "complete", "average", "mcquitty", "median", "centroid"
#	dist methods are
#	"euclidean", "maximum", "manhattan", "canberra", "binary", "minowski"
{

require(Biobase)
M <- Mtable
M[is.na(Mtable)] <- 0

if(dist.method == "correlation")
{
	atr <- hclust(dist(1-cor(M, method=cor.method)), method=clust.method)
	gtr <- hclust(dist(1-cor(t(M), method=cor.method)), method=clust.method)
	}else {
	atr <- hclust(dist(t(M), method=dist.method), method=clust.method)
	gtr <- hclust(dist((M), method=dist.method), method=clust.method)}

makeATR(at=atr, project=project)
makeGTR(gt=gtr, project=project)
makeCDT(eset=eset, gt=gtr, at=atr, Mtable=M, pheno=pheno, project=project)
atr$dist.method <- dist.method
atr$clust.method <- clust.method
atr
}


toTreeviewCGH <- function(cgh, Mtable=exprs(cgh), clust.method="ward", dist.method="euclidean", cor.method="pearson", chrom=NULL, project)
#	function to cluster aCGH samples and make spit ATR and CDT files for Treeview
{
	require(Biobase)
	cgh <- cgh[order(fData(cgh)$chrom, fData(cgh)$start),]
	M <- Mtable
	
	if(sum(is.na(Mtable)) != 0) stop("NAs in M table\nImpute values before clustering\n")
	if(!is.null(chrom)) {project <- paste(project,"chr",chrom, sep=".")}

	if(dist.method == "correlation") {atr <- hclust(dist(1-cor(M, method=cor.method)), method=clust.method)}else {
	atr <- hclust(dist(t(M), method=dist.method), method=clust.method)}

	makeATR(at=atr, project=project)

	M <- as.data.frame(M)
	names(M) <- sampleNames(cgh)
	M <- M[,atr$order]

	probes <- data.frame(featureNames=featureNames(cgh), fData(cgh), stringsAsFactors=F)
	head(probes)
	gweight <- as.data.frame(rep(1, nrow(M)))

	PM <- cbind(probes,gweight,M)
	if(!is.null(chrom)) {PM <- PM[which(fData(cgh)$chrom==chrom),]}
		
	lab.p <- c("AID",fvarLabels(cgh),"GWEIGHT")
	AIDs <- paste("ARRY", (atr$order - 1), "X", sep = "")
	labs.a <- c(lab.p, AIDs)

	s.lab.p <- c("Sample.id",rep("",ncol(probes)))
	s.lab.m <- sampleNames(cgh)[atr$order]
	s.labs.p <- c(s.lab.p, s.lab.m)
		
	e.lab.p <- c("EWEIGHT",rep("",ncol(probes)))
	e.lab.m <- (rep(1, ncol(M)))
	labs.e <- c(e.lab.p, e.lab.m)

#	header <- rbind.data.frame(labs.a,s.labs.p,labs.e)
#	names(header) <- names(PM)
	PGL <- rbind(labs.a,s.labs.p,labs.e,PM)

#	CDT <- as.data.frame(PGL)
	cdt.filename <- paste(project, "cdt", sep = ".")
	write.table(PGL, file = paste(project, "cdt", sep = "."), sep="\t", row.names=F, col.names=F)
}


repHybMean <- function(eset, rep.list)
{
	require(Biobase)
	eset.uni <- eset[,match(unique(rep.list), rep.list)]
	cat("Combining",ncol(eset.uni),"replicates\n")
	
	for(j in 1:length(assayDataElementNames(eset)))
	{
		AD.j <- assayData(eset)[[assayDataElementNames(eset)[j]]]
		AD.uni.j <- assayData(eset.uni)[[assayDataElementNames(eset.uni)[j]]]
			
			for (i in 1:ncol(eset.uni))
			{
			AD.i <- AD.j[,which(is.element(rep.list, unique(rep.list)[i])), drop=F]
			if(ncol(AD.i) > 1) AD.uni.j[,i] <- apply(AD.i, 1, mean)
			}
			
		assayDataElement(eset.uni, assayDataElementNames(eset.uni)[j]) <- AD.uni.j
	}
	
	sampleNames(eset.uni) <- unique(rep.list)
	
	if(!is.null(assayData(eset.uni)$flags)) {
		flags.uni <- assayData(eset.uni)$flags
		flags.uni[flags.uni != 1] <- 0
		assayDataElement(eset.uni, "flags") <- flags.uni
		}
		
	eset.uni	
	
}


dietSAMMulti <- function(eset, pheno, fData.ID="probeID", fData.genename="symbol", M.table=exprs(eset), q.value=5, local.FDR=NULL, logged2=TRUE, project, return.eset=T, nperms=100)
{
	require(Biobase)
	library(samr)
	set.seed(12345)
	
	M <- M.table[,which(!is.na(pheno))]
	eset <- eset[,which(!is.na(pheno))]
	pheno <- as.factor(pheno[which(!is.na(pheno))])
	
	geneid <-  fData(eset)[,match(fData.ID, names(fData(eset)))]
	genenames <- fData(eset)[,match(fData.genename, names(fData(eset)))]
	
	samr.data <- list(x=M,y=as.numeric(pheno), geneid=geneid,genenames=genenames, logged2=logged2)
	samr.result <-samr(samr.data,  resp.type="Multiclass", nperms=nperms)
	
	delta.table <- samr.compute.delta.table(samr.result)
	siggenes.table <- samr.compute.siggenes.table(samr.result, del=0, samr.data, delta.table, compute.localfdr=T, all.genes=T)
	
	sig.genes.sam <- siggenes.table$genes.up[which(as.numeric(siggenes.table$genes.up[,ncol(siggenes.table$genes.up)-1]) < q.value),,drop=F]

	if(!is.null(local.FDR)){
		sig.genes.sam <- sig.genes.sam[which(as.numeric(sig.genes.sam[,ncol(sig.genes.sam)])<local.FDR),,drop=F]
		}
	
	index <- as.numeric(sig.genes.sam[,1])-1
	
	sig.genes <- sig.genes.sam[,-1]
	
	if(nrow(sig.genes) <1) stop("No significant probes below FDR of ", local.FDR,"%\n")
	
	for (i in 1:length(levels(pheno))){
		
		exprs.i <- exprs(eset)[index,which(is.element(pheno, levels(pheno)[i]))]
		sig.genes <- cbind.data.frame(sig.genes, apply(exprs.i,1,mean))
		names(sig.genes)[ncol(sig.genes)] <- paste(levels(pheno)[i], "mean.expression", sep=".")
		
		sig.genes.sam <- cbind.data.frame(sig.genes.sam, apply(exprs.i,1,mean))
		names(sig.genes.sam)[ncol(sig.genes.sam)] <- paste(levels(pheno)[i], "mean.expression", sep=".")
		}
		
	sig.genes <- cbind.data.frame(sig.genes,fData(eset)[index,])
	
	names(sig.genes)[which(names(sig.genes)=="contrast-1"):
		(which(names(sig.genes)=="contrast-1")+length(levels(pheno))-1)] <- paste(levels(pheno), "contrast", sep=".")
	names(sig.genes.sam)[which(names(sig.genes.sam)=="contrast-1"):
		(which(names(sig.genes.sam)=="contrast-1")+length(levels(pheno))-1)] <- paste(levels(pheno), "contrast", sep=".")
		
	write.table(sig.genes, paste(project,"Multiclass.SAM.genes.xls",sep="."), sep="\t", row.names=F)
	
	cat(nrow(sig.genes), "significant probes with local FDR less than", local.FDR, "%\n")
	
	
	if(identical(return.eset,T)){
		eset.sig <- eset[index, ]
		fData(eset.sig) <- cbind.data.frame(fData(eset.sig), sig.genes.sam[,-1:-3])
		return(eset.sig)
	}
}


dietSAM2Class <- function(eset, pheno, project, fData.ID="probeID", fData.genename="symbol", resp.type="Two class unpaired", M.table=exprs(eset), nperms=100, q.value=5, local.FDR=NULL, fold.change=NULL, logged2=TRUE,  return.eset=T)
{
	
	require(Biobase)
	library(samr)
	set.seed(12345)
	
	M <- M.table[,which(!is.na(pheno))]
	eset <- eset[,which(!is.na(pheno))]
	pheno <- as.factor(pheno[which(!is.na(pheno))])
	geneid <-  fData(eset)[,match(fData.ID, names(fData(eset)))]
	genenames <- fData(eset)[,match(fData.genename, names(fData(eset)))]
	
	if(resp.type=="Two class paired") samr.data <- list(x=M,y=as.numeric(as.character(pheno)), geneid=geneid,genenames=genenames, logged2=logged2)
	if(resp.type=="Two class unpaired") samr.data <- list(x=M,y=as.numeric(pheno), geneid=geneid,genenames=genenames, logged2=logged2)

	samr.result <-samr(samr.data,  resp.type=resp.type, nperms=nperms)
	
	delta.table <- samr.compute.delta.table(samr.result)

	siggenes.table <- samr.compute.siggenes.table(samr.result, del=0, samr.data, delta.table, compute.localfdr=T, all.genes=T)
	
	response.2.genes <- siggenes.table$genes.up[as.numeric(siggenes.table$genes.up[,8]) < q.value,,drop=F]
	response.1.genes <- siggenes.table$genes.lo[as.numeric(siggenes.table$genes.lo[,8]) < q.value,,drop=F]
	
	response.2.genes[,7] <- as.numeric(response.2.genes[,7])
	response.1.genes[,7] <- 1/(as.numeric(response.1.genes[,7]))
	
	if(!is.null(local.FDR))
		{
		response.2.genes <- response.2.genes[as.numeric(response.2.genes[,9]) < local.FDR,,drop=F]
		response.1.genes <- response.1.genes[as.numeric(response.1.genes[,9]) < local.FDR,,drop=F]
		}
	
	if(!is.null(fold.change))
		{
		response.2.genes <- response.2.genes[response.2.genes[,7] > fold.change,,drop=F]
		response.1.genes <- response.1.genes[response.1.genes[,7] < fold.change,,drop=F]
		}
	
	response.2.genes.sam <- data.frame(response.2.genes, stringsAsFactors=F)
	names(response.2.genes.sam) <- c("Row", "Gene.symbol", "Gene.description", "score.d", "Numerator(r)", "Denominator(s+s0)", "Fold.change", "q.value", "local.fdr")
	genes.2.index <- as.numeric(response.2.genes.sam$Row)-1
	
	if(resp.type=="Two class paired"){
		response.2.exprs.1 <- exprs(eset)[genes.2.index,pheno<0]
		response.2.exprs.2 <- exprs(eset)[genes.2.index,pheno>0]
		}
	if(resp.type=="Two class unpaired"){
		if (length(genes.2.index)==1) {
			response.2.exprs.1 <- matrix(exprs(eset)[genes.2.index, which(is.element(pheno, levels(pheno)[1]))],nrow=1)
			response.2.exprs.2 <- matrix(exprs(eset)[genes.2.index,which(is.element(pheno, levels(pheno)[2]))],nrow=1)
		} else {
			response.2.exprs.1 <- exprs(eset)[genes.2.index,which(is.element(pheno, levels(pheno)[1]))]
			response.2.exprs.2 <- exprs(eset)[genes.2.index,which(is.element(pheno, levels(pheno)[2]))]
		}
	}

	expression.2 <- cbind.data.frame(apply(response.2.exprs.1,1,mean), apply(response.2.exprs.2,1,mean))
	names(expression.2) <- paste(levels(pheno), "expression", sep=".")
	response.2.genes.out <- cbind(response.2.genes.sam, expression.2)
	response.2.genes.out <- cbind(response.2.genes.out[,-1], fData(eset)[genes.2.index,])
	
	response.1.genes.sam <- data.frame(response.1.genes, stringsAsFactors=F)
	names(response.1.genes.sam) <- c("Row", "Gene.symbol", "Gene.description", "score.d", "Numerator(r)", "Denominator(s+s0)", "Fold.change", "q.value", "local.fdr")
	genes.1.index <- as.numeric(response.1.genes.sam$Row)-1
	
	if(resp.type=="Two class paired"){
		response.1.exprs.1 <- exprs(eset)[genes.1.index,pheno<0]
		response.1.exprs.2 <- exprs(eset)[genes.1.index,pheno>0]
		}
	if(resp.type=="Two class unpaired"){
		if (length(genes.1.index)==1) {
			response.1.exprs.1 <- matrix(exprs(eset)[genes.1.index,which(is.element(pheno, levels(pheno)[1]))],nrow=1)
			response.1.exprs.2 <- matrix(exprs(eset)[genes.1.index,which(is.element(pheno, levels(pheno)[2]))],nrow=1)
		} else { 	
			response.1.exprs.1 <- exprs(eset)[genes.1.index,which(is.element(pheno, levels(pheno)[1]))]
			response.1.exprs.2 <- exprs(eset)[genes.1.index,which(is.element(pheno, levels(pheno)[2]))]
		}
	}		
	
	expression.1 <- cbind.data.frame(apply(response.1.exprs.1,1,mean), apply(response.1.exprs.2,1,mean))
	names(expression.1) <- paste(levels(pheno), "expression", sep=".")
	response.1.genes.out <- cbind(response.1.genes.sam, expression.1)
	response.1.genes.out <- cbind(response.1.genes.out[,-1], fData(eset)[genes.1.index,])
	
	response.2.genes.out <- response.2.genes.out[rev(order(as.numeric(response.2.genes.out$Fold.change))),]
	response.1.genes.out <- response.1.genes.out[rev(order(as.numeric(response.1.genes.out$Fold.change))),]
	
	combined.response.genes.out <- rbind(response.2.genes.out, response.1.genes.out)
	combined.response.genes.sam <- rbind(response.2.genes.sam, response.1.genes.sam)
	
	if(resp.type=="Two class unpaired"){
		regulation <- c(rep(levels(pheno)[2], nrow(response.2.genes)), rep(levels(pheno)[1], nrow(response.1.genes)))
		combined.response.genes.out <- cbind(combined.response.genes.out, regulation)
		combined.response.genes.sam <- cbind(combined.response.genes.sam, regulation)
		
		if(nrow(combined.response.genes.out) <1) stop("No significant probes below Q value of ", q.value,"%\n")

		write.table(response.2.genes.out, paste(project,levels(pheno)[2],"significant.genes.xls",sep="."), sep="\t", row.names=F)
		write.table(response.1.genes.out, paste(project,levels(pheno)[1],"significant.genes.xls",sep="."), sep="\t", row.names=F)
		
		cat(nrow(response.2.genes.sam), levels(pheno)[2], "significant probes with Q value less than", q.value, "%\n")
		cat(nrow(response.1.genes.sam), levels(pheno)[1], "significant probes with Q value less than", q.value, "%\n")

	}
	
	if(resp.type=="Two class paired"){
		regulation <- c(rep("upregulated", nrow(response.2.genes)), rep("downregulated", nrow(response.1.genes)))
		combined.response.genes.out <- cbind(combined.response.genes.out, regulation)
		combined.response.genes.sam <- cbind(combined.response.genes.sam, regulation)
		
		if(nrow(combined.response.genes.out) <1) stop("No significant probes below Q value of ", q.value,"%\n")
		
		write.table(response.2.genes.out, paste(project,"significant.upregulated.genes.xls", sep="."), sep="\t", row.names=F)
		write.table(response.1.genes.out, paste(project,"significant.downregulated.genes.xls", sep="."), sep="\t", row.names=F)
		
		cat(nrow(response.2.genes.sam), "significant upregulated probes with Q value less than", q.value, "%\n")
		cat(nrow(response.1.genes.sam), "significant downregulated probes with Q value less than", q.value, "%\n")

	}
	
	write.table(combined.response.genes.out, paste(project,"combined.significant.genes.xls",sep="."), sep="\t", row.names=F)
	
	if(identical(return.eset,T)){
		index <- as.numeric(combined.response.genes.sam$Row)-1
		eset.sig <- eset[index, ]
		fData(eset.sig) <- cbind.data.frame(fData(eset.sig), combined.response.genes.sam[,-1:-3])
		head(fData(eset.sig))
		head(combined.response.genes.sam)
		return(eset.sig)
	}

}


dietLIMMA <- function(eset, pheno, return.eset=F, adj.pval.thresh=0.05, pval.thresh=NULL, project){
	
	require(Biobase)
	require(limma)
	eset <- eset[,!is.na(pheno)]
	pheno <- factor(pheno[!is.na(pheno)])
	
	limma.parameters <- cbind(intercept=rep(1, ncol(eset)), pheno)
	
	limma.fit <- lmFit(exprs(eset), limma.parameters)
	limma.fit <- eBayes(limma.fit)
	limma.table <- topTable(limma.fit, coef=2, number=nrow(eset), sort.by="none")
	
	means <- matrix(0, nrow=nrow(eset), ncol=length(levels(pheno)))
	for(i in 1:ncol(means)){
		means[,i] <- apply(exprs(eset)[,pheno==levels(pheno)[i]],1,mean)
		}
	means.df <- as.data.frame(means)
	names(means.df) <- as.character(levels(pheno))
	
	limma.table <- cbind.data.frame(limma.table, means.df, fData(eset))
	rownames(limma.table) <- featureNames(eset)
	
	up.significant <- limma.table[which(limma.table$logFC > 0),]
	up.significant <- up.significant[rev(order(up.significant$logFC)),]
	if(!is.null(pval.thresh)) {
		up.significant <- up.significant[which(up.significant$P.Value <pval.thresh),]
		up.significant <- up.significant[order(up.significant$P.Value),]
		}
	if(!is.null(adj.pval.thresh)){
		up.significant <- up.significant[which(up.significant$adj.P.Val <adj.pval.thresh),]
		up.significant <- up.significant[order(up.significant$adj.P.Val),]
		}

	down.significant <- limma.table[which(limma.table$logFC < 0),]
	down.significant <- down.significant[rev(order(down.significant$logFC)),]
	if(!is.null(pval.thresh)) {
		down.significant <- down.significant[which(down.significant$P.Value <pval.thresh),]
		down.significant <- down.significant[order(down.significant$P.Value),]
		}
	if(!is.null(adj.pval.thresh)){
		down.significant <- down.significant[which(down.significant$adj.P.Val <adj.pval.thresh),]
		down.significant <- down.significant[order(down.significant$adj.P.Val),]
		}
	
	
	write.table(up.significant, file=paste(project, levels(pheno)[2], "significant.genes.xls", sep="."), sep="\t", row.names=F, na="")
	write.table(down.significant, file=paste(project, levels(pheno)[1], "significant.genes.xls", sep="."), sep="\t", row.names=F, na="")
	
	cat(nrow(up.significant), levels(pheno)[2], "significant genes\n")
	cat(nrow(down.significant), levels(pheno)[1], "significant genes\n")
	
	eset.sig <- eset[featureNames(eset) %in% c(rownames(up.significant), rownames(down.significant)),]
	if(identical(return.eset,T)) return(eset.sig)

}

dietTwilight <- function(eset, pheno, pval.thresh=0.05, fdr=NULL, perms=100, project, return.eset=F){
	require(twilight)
	eset <- eset[,!is.na(pheno)]
	pheno <- pheno[!is.na(pheno)]
	as.numeric(factor(pheno))
	
	twilight.pval <- twilight.pval(exprs(eset), as.numeric(factor(pheno)), B=perms)
	
	eset.twig <- twilight(twilight.pval)
	
	eset.twig.result <- eset.twig$result[,-which(is.element(names(eset.twig$result), c("candidate", "mean.fdr", "upper.fdr", "lower.fdr")))]
	
	eset.twig.result <- data.frame(featureNames=featureNames(eset), eset.twig.result[featureNames(eset),], fData(eset))
	
	#	expression means
	mean.exprs.1 <- apply(exprs(eset)[,which(is.element(pheno, levels(pheno)[1]))],1,mean)
	mean.exprs.2 <- apply(exprs(eset)[,which(is.element(pheno, levels(pheno)[2]))],1,mean)
	
	eset.twig.result.out <- cbind(eset.twig.result, mean.exprs.1, mean.exprs.2)
	
	names(eset.twig.result.out)[(ncol(eset.twig.result.out)-1):ncol(eset.twig.result.out)] <- paste(levels(pheno), "mean.expression", sep=".")
	
	eset.twig.result.out <- eset.twig.result.out[order(eset.twig.result.out$pvalue),]
	
	write.table(eset.twig.result.out, paste(project,"twilight.all.probes.xls",sep="."), sep="\t", row.names=F, na="")
	
	if(!is.null(pval.thresh)){
		
		eset.twig.sig <- eset.twig.result.out[which(eset.twig.result.out$pvalue < pval.thresh),]
		if(!is.null(fdr)){
			eset.twig.sig <- eset.twig.sig[which(eset.twig.sig$fdr < fdr)]
			}
		
		eset.twig.sig1 <- eset.twig.sig[eset.twig.sig[,(ncol(eset.twig.sig)-1)] > eset.twig.sig[,(ncol(eset.twig.sig))],]
		eset.twig.sig2 <- eset.twig.sig[eset.twig.sig[,(ncol(eset.twig.sig)-1)] < eset.twig.sig[,(ncol(eset.twig.sig))],]
		
		cat(nrow(eset.twig.sig1), levels(pheno)[1], "significant probes with p value less than", pval.thresh, "%\n")
		cat(nrow(eset.twig.sig2), levels(pheno)[2], "significant probes with p value less than", pval.thresh, "%\n")
		write.table(eset.twig.sig1, paste(project,levels(pheno)[1],"Twilight.significant.genes.xls",sep="."), sep="\t", row.names=F, na="")
		write.table(eset.twig.sig2, paste(project,levels(pheno)[2],"Twilight.significant.genes.xls",sep="."), sep="\t", row.names=F, na="")
		}

	
}



calculateCentroids <- function(eset, pheno, project, return.eset=F)
{
	require(Biobase)
	pheno <- pheno[!is.na(pheno)]
	
	eset <- eset[, which(!is.na(pheno))]
	

	phenotypes <- as.factor(pheno)
	
	centroids <- matrix(0,nrow(eset), length(levels(phenotypes)))
	phenos <- levels(phenotypes)
	
	for (i in 1:length(phenos)){
		M.i <- exprs(eset)[,phenotypes == phenos[i]]
		centroids[,i] <- apply(M.i,1,mean, na.rm=T)/apply(M.i,1,sd, na.rm=T)
		}
	
	centroids <- as.data.frame(centroids)
	names(centroids) <- as.character(phenos)
	centroids.out <- cbind(centroids,fData(eset))
	
	write.table(centroids.out, file=paste(project,"centroids.txt", sep="."), sep="\t", row.names=F, na="", quote=F)
	
	if(identical(return.eset,T)){
		names(centroids) <- paste(as.character(phenos),"centroid", sep=".")
		fData(eset) <- cbind(fData(eset),centroids)
		return(eset)
	}
	
}


weightedAverageSignature <- function (eset, eset.pheno, project, gene.list.file, gene.list.metric, main="Weighted average signature", eset.id="symbol", device="PDF", pheno.colours=NULL, return.eset=F, eset.pheno.label)
{
	require(RColorBrewer)
	if(!is.null(eset.pheno)){
		eset <- eset[,!is.na(eset.pheno)]
		eset.pheno <- factor(eset.pheno[!is.na(eset.pheno)])
		eset <- eset[,order(eset.pheno)]
		eset.pheno <- eset.pheno[order(eset.pheno)]
		}
	
	#	read in the gene list
	gene.list <- read.delim(gene.list.file, header=T, sep="\t", stringsAsFactors=F, na.strings=c("", " ","NA", "N/A", "#N/A" ))
	gene.list <- gene.list[!is.na(gene.list[,eset.id]),]
	gene.list <- gene.list[order(gene.list[,gene.list.metric], decreasing=T),]
	#	remove duplicates
	gene.list <- gene.list[!duplicated(gene.list[,eset.id]),]
	
	#	match gene list and eset
	eset.uni <- removeReplicateProbesByVariance(eset, eset.id=eset.id, verbose=F)
	eset.match <- eset.uni[fData(eset.uni)[,eset.id] %in% gene.list[,eset.id],]
	gene.list.match <- gene.list[gene.list[,eset.id] %in% fData(eset.match)[,eset.id] ,]
	eset.match <- eset.match[match(gene.list.match[,eset.id], fData(eset.match)[,eset.id]),]

	cat(sum(fData(eset.match)[,eset.id] == gene.list.match[,eset.id]), "matching identifiers in Gene List and Expression Set\n")
	
	#	calculate weighted averages
	weight.matrix <- matrix(gene.list.match[,gene.list.metric], nrow(eset.match), ncol(eset.match))
	weighted.average <- exprs(eset.match)*weight.matrix
	signat.scores <- apply(weighted.average, 2, mean)
	signat.scores.mad <- apply(weighted.average, 2, mad)
	signat.scores.mad[signat.scores < 0] <- 0-signat.scores.mad[signat.scores < 0]

	
	#	write out results table
	scores.table <- data.frame(sampleNames=sampleNames(eset.match), Signature.scores=signat.scores, pData(eset.match))
	scores.table <- scores.table[order(signat.scores, decreasing=T),]
	write.table(scores.table, file=paste(project, ".signature.scores.xls", sep=""), sep="\t",row.names=F)
	weighted.average.table <- data.frame(featureNames=featureNames(eset.match), weighted.average, fData(eset.match), gene.list.match)
	write.table(weighted.average.table, file=paste(project, ".weighted.average.expression.xls", sep=""), sep="\t",row.names=F)
	expression.table <- data.frame(featureNames=featureNames(eset.match), exprs(eset.match), fData(eset.match), gene.list.match)
	write.table(weighted.average.table, file=paste(project, ".expression.xls", sep=""), sep="\t",row.names=F)
	
	if(!is.element("barplots", list.files())) dir.create("barplots")
	if(!is.element("boxplots", list.files())) dir.create("boxplots")
	
	
	#	barplot waterfall coloured for eset pheno
	if(length(device) !=0){
	
	
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("barplots/",project,".barplots.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("barplots/",project,".barplots.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("barplots/",project,".barplots.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("barplots/",project,".barplots.ps",sep=""),width=9, height=6, pointsize=12)}}
		
	signat.scores.sort <-sort(signat.scores, index.return = TRUE)
	if(is.null(pheno.colours)) pheno.colours <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
	pheno.bar.colours <- pheno.colours[as.numeric(eset.pheno)]
	
	barplot2(signat.scores[signat.scores.sort$ix],beside=TRUE, main=main, names.arg=sampleNames(eset.match)[signat.scores.sort$ix],las=2, cex.names=0.5,col=pheno.bar.colours[signat.scores.sort$ix], plot.ci=T, ci.u=signat.scores[signat.scores.sort$ix]+signat.scores.mad[signat.scores.sort$ix], ci.l=signat.scores[signat.scores.sort$ix]-signat.scores.mad[signat.scores.sort$ix])
	
	legend("topleft",legend=unique(eset.pheno[signat.scores.sort$ix]),fill=unique(pheno.bar.colours[signat.scores.sort$ix]),cex=0.8, bty="n")
	
	
	#	boxplot of average expression
	if(length(device) != 0 && device != "quartz") dev.off()
	if(length(device) !=0){
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("boxplots/",project,".expression.up.boxplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("boxplots/",project,".expression.up.boxplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("boxplots/",project,".expression.up.boxplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("boxplots/",project,".expression.up.boxplot.ps",sep=""),width=9, height=6, pointsize=12)}}
		
	boxplot(exprs(eset.match)[gene.list.match[,gene.list.metric] > 0,], range=1, outline=F, col=pheno.bar.colours, main=paste(main, "upregulated expression"))
	for(i in 1:length(eset.pheno)) stripchart(exprs(eset.match)[gene.list.match[,gene.list.metric] > 0,i], col="black", add=T,method="jitter", jitter=0.2, vertical=TRUE, pch="+", at=i, cex=0.25)
	legend("topleft",legend=unique(eset.pheno),fill=unique(pheno.bar.colours),cex=0.8, bty="n")
	
	if(length(device) != 0 && device != "quartz") dev.off()
	

	if(length(device) !=0){
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("boxplots/",project,".expression.down.boxplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("boxplots/",project,".expression.down.boxplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("boxplots/",project,".expression.down.boxplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("boxplots/",project,".expression.down.boxplot.ps",sep=""),width=9, height=6, pointsize=12)}}
		
	boxplot(exprs(eset.match)[gene.list.match[,gene.list.metric] < 0,], range=1, outline=F, col=pheno.bar.colours, main=paste(main, "downregulated expression"))
	for(i in 1:length(eset.pheno)) stripchart(exprs(eset.match)[gene.list.match[,gene.list.metric] < 0,i], col="black", add=T,method="jitter", jitter=0.2, vertical=TRUE, pch="+", at=i, cex=0.25)
	legend("topleft",legend=unique(eset.pheno),fill=unique(pheno.bar.colours),cex=0.8, bty="n")
	
	if(length(device) != 0 && device != "quartz") dev.off()
	
	
	#	boxplot of weighted expression
	if(length(device) !=0){
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("boxplots/",project,".weighted.boxplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("boxplots/",project,".weighted.boxplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("boxplots/",project,".weighted.boxplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("boxplots/",project,".weighted.boxplot.ps",sep=""),width=9, height=6, pointsize=12)}}
	
	boxplot(weighted.average, range=1, outline=F, col=pheno.bar.colours, main=paste(main, "combined weighted expression"))
	for(i in 1:length(eset.pheno)) stripchart(weighted.average[,i], col="black", add=T,method="jitter", jitter=0.2, vertical=TRUE, pch="+", at=i, cex=0.5)
	legend("topleft",legend=unique(eset.pheno),fill=unique(pheno.bar.colours),cex=0.8, bty="n")
	
	if(length(device) != 0 && device != "quartz") dev.off()
	
	#	boxplot of weighted up expression
	if(length(device) !=0){
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("boxplots/",project,".weighted.up.boxplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("boxplots/",project,".weighted.up.boxplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("boxplots/",project,".weighted.up.boxplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("boxplots/",project,".weighted.up.boxplot.ps",sep=""),width=9, height=6, pointsize=12)}}
		
	weighted.average.up <- weighted.average[gene.list.match[,gene.list.metric] > 0,]
	boxplot(weighted.average.up, range=1, outline=F, col=pheno.bar.colours, main=paste(main, "upregulated weighted expression"))
	for(i in 1:length(eset.pheno)) stripchart(weighted.average.up[,i], col="black", add=T,method="jitter", jitter=0.2, vertical=TRUE, pch="+", at=i, cex=0.5)
	if(length(device) != 0 && device != "quartz") dev.off()
	
	#	boxplot of weighted up expression
	if(length(device) !=0){
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("boxplots/",project,".weighted.down.boxplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("boxplots/",project,".weighted.down.boxplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("boxplots/",project,".weighted.down.boxplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("boxplots/",project,".weighted.down.boxplot.ps",sep=""),width=9, height=6, pointsize=12)}}
	
	weighted.average.down <- 0-weighted.average[gene.list.match[,gene.list.metric] < 0,]
	boxplot(weighted.average.down, range=1, outline=F, col=pheno.bar.colours, main=paste(main, "downregulated weighted expression"))
	for(i in 1:length(eset.pheno)) stripchart(weighted.average.down[,i], col="black", add=T,method="jitter", jitter=0.2, vertical=TRUE, pch="+", at=i, cex=0.5)
	if(length(device) != 0 && device != "quartz") dev.off()
	
	if(is.null(pData(eset.match))) {
			p.data <- data.frame(signat.scores)
			row.names(p.data) <- sampleNames(eset.match)
			pData(eset.match) <- p.data
			}else{
			pData(eset.match) <- cbind(pData(eset.match),signif(signat.scores,3))
			}
	if(return.eset==T) return(eset.match)
	
}

expressionSignature <- function (eset, eset.pheno, project, gene.list.file, gene.list.metric, main=NULL, eset.id="symbol", device="PDF", pheno.colours=NULL, return.eset=F, cex.axis=0.5, plot.jitter=F, cex.jitter=0.4)
{
	require(RColorBrewer)
	require(gplots)
	require(genefilter)
	if(!is.null(eset.pheno)){
		eset <- eset[,!is.na(eset.pheno)]
		eset.pheno <- factor(eset.pheno[!is.na(eset.pheno)])
		eset <- eset[,order(eset.pheno)]
		eset.pheno <- eset.pheno[order(eset.pheno)]
		}
	
	#	read in the gene list
	gene.list <- read.delim(gene.list.file, header=T, sep="\t", stringsAsFactors=F, na.strings=c("", " ","NA", "N/A", "#N/A" ))
	gene.list <- gene.list[!is.na(gene.list[,eset.id]),]
	gene.list <- gene.list[order(gene.list[,gene.list.metric], decreasing=T),]
	#	remove duplicates
	gene.list <- gene.list[!duplicated(gene.list[,eset.id]),]
	
	#	match gene list and eset
	eset.uni <- removeReplicateProbesByVariance(eset, eset.id=eset.id, verbose=F)
	eset.match <- eset.uni[fData(eset.uni)[,eset.id] %in% gene.list[,eset.id],]
	gene.list.match <- gene.list[gene.list[,eset.id] %in% fData(eset.match)[,eset.id] ,]
	eset.match <- eset.match[match(gene.list.match[,eset.id], fData(eset.match)[,eset.id]),]

	cat(sum(fData(eset.match)[,eset.id] == gene.list.match[,eset.id]), "matching identifiers in Gene List and Expression Set\n")
	
	eset.match.up <- eset.match[gene.list.match[,gene.list.metric] > 0,]
	eset.match.down <- eset.match[gene.list.match[,gene.list.metric] <= 0,]
	
	expression.scores.up <- apply(exprs(eset.match.up), 2, mean)
	eset.match.up <- eset.match.up[,order(expression.scores.up)]
	eset.pheno.up <- eset.pheno[order(expression.scores.up)]
	expression.scores.up <- expression.scores.up[order(expression.scores.up)]
	expression.scores.up.mad <- apply(exprs(eset.match.up), 2, mad)
	expression.scores.up.mad[expression.scores.up < 0] <- 0-expression.scores.up.mad[expression.scores.up < 0]
	
	upperMAD.up <- lowerMAD.up <- expression.scores.up.mad
	upperMAD.up[expression.scores.up < 0] <- 0
	lowerMAD.up[expression.scores.up >= 0] <- 0
	upperCI.up <- expression.scores.up+upperMAD.up
	lowerCI.up <- expression.scores.up+lowerMAD.up
	
	expression.scores.down <- apply(exprs(eset.match.down), 2, mean)
	eset.match.down <- eset.match.down[,order(expression.scores.down)]
	eset.pheno.down <- eset.pheno[order(expression.scores.down)]
	expression.scores.down <- expression.scores.down[order(expression.scores.down)]
	expression.scores.down.mad <- apply(exprs(eset.match.down), 2, mad)
	expression.scores.down.mad[expression.scores.down < 0] <- 0-expression.scores.down.mad[expression.scores.down < 0]
	
	upperMAD.down <- lowerMAD.down <- expression.scores.down.mad
	upperMAD.down[expression.scores.down < 0] <- 0
	lowerMAD.down[expression.scores.down >= 0] <- 0
	upperCI.down <- expression.scores.down+upperMAD.down
	lowerCI.down <- expression.scores.down+lowerMAD.down
	
	if(is.null(pheno.colours)) pheno.colours <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
	pheno.bar.colours <- pheno.colours[as.numeric(eset.pheno)]
	pheno.bar.colours.up <- pheno.colours[as.numeric(eset.pheno.up)]
	pheno.bar.colours.down <- pheno.colours[as.numeric(eset.pheno.down)]
	
	t.up.matrix <- matrix(NA, ncol(eset.match.up), ncol(eset.match.up))
	for(i in 1:ncol(eset.match.up)){
		for(j in 1:ncol(eset.match.up)){
			t.up.matrix[i,j] <- signif(t.test(exprs(eset.match.up)[,i], exprs(eset.match.up)[,j])$p.value, 5)
			}
		}
	colnames(t.up.matrix) <- sampleNames(eset.match.up)
	t.up.matrix <- t.up.matrix[match(sampleNames(eset.match), sampleNames(eset.match.up)),match(sampleNames(eset.match), sampleNames(eset.match.up))]
	colnames(t.up.matrix) <- paste("upregulated.ttest.p.vs", colnames(t.up.matrix), sep=".")
	
	t.down.matrix <- matrix(NA, ncol(eset.match.down), ncol(eset.match.down))
	for(i in 1:ncol(eset.match.down)){
		for(j in 1:ncol(eset.match.down)){
			t.down.matrix[i,j] <- signif(t.test(exprs(eset.match.down)[,i], exprs(eset.match.down)[,j])$p.value, 5)
			}
		}
	colnames(t.down.matrix) <- sampleNames(eset.match.down)
	t.down.matrix <- t.down.matrix[match(sampleNames(eset.match), sampleNames(eset.match.down)),match(sampleNames(eset.match), sampleNames(eset.match.down))]
	colnames(t.down.matrix) <- paste("downregulated.ttest.p.vs", colnames(t.down.matrix), sep=".")
	
	
	if(!is.element("signature.barplots", list.files())) dir.create("signature.barplots")
	if(!is.element("signature.boxplots", list.files())) dir.create("signature.boxplots")
	if(!is.element("signature.tables", list.files())) dir.create("signature.tables")
	
	expression.scores.up[match(sampleNames(eset.match), sampleNames(eset.match.up))]
	
	#	write out results table
	scores.table <- data.frame(sampleNames=sampleNames(eset.match), Expression.scores.up=expression.scores.up, Expression.scores.up.MAD=abs(expression.scores.up.mad), t.up.matrix, Expression.scores.down=expression.scores.down, Expression.scores.down.MAD=abs(expression.scores.down.mad), t.down.matrix, pData(eset.match))
	write.table(scores.table, file=paste("signature.tables/",paste(project, ".expression.scores.xls", sep=""), sep=""), sep="\t",row.names=F, na="")
	
	expression.table <- data.frame(featureNames=featureNames(eset.match), exprs(eset.match), fData(eset.match), gene.list.match)
	expression.up.table <- data.frame(featureNames=featureNames(eset.match.up), exprs(eset.match.up), fData(eset.match.up), gene.list.match[gene.list.match[,gene.list.metric] > 0,])
	expression.down.table <- data.frame(featureNames=featureNames(eset.match.down), exprs(eset.match.down), fData(eset.match.down), gene.list.match[gene.list.match[,gene.list.metric] <= 0,])
	
	write.table(expression.table, file=paste("signature.tables/",paste(project, ".expression.xls", sep="") ,sep=""), sep="\t",row.names=F, na="")
	write.table(expression.up.table, file=paste("signature.tables/",paste(project, ".upregulated.expression.xls", sep="") ,sep=""), sep="\t",row.names=F, na="")
	write.table(expression.down.table, file=paste("signature.tables/",paste(project, ".downregulated.expression.xls", sep="") ,sep=""), sep="\t",row.names=F, na="")
	
	
	if(is.null(main)) main <- project
	#	barplot waterfalls coloured for eset pheno
	if(length(device) !=0){
	
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("signature.barplots/",project,".expression.up.barplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("signature.barplots/",project,".expression.up.barplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("signature.barplots/",project,".expression.up.barplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("signature.barplots/",project,".expression.up.barplot.ps",sep=""),width=9, height=6, pointsize=12)}}
	
	barplot2(expression.scores.up,beside=TRUE, main=paste(main, "\nupregulated expression"), names.arg=sampleNames(eset.match.up),las=2, cex.names=0.5,col=pheno.bar.colours.up, plot.ci=T, ci.u=upperCI.up, ci.l=lowerCI.up)

legend("topleft",legend=unique(eset.pheno),fill=unique(pheno.bar.colours),cex=0.8, bty="n")

	if(length(device) != 0 && device != "quartz") dev.off()
	
	if(length(device) !=0){
	
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("signature.barplots/",project,".expression.down.barplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("signature.barplots/",project,".expression.down.barplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("signature.barplots/",project,".expression.down.barplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("signature.barplots/",project,".expression.down.barplot.ps",sep=""),width=9, height=6, pointsize=12)}}
		
	barplot2(expression.scores.down,beside=TRUE, main=paste(main, "\ndownregulated expression"), names.arg=sampleNames(eset.match.down),las=2, cex.names=0.5,col=pheno.bar.colours.down, plot.ci=T, ci.u=upperCI.down, ci.l=lowerCI.down)

legend("topleft",legend=unique(eset.pheno),fill=unique(pheno.bar.colours),cex=0.8, bty="n")

if(length(device) != 0 && device != "quartz") dev.off()

	
	
	#	boxplot of average expression
	
	if(length(device) !=0){
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("signature.boxplots/",project,".expression.up.boxplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("signature.boxplots/",project,".expression.up.boxplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("signature.boxplots/",project,".expression.up.boxplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("signature.boxplots/",project,".expression.up.boxplot.ps",sep=""),width=9, height=6, pointsize=12)}}
		
	boxplot(exprs(eset.match.up), range=1, outline=F, col=pheno.bar.colours.up, main=paste(main, "\nupregulated expression"), las=2, cex.axis=cex.axis)
	if(plot.jitter) for(i in 1:length(eset.pheno.up)) stripchart(exprs(eset.match.up)[,i], col="black", add=T,method="jitter", jitter=0.2, vertical=TRUE, pch="+", at=i, cex=cex.jitter)
	legend("topleft",legend=unique(eset.pheno.up),fill=unique(pheno.bar.colours.up),cex=0.8, bty="n")
	
	if(length(device) != 0 && device != "quartz") dev.off()
	

	if(length(device) !=0){
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("signature.boxplots/",project,".expression.down.boxplot.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("signature.boxplots/",project,".expression.down.boxplot.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("signature.boxplots/",project,".expression.down.boxplot.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("signature.boxplots/",project,".expression.down.boxplot.ps",sep=""),width=9, height=6, pointsize=12)}}
		
	boxplot(exprs(eset.match.down), range=1, outline=F, col=pheno.bar.colours.down, main=paste(main, "\ndownregulated expression"), las=2, cex.axis=cex.axis)
	if(plot.jitter) for(i in 1:length(eset.pheno.down)) stripchart(exprs(eset.match.down)[,i], col="black", add=T,method="jitter", jitter=0.2, vertical=TRUE, pch="+", at=i, cex=cex.jitter)
	legend("topleft",legend=unique(eset.pheno.down),fill=unique(pheno.bar.colours.down),cex=0.8, bty="n")
	
	if(length(device) != 0 && device != "quartz") dev.off()
	
			
	if(return.eset==T) return(eset.match)
	
}

centroidCorrelation <- function(eset, centroid.file, centroids, centroid.id="symbol", eset.id="symbol", project, cor.method="spearman", return.eset=T, pheno.label=NULL, return.exprs=F)
{
	require(Biobase)
	require(Biobase)
	if(!is.element("centroid.correlations", list.files())) dir.create("centroid.correlations")
	
	centroid <- read.delim(centroid.file, header=T, sep="\t", stringsAsFactors=F, na.strings=c("", " ","NA", "N/A", "#N/A" ))
	
	centroid.index <- which(names(centroid) == centroid.id)
	eset.index <- which(fvarLabels(eset) == eset.id)
	
	centroid.na <- centroid[!is.na(centroid[,centroid.index]),]
	eset.na <- eset[!is.na(fData(eset)[,eset.index]),]
	
	eset.probes <- fData(eset.na)[,eset.index]
	centroid.probes <- centroid.na[,centroid.index]
	
	eset.sub <- eset.na[eset.probes %in% centroid.probes,]
	
	eset.variance <- apply(exprs(eset.sub),1,var,na.rm=T)
	eset.sorted <- eset.sub[rev(order(eset.variance)),]
	eset.sorted.ids <- fData(eset.sorted)[,eset.index]
	eset.sub.uni <- eset.sorted[!duplicated(eset.sorted.ids),]
	
	centroid.sub <- centroid.na[centroid.probes %in% eset.probes,]
	centroid.variance <- apply(centroid.sub[,1:centroids],1,var,na.rm=T)
	centroid.sub.sorted <- centroid.sub[rev(order(centroid.variance)),]
	centroid.sub.uni <- centroid.sub.sorted[match(unique(centroid.sub.sorted[,centroid.index]),centroid.sub.sorted[,centroid.index]),]
	
	centroid.sub.uni <- centroid.sub.uni[order(centroid.sub.uni[,centroid.index]),]
	eset.sub.uni <- eset.sub.uni[order(fData(eset.sub.uni)[,eset.index]),]
	
	cat(nrow(eset.sub.uni), "matching identifiers\n")
	
	centroid.values <- centroid.sub.uni[,1:centroids]
	names(centroid.values) <- paste(names(centroid.values),"centroid", sep=".")
	
	if(return.exprs == T){
	centroid.out <- cbind(centroid.values, centroid.sub.uni[,centroid.index], fData(eset.sub.uni), exprs(eset.sub.uni))
	write.table(centroid.out, file=paste("centroid.correlations/", project, ".centroid.M.table.xls", sep=""), sep="\t",row.names=F, na="")
	}
	
	centroid.correlations <- matrix(NA, ncol(exprs(eset)), centroids+3)
	for(i in 1:centroids)
	{
		for(j in 1:ncol(eset)){
			centroid.correlations[j,i] <- signif(cor(centroid.sub.uni[,i], exprs(eset.sub.uni)[,j], method=cor.method, use="complete.obs"),3)
			}
	}
	
	centroid.correlations <- as.data.frame(centroid.correlations)
	
	for(i in 1:nrow(centroid.correlations))
		{
		centroid.correlations[i,(centroids+1)] <- names(centroid.sub.uni[,(1:centroids)])[which(centroid.correlations[i,(1:centroids)] == max(centroid.correlations[i,1:centroids]))][1]
		centroid.correlations[i,(centroids+2)] <- centroid.correlations[i,which(centroid.correlations[i,(1:centroids)] == max(centroid.correlations[i,(1:centroids)]))][1]
		}
	
	centroid.correlations[centroid.correlations[,centroids+2]>0.1,centroids+3] <- centroid.correlations[centroid.correlations[,centroids+2]>0.1,centroids+1]
	names(centroid.correlations)[1:centroids]  <- paste(names(centroid.sub)[1:centroids],"correlation", sep=".")
	rownames(centroid.correlations) <- sampleNames(eset)
	if(is.null(pheno.label)) pheno.label <- project
	names(centroid.correlations)[centroids+1] <- paste(pheno.label,"Nearest.centroid", sep=".")
	names(centroid.correlations)[centroids+2] <- paste(pheno.label,paste(cor.method,"cor", sep="."), sep=".")
	names(centroid.correlations)[centroids+3] <- paste(pheno.label,"Nearest.centroid.over.0.1", sep=".")
	
	centroid.correlations.out <- data.frame(sample.id=sampleNames(eset.sub.uni), centroid.correlations)
	centroid.correlations.out[is.na(centroid.correlations.out[,ncol(centroid.correlations.out)]),ncol(centroid.correlations.out)] <- NA
	
	centroid.correlations.out <- cbind(centroid.correlations.out, pData(eset.sub.uni))
	rownames(centroid.correlations.out) <- sampleNames(eset)
	write.table(centroid.correlations.out, file=paste("centroid.correlations/",project, ".nearest.centroids.xls", sep=""), sep="\t",row.names=F)
	
	if(is.null(pData(eset))) {
			p.data <- data.frame(sampleNames=sampleNames(est), centroid.correlations)
			row.names(p.data) <- sampleNames(eset)
			pData(eset) <- p.data
			}else{
			pData(eset) <- cbind(pData(eset),centroid.correlations[,(ncol(centroid.correlations)-2):ncol(centroid.correlations)])
			}
			validObject(eset)
			rm(centroid.correlations)
	
	if(return.eset) return(eset)
		
}

removeReplicateProbesByVariance <- function(eset, eset.id="symbol", verbose=T)
{
#	commands <- c(eset$commands,deparse(match.call()))
	require(Biobase)
	probe.col <- which(names(fData(eset)) == eset.id)
	eset <- eset[!is.na(fData(eset)[,probe.col]),]
	eset.variance <- apply(exprs(eset),1,var,na.rm=T)
	
	eset.sorted <- eset[rev(order(eset.variance)),]
	sorted.ids <- fData(eset.sorted)[,probe.col]
	eset.out <- eset.sorted[!duplicated(sorted.ids),]
	
	if(identical(verbose, T)) cat("",nrow(eset.out),"probes remaining\n","Done\n")
	eset.out
}

removeReplicateProbesByMAD <- function(eset, eset.id="symbol", verbose=T)
{
#	commands <- c(eset$commands,deparse(match.call()))
	require(Biobase)
	probe.col <- which(names(fData(eset)) == eset.id)
	eset <- eset[!is.na(fData(eset)[,probe.col]),]
	
	eset.mad <- apply(exprs(eset),1,mad,na.rm=T)
	
	eset.sorted <- eset[rev(order(eset.mad)),]
	sorted.ids <- fData(eset.sorted)[,probe.col]
	eset.out <- eset.sorted[!duplicated(sorted.ids),]
	
	if(identical(verbose, T)) cat("",nrow(eset.out),"probes remaining\n","Done\n")
	eset.out
}

removeReplicateProbesByCorrelation <- function(eset, eset.cor="pearson.adjp", eset.id="symbol", verbose=T)
{
#	commands <- c(ma$commands,deparse(match.call()))
	require(Biobase)
	probe.col <- which(names(fData(eset)) == eset.id)
	eset <- eset[!is.na(fData(eset)[,probe.col]),]
	cor.col <- which(names(fData(eset)) == eset.cor)
	
	eset.sorted <- eset[rev(order(fData(eset)[,cor.col])),]
	sorted.ids <- fData(eset.sorted)[,probe.col]
	eset.out <- eset.sorted[!duplicated(sorted.ids),]
	
	if(identical(verbose, T)) cat("",nrow(eset.out),"probes remaining\n","Done\n")
	eset.out
}

removeReplicateProbesByRange <- function (eset, eset.id = "symbol", verbose = T) 
{
     require(Biobase)
     probe.col <- which(names(fData(eset)) == eset.id)
     eset <- eset[!is.na(fData(eset)[, probe.col]), ]
     drange <- unlist(apply(exprs(eset),1,function(x) {max(x,na.rm=T)-min(x,na.rm=T)}))
     ss <- tapply(drange, as.factor(fData(eset)[, probe.col]), function(x) {max(x)})
     dd <- tapply(1:length(fData(eset)[, probe.col]), as.factor(fData(eset)[, probe.col]), function(x) {x[which.max(drange[x])]})
     eset.out <- eset[dd, ]
     if (identical(verbose, T)) 
         cat("", nrow(eset.out), "probes remaining\n", "Done\n")
     eset.out
}



multiHeatmap <- function(eset, cluster=exprs(eset), chrom.order=F, heatmap.matrix=exprs(eset), project, cex.label=NULL, heatmap.scale=2, plot.symbols=T, main=NULL, device=NULL, dist.method="correlation", clust.method="ward", cor.method="pearson", col.order=NULL, row.order=NULL, sample.names=sampleNames(eset), heatmap.cols=c("green", "black", "red"), m.pal=NULL, single.colour=F, phenotypes=NULL, pheno.colours=NULL, plot.sample.names=T, pb.height=0.1)
{
	require(Biobase)
	require(lattice)
	require(latticeExtra)
	require(grid)
	require(marray)
	
	heatmap.matrix.name <- rev(unlist(strsplit(deparse(match.call()$heatmap.matrix), "\\$")))[1]
	
	if(identical(chrom.order, T)){
		if(is.null(fData(eset)$chrom)) stop("Object not mapped\nRecord chrom/start/end in feature Data\n")
		eset <- eset[order(fData(eset)$chrom, fData(eset)$start),]
		}
	
	
	if(!is.null(cluster) && !identical(cluster,F)){
		
		if(dist.method == "correlation"){
			atr <- hclust(dist(1-cor(cluster, method=cor.method)), method=clust.method)
			if(identical(chrom.order, F)) gtr <- hclust(dist(1-cor(t(cluster), method=cor.method)), method=clust.method)
			}else {
			atr <- hclust(dist(t(cluster), method=dist.method), method=clust.method)
			if(identical(chrom.order, F)) gtr <- hclust(dist(cluster, method=dist.method), method=clust.method)
			}
		atr.dg <- as.dendrogram(atr)
		if(identical(chrom.order, F)) {
			gtr.dg <- as.dendrogram(gtr)
			row.order <- gtr$order
			gtr.grob <- dendrogramGrob(gtr.dg, side="right", size=3)
			}
		col.order <- atr$order
		atr.grob <- dendrogramGrob(atr.dg, side="top", size=3)
		
	}
	
	if(is.null(row.order)) row.order <- 1:nrow(fData(eset))
	if(is.null(col.order)) col.order <- 1:nrow(eset)
	
	if(!is.null(phenotypes)){
		
		#	phenotypes
		pheno.labels <- as.matrix(phenotypes)
		
		pheno.matrix <- matrix(1, nrow(pheno.labels), ncol(pheno.labels))
	
		phenos <- unique(pheno.labels[!is.na(pheno.labels[,1]),1])
		phenos <- phenos[order(phenos)]
		if(ncol(pheno.labels) > 1)
		{
			for(i in 2:ncol(pheno.labels))
			{
			phenos.i <- unique(pheno.labels[!is.na(pheno.labels[,i]),i])
			phenos.i <- phenos.i[order(phenos.i)]
			phenos <- c(phenos,phenos.i[which(!is.element(phenos.i,phenos))])
			}
		}
	
	pheno.matrix <- matrix(as.numeric(match(pheno.labels, phenos)),nrow(pheno.labels), ncol(pheno.labels))
	
	pheno.matrix <- pheno.matrix[col.order,,drop=F]
	pheno.labels <- pheno.labels[col.order,,drop=F]
	
	leg.colours <- unique(as.vector(pheno.matrix[!is.na(pheno.matrix)]))
	leg.labels <- unique(as.vector(pheno.labels[!is.na(pheno.labels)]))
	
	display.brewer.all()
	if(is.null(pheno.colours)) {leg.pal <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2"), brewer.pal(8,"Dark2"))[1:(length(leg.colours))]}else{
		leg.pal <- pheno.colours}
	leg.pal <- leg.pal[1:length(leg.colours)]	
	
	if(ncol(pheno.matrix) ==1 ){pheno.matrix <- cbind(pheno.matrix, pheno.matrix)}
	pheno.mat.grob <- levelplot(pheno.matrix, col.regions=leg.pal, at=0:length(leg.colours)+0.5,colorkey=F)
	
	# reverse legend labels - ie alphabetical downwards
	leg.mat <- as.matrix(leg.colours[rev(order(leg.colours))],length(leg.colours),1)
	leg.grob <- levelplot(t(leg.mat), col.regions=leg.pal, at=0:length(leg.colours)+0.5, colorkey=F)
	leg.labels <- leg.labels[rev(order(leg.colours))]
		
		}
	
	heat.matrix <- heatmap.matrix[rev(row.order),col.order]
	maxM <- (as.integer(max(abs(heat.matrix[!is.na(heat.matrix)]))))+1
	
	if(!is.null(heatmap.scale)){
		if(maxM > heatmap.scale) {
		heat.matrix[heat.matrix > heatmap.scale] <- heatmap.scale
		heat.matrix[heat.matrix < -heatmap.scale] <- -heatmap.scale
		}
		}

	if(heatmap.matrix.name == "GL"){
		split.levels=seq(1:6)-3.5
		m.pal <- maPalette(low = "blue", mid="white", high = "red", k =6)
		}else{
		split.levels <- seq(min(heat.matrix, na.rm=T)-0.1,max(heat.matrix, na.rm=T)+0.1,by=0.1)
		if(is.null(m.pal)){
			if(length(heatmap.cols) == 2){
				m.pal <- maPalette(low = heatmap.cols[1], high = heatmap.cols[2], k =length(split.levels))}else{
				m.pal <- maPalette(low = heatmap.cols[1], high = heatmap.cols[3], mid=heatmap.cols[2], k =length(split.levels))}
		}
	}
	 
	heat.mat <- as.matrix(heat.matrix, nrow(heat.matrix), ncol(heat.matrix))
	
	if(single.colour==T){heat.mat.grob <- levelplot(t(heat.mat), axes=F, col.regions=c("white", heatmap.cols[1]), at=c(0,0.5,1))}else{
		heat.mat.grob <- levelplot(t(heat.mat), axes=F, col.regions=m.pal, at=split.levels)
		}
         
	s.labels <- sample.names[col.order]
	
	
	if(length(device) !=0){
	if(!is.element("heatmaps", list.files())) dir.create("heatmaps")
			
	if(device == "quartz") {quartz(title=main,width=6, height=4)}
	if(device == "PNG")
		{png(file=paste("heatmaps/",project,".heatmap.png",sep=""),width=600, height=900, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("heatmaps/",project,".heatmap.jpeg",sep=""),width=600, height=900, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("heatmaps/",project,".heatmap.pdf",sep=""),width=6, height=9, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("heatmaps/",project,".heatmap.ps",sep=""),width=6, height=9, pointsize=12)}}
		
	grid.newpage()

	if(!is.null(cluster) && !identical(cluster,F)){
	
	if(!is.null(phenotypes)) {
		phenobar.height <- pb.height
		
		
		if(ncol(pheno.labels) ==1 ) {phenobar.legend.height <- phenobar.height/2*ncol(pheno.labels)}else{
			phenobar.legend.height <- phenobar.height/ncol(pheno.labels)}
			
		phenobar.legend <-
   			viewport(x=0 + phenobar.legend.height/2, y = unit(0.8-phenobar.height+phenobar.legend.height/2,  "npc"),
   			width = 1, height = unit(phenobar.legend.height*length(leg.labels), "npc"),
			layout = # necessary to fix aspect ratio
			grid.layout(1, 4,
			widths = unit(c(0.2,0.6,phenobar.legend.height, 0.2-phenobar.legend.height),"npc"),
			heights = rep(phenobar.legend.height*length(leg.labels),3),
			respect = TRUE),
			just = c("left", "bottom"),
			name = 'phenobar.legend')
		
		
		
		pushViewport(phenobar.legend)
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3,
                      xscale = leg.grob$x.limits,
                      yscale = leg.grob$y.limits))
		do.call("panel.levelplot", trellis.panelArgs(leg.grob, 1))
		popViewport()
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4,
                      xscale = leg.grob$x.limits,
                      yscale = leg.grob$y.limits))
		
		leg.fontsize <- min(8,(300/length(leg.labels)))
		for(i in 1:length(leg.labels)){
		grid.text(leg.labels[i],x=0.1, y=(i/length(leg.labels)-0.5/length(leg.labels)), 		just="left",gp=gpar(fontsize=leg.fontsize))}
		popViewport()
		
		popViewport()

		
			phenobar.vp <-
   		viewport(x=0, y = unit(0.8-phenobar.height,  "npc"),
   		width = 1, height = phenobar.height,
		layout = # necessary to fix aspect ratio
		grid.layout(1, 3,
			widths = unit(c(0.2,0.6,0.2),"npc"),
			heights = rep(phenobar.height,3),
			respect = TRUE),
			just = c("left", "bottom"),
			name = 'phenobar.vp',
			xscale = pheno.mat.grob$x.limits,
            yscale = pheno.mat.grob$y.limits)
		
			
		pushViewport(phenobar.vp)
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
		
		phenotypes.fontsize <- min(8,(300/ncol(pheno.labels)))
		for(i in 1:ncol(pheno.labels)){
		grid.text(names(phenotypes)[i],x=0.9, y=(i/ncol(pheno.labels)-0.5/ncol(pheno.labels)), 		just=c("right","centre"),gp=gpar(fontsize=phenotypes.fontsize))}
		popViewport()
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                      xscale = pheno.mat.grob$x.limits,
                      yscale = pheno.mat.grob$y.limits))
   		grid.rect(gp = gpar(col = 'black', lwd=1))
   		do.call("panel.levelplot", trellis.panelArgs(pheno.mat.grob, 1))
		
		#	plot horoizontal lines
	#	line.breaks <- seq(0,1,by=1/ncol(pheno.labels))
	#	for(i in 1:length(line.breaks)) grid.lines(y=line.breaks[i], gp = gpar(col = 'black', lwd=1))

		popViewport()
		popViewport()
		
				
		}else{phenobar.height <- 0}
		
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.9,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "bottom"))
	pushViewport(title.vp)
	grid.text(main, just="top")
	popViewport()
	
	atr.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.8,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(3, "lines"),
	just = c("left", "bottom"))
	pushViewport(atr.vp)
	grid.draw(atr.grob)
	popViewport()
	
	if(plot.sample.names==T){
	sample.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "top"))
	pushViewport(sample.vp)
	for(i in 1:length(s.labels)){
		grid.text(s.labels[i],x=(i/length(s.labels)-0.5/length(s.labels)),
		y=0.95, just=c("right","centre"), gp=gpar(cex=min(0.5,15/length(s.labels))), rot=90)}
	popViewport()
	}
	
	symbol.vp <- viewport(x = unit(0.8, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.2, "npc"),
	height = unit(0.7-phenobar.height, "npc"),
	just = c("left", "bottom"))
	pushViewport(symbol.vp)
	if(identical(chrom.order, F)) {
		for(i in 1:nrow(fData(eset))){
		if(identical(plot.symbols,T) && nrow(heat.mat) < 100){
			grid.text(rev(fData(eset)$symbol[row.order])[i],x=0.05,
			y=(i/length(fData(eset)$symbol)-0.5/length(fData(eset)$symbol)),
			just=c("left","centre"), gp=gpar(cex=25/nrow(heat.mat)))
		}}}else{
			chrom.lengths <- rev(as.vector(table(fData(eset)$chrom))/length(fData(eset)$chrom))
			i <- 1
			chrom.start <- 0
			while(i <= length(chrom.lengths)){
			grid.rect(0.7,chrom.start,width=0.2,height=chrom.lengths[i], just=c("centre", "bottom"), gp=gpar(fill=i%%2))
			chroms <- unique(fData(eset)$chrom)
			chroms[chroms == 23] <- "X"
			chroms[chroms == 24] <- "Y"
			grid.text(rev(chroms)[i], 0.5,chrom.start+(chrom.lengths[i]/2), gp=gpar(cex=min(1,12/length(chroms))))
			chrom.start <- chrom.start+chrom.lengths[i]
			i <- i+1
			}}
	popViewport()
	
	gtr.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(3, "lines"),
	height = unit(0.7-phenobar.height, "npc"),
	just = c("left", "top"),
	angle=180)
	pushViewport(gtr.vp)

	if(identical(chrom.order, F)) {grid.draw(gtr.grob)}
	popViewport()
	
	matrix.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(0.7-phenobar.height, "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,ncol(heat.mat)+0.5),
	yscale=c(0.5,nrow(heat.mat)+0.5))
	pushViewport(matrix.vp)
	grid.rect()
	do.call("panel.levelplot", trellis.panelArgs(heat.mat.grob, 1))
	popViewport()

	}
	
	if(is.null(cluster)|identical(cluster,F)){
		
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.9,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "bottom"))
	pushViewport(title.vp)
	grid.text(main, just="top")
	popViewport()
	
	if(plot.sample.names==T){
	sample.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "top"))
	pushViewport(sample.vp)
	for(i in 1:length(s.labels)){
		grid.text(s.labels[i],x=(i/length(s.labels)-0.5/length(s.labels)),
		y=0.95, just=c("right","centre"), gp=gpar(cex=min(0.5,15/length(s.labels))), rot=90)}
	popViewport()
	}
	
	symbol.vp <- viewport(x = unit(0.8, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.2, "npc"),
	height = unit(0.8, "npc"),
	just = c("left", "bottom"))
	pushViewport(symbol.vp)
	if(identical(chrom.order, F)) {
		for(i in 1:nrow(fData(eset))){
		if(identical(plot.symbols,T) && nrow(heat.mat) < 100){
			grid.text(rev(fData(eset)$symbol[row.order])[i],x=0.05,
			y=(i/length(fData(eset)$symbol)-0.5/length(fData(eset)$symbol)),
			just=c("left","centre"), gp=gpar(cex=25/nrow(heat.mat)))
		}}}else{
			chrom.lengths <- rev(as.vector(table(fData(eset)$chrom))/length(fData(eset)$chrom))
			i <- 1
			chrom.start <- 0
			while(i <= length(chrom.lengths)){
			grid.rect(0.7,chrom.start,width=0.2,height=chrom.lengths[i], just=c("centre", "bottom"), gp=gpar(fill=i%%2))
			chroms <- unique(fData(eset)$chrom)
			chroms[chroms == 23] <- "X"
			chroms[chroms == 24] <- "Y"
			grid.text(rev(chroms)[i], 0.5,chrom.start+(chrom.lengths[i]/2), gp=gpar(cex=min(1,12/length(chroms))))
			chrom.start <- chrom.start+chrom.lengths[i]
			i <- i+1
			}}
	
	popViewport()
	
	matrix.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.8, "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,ncol(heat.mat)+0.5),
	yscale=c(0.5,nrow(heat.mat)+0.5))
	pushViewport(matrix.vp)
	grid.rect()
	do.call("panel.levelplot", trellis.panelArgs(heat.mat.grob, 1))
	popViewport()
		
		
		
		
	}
	
	if(length(device) != 0 && device != "quartz") dev.off()
	
}


cghHeatmap <- function(cgh, cluster="GL", heatmap.matrix="GL", project, cex.labels=1, heatmap.scale=NULL, plot.symbols=T, main=NULL, device=NULL, dist.method="correlation", clust.method="ward", cor.method="pearson", phenotypes=NULL, pheno.colours=NULL, plot.sample.names=T, plot.matrix=T)
{
	require(Biobase)
	require(lattice)
	require(latticeExtra)
	require(grid)
	require(marray)
	
	cgh <- cgh[order(fData(cgh)$chrom, fData(cgh)$start),]
	
	heatmap.matrix.name <- heatmap.matrix
	heatmap.matrix <- assayDataElement(cgh, heatmap.matrix)
	
	if(!is.null(cluster) && !identical(cluster,F)){
		
		cluster <- assayDataElement(cgh, cluster)
		
		if(dist.method == "correlation"){
			atr <- hclust(as.dist(1-cor(cluster, method=cor.method)), method=clust.method)
			}else{atr <- hclust(dist(t(cluster), method=dist.method), method=clust.method)}
		atr.dg <- as.dendrogram(atr)
		col.order <- atr$order
		atr.grob <- dendrogramGrob(atr.dg, side="top", size=3)
		
	}else{col.order <- 1:ncol(cgh)}
	
	if(!is.null(phenotypes)){
		
		#	phenotypes
		pheno.labels <- as.matrix(phenotypes)
		
		pheno.matrix <- matrix(1, nrow(pheno.labels), ncol(pheno.labels))
	
		phenos <- unique(pheno.labels[!is.na(pheno.labels[,1]),1])
		phenos <- phenos[order(phenos)]
		if(ncol(pheno.labels) > 1)
		{
			for(i in 2:ncol(pheno.labels))
			{
			phenos.i <- unique(pheno.labels[!is.na(pheno.labels[,i]),i])
			phenos.i <- phenos.i[order(phenos.i)]
			phenos <- c(phenos,phenos.i[which(!is.element(phenos.i,phenos))])
			}
		}
	
		pheno.matrix <- matrix(as.numeric(match(pheno.labels, phenos)),nrow(pheno.labels), ncol(pheno.labels))
	
		pheno.matrix <- pheno.matrix[col.order,,drop=F]
		pheno.labels <- pheno.labels[col.order,,drop=F]
	
		leg.colours <- unique(as.vector(pheno.matrix[!is.na(pheno.matrix)]))
		leg.labels <- unique(as.vector(pheno.labels[!is.na(pheno.labels)]))
	
	
		if(is.null(pheno.colours)) {leg.pal <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2"), brewer.pal(8,"Dark2"))[1:(length(leg.colours))]}else{
		leg.pal <- pheno.colours}
		leg.pal <- leg.pal[1:length(leg.colours)]	
	
		if(ncol(pheno.matrix) ==1 ){pheno.matrix <- cbind(pheno.matrix, pheno.matrix)}
		pheno.mat.grob <- levelplot(pheno.matrix, col.regions=leg.pal, at=0:length(leg.colours)+0.5,colorkey=F)
	
		# reverse legend labels - ie alphabetical downwards
		leg.mat <- as.matrix(leg.colours[rev(order(leg.colours))],length(leg.colours),1)
		leg.grob <- levelplot(t(leg.mat), col.regions=leg.pal, at=0:length(leg.colours)+0.5, colorkey=F)
		leg.labels <- leg.labels[rev(order(leg.colours))]
		
		}
	
	heat.matrix <- heatmap.matrix[rev(1:nrow(cgh)),col.order]
	maxM <- (as.integer(max(abs(heat.matrix[!is.na(heat.matrix)]))))+1
	
	if(heatmap.matrix.name == "GL") heatmap.scale <- 2
	if(is.null(heatmap.scale)) heatmap.scale <- min(abs(range(heat.matrix)))
	
	if(maxM > heatmap.scale) {
		heat.matrix[heat.matrix > heatmap.scale] <- heatmap.scale
		heat.matrix[heat.matrix < -heatmap.scale] <- -heatmap.scale
		}

	if(heatmap.matrix.name == "GL"){
		split.levels=seq(1:6)-3.5
		m.pal <- maPalette(low = "blue", mid="white", high = "red", k =6)
		}else{
			split.levels <- seq(min(heat.matrix, na.rm=T),max(heat.matrix, na.rm=T),by=0.1)
			m.pal <- maPalette(low = "blue", mid="white", high = "red", k =length(split.levels))
		}
	 
	heat.mat <- as.matrix(heat.matrix, nrow(heat.matrix), ncol(heat.matrix))
	heat.mat.grob <- levelplot(t(heat.mat), axes=F, col.regions=m.pal, at=split.levels)
	s.labels <- sampleNames(cgh)[col.order]
	
	colorkey.mat <- matrix(split.levels, nrow=2, ncol=length(split.levels), byrow=T)
	colorkey.grob <- levelplot(colorkey.mat[,-1], axes=F, col.regions=m.pal, at=split.levels)
	
	if(length(device) !=0){
	if(!is.element("heatmaps", list.files())) dir.create("heatmaps")
			
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("heatmaps/",project,".heatmap.png",sep=""),width=600, height=900, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("heatmaps/",project,".heatmap.jpeg",sep=""),width=600, height=900, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("heatmaps/",project,".heatmap.pdf",sep=""),width=6, height=9, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("heatmaps/",project,".heatmap.ps",sep=""),width=6, height=9, pointsize=12)}}
	
	grid.newpage()
	
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.9,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "bottom"))
	pushViewport(title.vp)
	grid.text(main, just="top", gp=gpar(cex=1))
	popViewport()
	
	colorkey.vp <- viewport(x = unit(0.825, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.025,  "npc"),
	height = unit(0.1,  "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,(ncol(colorkey.mat)/2)-0.5),
	yscale=c(0.5,ncol(colorkey.mat)-0.5))
	pushViewport(colorkey.vp)
	do.call("panel.levelplot", trellis.panelArgs(colorkey.grob, 1))
	popViewport()
	
	colorkey.leg.vp <- viewport(x = unit(0.85, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.1,  "npc"),
	height = unit(0.1,  "npc"),
	just = c("left", "bottom"))
	pushViewport(colorkey.leg.vp)
	
		if(heatmap.matrix.name == "GL"){
			glad.leg <- c("Del", "Loss", "NC", "Gain", "Amp")
			for(i in 1:5) grid.text(glad.leg[i],y=(i/5-0.5/5),
		x=0.05, just=c("left","centre"), gp=gpar(cex=min(0.5,15/5)))
			}else{
				for(i in 1:length(split.levels)){
			grid.text(split.levels[i],y=(i/length(split.levels)-0.5/length(split.levels)),
		x=0.05, just=c("left","centre"), gp=gpar(cex=min(0.5,15/length(split.levels))))
			}
		}
	popViewport()
	
	if(!is.null(phenotypes)) {
		phenobar.height <- 0.05
		
		
		if(ncol(pheno.labels) ==1 ) {phenobar.legend.height <- phenobar.height/2*ncol(pheno.labels)}else{
			phenobar.legend.height <- phenobar.height/ncol(pheno.labels)}
			
		phenobar.legend <-
   			viewport(x=0 + phenobar.legend.height/2, y = unit(0.8-phenobar.height+phenobar.legend.height/2,  "npc"),
   			width = 1, height = unit(phenobar.legend.height*length(leg.labels), "npc"),
			layout = # necessary to fix aspect ratio
			grid.layout(1, 4,
			widths = unit(c(0.2,0.6,phenobar.legend.height, 0.2-phenobar.legend.height),"npc"),
			heights = rep(phenobar.legend.height*length(leg.labels),3),
			respect = TRUE),
			just = c("left", "bottom"),
			name = 'phenobar.legend')
		
		
		
		pushViewport(phenobar.legend)
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3,
                      xscale = leg.grob$x.limits,
                      yscale = leg.grob$y.limits))
		do.call("panel.levelplot", trellis.panelArgs(leg.grob, 1))
		popViewport()
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4,
                      xscale = leg.grob$x.limits,
                      yscale = leg.grob$y.limits))
		
		leg.fontsize <- min(8,(300/length(leg.labels)))
		for(i in 1:length(leg.labels)){
		grid.text(leg.labels[i],x=0.1, y=(i/length(leg.labels)-0.5/length(leg.labels)), 		just="left",gp=gpar(fontsize=leg.fontsize))}
		popViewport()
		
		popViewport()

		
		phenobar.vp <-
   		viewport(x=0, y = unit(0.8-phenobar.height,  "npc"),
   		width = 1, height = phenobar.height,
		layout = # necessary to fix aspect ratio
		grid.layout(1, 3,
			widths = unit(c(0.2,0.6,0.2),"npc"),
			heights = rep(phenobar.height,3),
			respect = TRUE),
			just = c("left", "bottom"),
			name = 'phenobar.vp',
			xscale = pheno.mat.grob$x.limits,
            yscale = pheno.mat.grob$y.limits)
		
			
		pushViewport(phenobar.vp)
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
		
		phenotypes.fontsize <- min(8,(300/ncol(pheno.labels)))
		for(i in 1:ncol(pheno.labels)){
		grid.text(names(phenotypes)[i],x=0.9, y=(i/ncol(pheno.labels)-0.5/ncol(pheno.labels)), 		just=c("right","centre"),gp=gpar(fontsize=phenotypes.fontsize))}
		popViewport()
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                      xscale = pheno.mat.grob$x.limits,
                      yscale = pheno.mat.grob$y.limits))
   		grid.rect(gp = gpar(col = 'black', lwd=1))
   		do.call("panel.levelplot", trellis.panelArgs(pheno.mat.grob, 1))
		
	#	plot horoizontal lines
	#	line.breaks <- seq(0,1,by=1/ncol(pheno.labels))
	#	for(i in 1:length(line.breaks)) grid.lines(y=line.breaks[i], gp = gpar(col = 'black', lwd=1))

		popViewport()
		popViewport()
		
				
		}else{phenobar.height <- 0}
		
		
	symbol.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.2, "npc"),
	height = unit(0.7-phenobar.height, "npc"),
	just = c("left", "bottom"))
	pushViewport(symbol.vp)
	
			chrom.lengths <- rev(as.vector(table(fData(cgh)$chrom))/length(fData(cgh)$chrom))
			i <- 1
			chrom.start <- 0
			while(i <= length(chrom.lengths)){
			grid.rect(0.7,chrom.start,width=0.2,height=chrom.lengths[i], just=c("centre", "bottom"), gp=gpar(fill=i%%2))
			chroms <- unique(fData(cgh)$chrom)
			chroms[chroms == 23] <- "X"
			chroms[chroms == 24] <- "Y"
			grid.text(rev(chroms)[i], 0.5,chrom.start+(chrom.lengths[i]/2), gp=gpar(cex=min(1,12/length(chroms))))
			chrom.start <- chrom.start+chrom.lengths[i]
			i <- i+1
			}
	popViewport()
	
	if(plot.sample.names==T){
	sample.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "top"))
	pushViewport(sample.vp)
	for(i in 1:length(s.labels)){
		grid.text(s.labels[i],x=(i/length(s.labels)-0.5/length(s.labels)),
		y=0.95, just=c("right","centre"), gp=gpar(cex=cex.labels*min(1,30/length(s.labels))), rot=90)}
	popViewport()
	}
	
	matrix.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(0.7-phenobar.height, "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,ncol(heat.mat)+0.5),
	yscale=c(0.5,nrow(heat.mat)+0.5))
	pushViewport(matrix.vp)
	grid.rect()
	if(plot.matrix == T) do.call("panel.levelplot", trellis.panelArgs(heat.mat.grob, 1))
	popViewport()
	
	if(!is.null(cluster) && !identical(cluster,F)){
	
	atr.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.8,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(3, "lines"),
	just = c("left", "bottom"))
	pushViewport(atr.vp)
	grid.draw(atr.grob)
	popViewport()

	}
	if(length(device) != 0 && device != "quartz") dev.off()
}


plotPhenoBar<- function(eset, cluster=exprs(eset), main=NULL, dist.method="correlation", clust.method="ward", cor.method="pearson", phenotypes=NULL, pheno.colours=NULL, plot.sample.names=T, cex.labels=1, device=NULL, project, col.order=NULL, plot.legend=T)
#	Function to plot phenotype infomation to accompany a dendrogram in a tree order file
{

	require(Biobase)
	require(lattice)
	require(latticeExtra)
	require(grid)
	require(marray)

	if(!is.null(cluster) && !identical(cluster,F)){
		
		if(dist.method == "correlation"){
			atr <- hclust(as.dist(1-cor(cluster, method=cor.method)), method=clust.method)
			}else{atr <- hclust(dist(t(cluster), method=dist.method), method=clust.method)}
		atr.dg <- as.dendrogram(atr)
		col.order <- atr$order
		atr.grob <- dendrogramGrob(atr.dg, side="top", size=6)
		
	}
	
	if(is.null(col.order)) col.order <- 1:ncol(eset)
	
	pheno.labels <- as.matrix(phenotypes)
		
	pheno.matrix <- matrix(1, nrow(pheno.labels), ncol(pheno.labels))
	
	phenos <- unique(pheno.labels[!is.na(pheno.labels[,1]),1])
	phenos <- phenos[order(phenos)]
	if(ncol(pheno.labels) > 1)
		{
			for(i in 2:ncol(pheno.labels))
			{
			phenos.i <- unique(pheno.labels[!is.na(pheno.labels[,i]),i])
			phenos.i <- phenos.i[order(phenos.i)]
			phenos <- c(phenos,phenos.i[which(!is.element(phenos.i,phenos))])
			}
		}
	
	pheno.matrix <- matrix(as.numeric(match(pheno.labels, phenos)),nrow(pheno.labels), ncol(pheno.labels))
	
	pheno.matrix <- pheno.matrix[col.order,,drop=F]
	pheno.labels <- pheno.labels[col.order,,drop=F]
	
	leg.colours <- unique(as.vector(pheno.matrix[!is.na(pheno.matrix)]))
	leg.labels <- unique(as.vector(pheno.labels[!is.na(pheno.labels)]))
	
	
	if(is.null(pheno.colours)) {leg.pal <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2"), brewer.pal(8,"Dark2"))[1:(length(leg.colours))]}else{
		leg.pal <- pheno.colours}
	leg.pal <- leg.pal[1:length(leg.colours)]	
	
	if(ncol(pheno.matrix) ==1 ){pheno.matrix <- cbind(pheno.matrix, pheno.matrix)}
	pheno.mat.grob <- levelplot(pheno.matrix, col.regions=leg.pal, at=0:length(leg.colours)+0.5,colorkey=F)
	
	leg.mat <- as.matrix(leg.colours[rev(order(leg.colours))],length(leg.colours),1)
	leg.grob <- levelplot(t(leg.mat), col.regions=leg.pal, at=0:length(leg.colours)+0.5, colorkey=F)
	leg.labels <- leg.labels[rev(order(leg.colours))]
	
	if(length(device) !=0){
	if(!is.element("phenobars", list.files())) dir.create("phenobars")
			
	if(device == "quartz") {quartz(title=main,width=6, height=4)}
	if(device == "PNG")
		{png(file=paste("phenobars/",project,".phenobar.png",sep=""),width=900, height=600, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("phenobars/",project,".phenobar.jpeg",sep=""),width=900, height=600, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("phenobars/",project,".phenobar.pdf",sep=""),width=9, height=6, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("phenobars/",project,".phenobar.ps",sep=""),width=9, height=6, pointsize=12)}}
	
	if(!is.null(cluster) && !identical(cluster,F)){	
	grid.newpage()
	
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.8,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "bottom"))
	pushViewport(title.vp)
	grid.text(main, just="top", gp=gpar(cex=2))
	popViewport()

	phenobar.height <- 0.2
		
	if(ncol(pheno.labels) ==1 ) {phenobar.legend.height <- phenobar.height/2*ncol(pheno.labels)}else{
		phenobar.legend.height <- phenobar.height/ncol(pheno.labels)}
			
	phenobar.vp <- viewport(x=0, y = unit(0.5-phenobar.height,  "npc"),
   		width = 1, height = phenobar.height,
		layout = # necessary to fix aspect ratio
		grid.layout(1, 3,
			widths = unit(c(0.15,0.6, 0.25),"npc"),
			heights = rep(phenobar.height,3),
			respect = TRUE),
			just = c("left", "bottom"),
			name = 'phenobar.vp',
			xscale = pheno.mat.grob$x.limits,
            yscale = pheno.mat.grob$y.limits)
		
	pushViewport(phenobar.vp)
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
		
	phenotypes.fontsize <- min(8,(300/ncol(pheno.labels)))
	for(i in 1:ncol(pheno.labels)){
		grid.text(names(phenotypes)[i],x=0.95, y=(i/ncol(pheno.labels)-0.5/ncol(pheno.labels)), 		just=c("right","centre"),gp=gpar(fontsize=phenotypes.fontsize))}
	popViewport()
		
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                      xscale = pheno.mat.grob$x.limits,
                      yscale = pheno.mat.grob$y.limits))
	do.call("panel.levelplot", trellis.panelArgs(pheno.mat.grob, 1))
   	grid.rect(gp = gpar(col = 'black', lwd=1))
	#	plot horoizontal lines
	#	line.breaks <- seq(0,1,by=1/ncol(pheno.labels))
	#	for(i in 1:length(line.breaks)) grid.lines(y=line.breaks[i], gp = gpar(col = 'black', lwd=1))
	popViewport()
	popViewport()
	
	if(plot.legend==T){
	legend.vp <- viewport(x=0.8, y = unit(0.5,  "npc"),
   		width = 0.2, height = min(phenobar.legend.height*length(leg.labels),0.8),
		layout = # necessary to fix aspect ratio
		grid.layout(1, 2,
			widths = unit(c(0.2,0.8),"npc"),
			heights = rep(phenobar.legend.height*length(leg.labels),2),
			respect = TRUE),
			just = c("left", "center"),
			name = 'legend.vp',
			xscale = pheno.mat.grob$x.limits,
            yscale = pheno.mat.grob$y.limits)
	pushViewport(legend.vp)
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1,
                      xscale = leg.grob$x.limits,
                      yscale = leg.grob$y.limits))
	do.call("panel.levelplot", trellis.panelArgs(leg.grob, 1))
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                      xscale = leg.grob$x.limits,
                      yscale = leg.grob$y.limits))
		
	leg.fontsize <- min(8,(300/length(leg.labels)))
	for(i in 1:length(leg.labels)){
		grid.text(leg.labels[i],x=0.1, y=(i/length(leg.labels)-0.5/length(leg.labels)), 		just="left",gp=gpar(fontsize=leg.fontsize))}
	popViewport()
	popViewport()
	}
			
	atr.vp <- viewport(x = unit(0.15, "npc"),
		y = unit(0.5,  "npc"),
		width = unit(0.6, "npc"),
		height = unit(6, "lines"),
		just = c("left", "bottom"))
	pushViewport(atr.vp)
	grid.draw(atr.grob)
	popViewport()
		
	if(plot.sample.names){
		sample.vp <- viewport(x = unit(0.15, "npc"),
		y = unit(0.5-phenobar.height,  "npc"),
		width = unit(0.6, "npc"),
		height = unit(0.5-phenobar.height, "npc"),
		just = c("left", "top"))
	pushViewport(sample.vp)
	s.labels <- sampleNames(eset)[col.order]
	for(i in 1:length(sampleNames(eset))){
			grid.text(s.labels[i],x=(i/length(sampleNames(eset))-0.5/length(sampleNames(eset))),
			y=0.95, just=c("right","centre"), gp=gpar(cex=cex.labels*min(0.5,15/length(sampleNames(eset)))), rot=90)}
	popViewport()
	}
	
	}else{
	
	grid.newpage()
 	
 	title.vp <- viewport(x = unit(0, "npc"),
 	y = unit(0.8,  "npc"),
 	width = unit(1, "npc"),
 	height = unit(0.1, "npc"),
 	just = c("left", "bottom"))
 	pushViewport(title.vp)
 	grid.text(main, just="top", gp=gpar(cex=2))
 	popViewport()
 
 	phenobar.height <- 0.2
 		
 	if(ncol(pheno.labels) ==1 ) {phenobar.legend.height <- phenobar.height/2*ncol(pheno.labels)}else{
 		phenobar.legend.height <- phenobar.height/ncol(pheno.labels)}
 			
 	phenobar.vp <- viewport(x=0, y = unit(0.5,  "npc"),
 	    		width = 1, height = phenobar.height,
 		layout = # necessary to fix aspect ratio
 		grid.layout(1, 3,
 			widths = unit(c(0.15,0.6, 0.25),"npc"),
 			heights = rep(phenobar.height,3),
 			respect = TRUE),
 			just = c("left", "centre"),
 			name = 'phenobar.vp',
 			xscale = pheno.mat.grob$x.limits,
            yscale = pheno.mat.grob$y.limits)
 		
 	pushViewport(phenobar.vp)
 	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
 		
 	phenotypes.fontsize <- min(8,(300/ncol(pheno.labels)))
 	for(i in 1:ncol(pheno.labels)){
 		grid.text(names(phenotypes)[i],x=0.95, y=(i/ncol(pheno.labels)-0.5/ncol(pheno.labels)), 		just=c("right","centre"),gp=gpar(fontsize=phenotypes.fontsize))}
 	popViewport()
 		
 	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                      xscale = pheno.mat.grob$x.limits,
                      yscale = pheno.mat.grob$y.limits))
 	do.call("panel.levelplot", trellis.panelArgs(pheno.mat.grob, 1))
    	grid.rect(gp = gpar(col = 'black', lwd=1))
 	#	plot horoizontal lines
 	#	line.breaks <- seq(0,1,by=1/ncol(pheno.labels))
 	#	for(i in 1:length(line.breaks)) grid.lines(y=line.breaks[i], gp = gpar(col = 'black', lwd=1))
 	popViewport()
 	popViewport()
 	
 	if(plot.legend==T){
 	legend.vp <- viewport(x=0.8, y = unit(0.5,  "npc"),
    		width = 0.2, height = phenobar.legend.height*length(leg.labels),
 		layout = # necessary to fix aspect ratio
 		grid.layout(1, 2,
 			widths = unit(c(0.2,0.8),"npc"),
 			heights = rep(phenobar.legend.height*length(leg.labels),2),
 			respect = TRUE),
 			just = c("left", "center"),
 			name = 'legend.vp',
 			xscale = pheno.mat.grob$x.limits,
            yscale = pheno.mat.grob$y.limits)
 	pushViewport(legend.vp)
 	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1,
                     xscale = leg.grob$x.limits,
                     yscale = leg.grob$y.limits))
 	do.call("panel.levelplot", trellis.panelArgs(leg.grob, 1))
 	popViewport()
 	
 	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                   xscale = leg.grob$x.limits,
                   yscale = leg.grob$y.limits))
 		
 	leg.fontsize <- min(8,(300/length(leg.labels)))
 	for(i in 1:length(leg.labels)){
 		grid.text(leg.labels[i],x=0.1, y=(i/length(leg.labels)-0.5/length(leg.labels)), 		just="left",gp=gpar(fontsize=leg.fontsize))}
 	popViewport()
 	popViewport()
 	}
 		
 	if(plot.sample.names){
 		sample.vp <- viewport(x = unit(0.15, "npc"),
 		y = unit(0.5-phenobar.height,  "npc"),
 		width = unit(0.6, "npc"),
 		height = unit(phenobar.height, "npc"),
 		just = c("left", "center"))
 	pushViewport(sample.vp)
 	s.labels <- sampleNames(eset)[col.order]
 	for(i in 1:length(sampleNames(eset))){
 			grid.text(s.labels[i],x=(i/length(sampleNames(eset))-0.5/length(sampleNames(eset))),
 			y=0.95, just=c("right","centre"), gp=gpar(cex=cex.labels*min(0.5,15/length(sampleNames(eset)))), rot=90)}
 	popViewport()
 	}	
		
	}
	
	if(length(device) != 0 && device != "quartz") dev.off()
}



sampleDendrogram <- function(eset, cluster=exprs(eset), dist.method="correlation", clust.method="ward", cor.method="pearson", main=NULL, device=NULL, project, phenotypes=NULL, pheno.colours=NULL, cex.label=0.5, return.atr=T)
#	Function to plot phenotype infomation to accompany a dendrogram in a tree order file
{
	require(Biobase)
	require(Biobase)
	require(lattice)
	require(grid)
	require(latticeExtra)
	
	if(length(device) !=0){
	if(!is.element("dendrograms", list.files())) dir.create("dendrograms")
			
	if(device == "quartz") {quartz(title=main,width=12, height=6)}
	if(device == "PNG")
		{png(file=paste("dendrograms/",project,".dendrogram.png",sep=""),width=1600, height=800, pointsize=min(18,900/ncol(eset)))}
	if(device == "PDF")
		{pdf(file=paste("dendrograms/",project,".dendrogram.pdf",sep=""),width=8, height=4, pointsize=min(8,400/ncol(eset)))}
	if(device == "PS")
		{postscript(file=paste("dendrograms/",project,".dendrogram.ps",sep=""),width=6, height=4, pointsize=min(8,400/ncol(eset)))}}
	

	if(!is.null(cluster) && !identical(cluster,F))
			{
			cluster[is.na(cluster)] <- 0
			
			if(dist.method == "correlation")
				{
					atr <- hclust(dist(1-cor(cluster, method=cor.method)), method=clust.method)
				}else{
					atr <- hclust(dist(t(cluster), method=dist.method), method=clust.method)
					}

		eset <- eset[,atr$order]
		atr.dg <- as.dendrogram(atr)
		atr.grob <- dendrogramGrob(atr.dg, side="top", size=6)
		if(!is.null(phenotypes)){phenotypes <- data.frame(phenotypes)[atr$order,]}
		}

	if(is.null(phenotypes)){
		atr$call <- NULL
		plclust(atr, main=main,axes=T,hang=-1,
		ylab=paste(dist.method," distance"),
		xlab=paste(nrow(eset),"probes - Cluster method = ",atr$method,"with",dist.method, "distance"))
		}else{
	
	pheno.labels <- as.matrix(phenotypes)
	pheno.matrix <- matrix(1, nrow(pheno.labels), ncol(pheno.labels))
	
	phenos <- unique(pheno.labels[!is.na(pheno.labels[,1]),1])
	phenos <- phenos[order(phenos)]
	
	if(ncol(pheno.labels) > 1)
	{
		for(i in 2:ncol(pheno.labels))
		{
		phenos.i <- unique(pheno.labels[!is.na(pheno.labels[,i]),i])
		phenos.i <- phenos.i[order(phenos.i)]
		phenos <- c(phenos,phenos.i[which(!is.element(phenos.i,phenos))])
		}
	}
	
	pheno.matrix <- matrix(as.numeric(match(pheno.labels, phenos)),nrow(pheno.labels), ncol(pheno.labels))
	
	leg.colours <- unique(as.vector(pheno.matrix[!is.na(pheno.matrix)]))
	leg.labels <- unique(as.vector(pheno.labels[!is.na(pheno.labels)]))
	
	
	if(is.null(pheno.colours)) {m.pal <- rep(brewer.pal(12,"Paired"),3)[1:(length(leg.colours))]}else{
		m.pal <- pheno.colours}
	m.pal <- m.pal[1:length(leg.colours)]	
	

	if(ncol(phenotypes) == 1){pheno.matrix <- cbind(pheno.matrix, pheno.matrix)}
	
	pheno.mat.grob <- levelplot(pheno.matrix, col.regions=m.pal, at=0:length(leg.colours)+0.5,colorkey=F)
	
	leg.mat <- as.matrix(leg.colours[order(leg.colours)],length(leg.colours),1)
	leg.grob <- levelplot(t(leg.mat), col.regions=m.pal, at=0:length(leg.colours)+0.5, colorkey=F)
	leg.labels <- leg.labels[order(leg.colours)]
	
	
#	plot the components - GRID
	
	if(!is.null(cluster) && !identical(cluster,F)){
	
	grid.newpage()
	
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.9,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "center"))
	pushViewport(title.vp)
	grid.text(main, just="top")
	popViewport()
	
	atr.vp <- viewport(x = unit(0.25, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(6.1, "lines"),
	just = c("left", "bottom"))
	pushViewport(atr.vp)
	grid.draw(atr.grob)
	popViewport()
	
	phenobar.height <- (1.6/ncol(eset))*ncol(phenotypes)
	if(ncol(phenotypes) == 1){phenobar.height <- phenobar.height/2}
	
	phenobar.vp <- viewport(x = unit(0.25, "npc"),
	y = unit((0.5-phenobar.height),  "npc"),
	width = unit(0.6, "npc"),
	height = unit(phenobar.height, "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,nrow(pheno.matrix)+0.5),
	yscale=c(0.5,ncol(pheno.matrix)+0.5))
	pushViewport(phenobar.vp)
	do.call("panel.levelplot", trellis.panelArgs(pheno.mat.grob, 1))
	popViewport()
	
	phenobar.labels.vp <- viewport(x = unit(0.85, "npc"),
	y = unit((0.5-phenobar.height),  "npc"),
	width = unit(0.15, "npc"),
	height = unit(phenobar.height, "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,nrow(pheno.matrix)+0.5),
	yscale=c(0.5,ncol(pheno.matrix)+0.5))
	pushViewport(phenobar.labels.vp)
	if(ncol(phenotypes) == 1) {grid.text(names(phenotypes)[1], x=0.05, y=0.5, just=c("left","centre"))}else{
		for(i in 1:ncol(phenotypes)){
		grid.text(names(phenotypes)[i],x=0.05, y=(i/ncol(phenotypes)-0.5/ncol(phenotypes)), just="left", gp=gpar(cex=cex.label*1.5))}
		}
	popViewport()
	
	
	sample.vp <- viewport(x = unit(0.25, "npc"),
	y = unit(0.5-phenobar.height,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "top"))
	pushViewport(sample.vp)
	for(i in 1:ncol(eset)){
		grid.text(sampleNames(eset)[i],x=(i/ncol(eset)-0.5/ncol(eset)),
		y=0.95, just=c("right","centre"), gp=gpar(cex=cex.label), rot=90)}
	popViewport()
	
	legend.height <- (1.6/ncol(eset))*length(leg.labels)
	
	legend.vp <- viewport(x = unit(0.05, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(legend.height/(2*length(leg.labels)), "npc"),
	height = unit(min(legend.height,0.8), "npc"),
	just = c("left", "centre"),
	xscale = leg.grob$x.limits,
    yscale = leg.grob$y.limits)
	pushViewport(legend.vp)
	do.call("panel.levelplot", trellis.panelArgs(leg.grob, 1))
	popViewport()
	
	legend.labels.vp <- viewport(x = unit((legend.height/(2*length(leg.labels)))+0.05, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(0.19-(legend.height/(2*length(leg.labels))), "npc"),
	height = unit(min(legend.height,0.8), "npc"),
	just = c("left", "centre"),
	xscale = leg.grob$x.limits,
    yscale = leg.grob$y.limits)
	pushViewport(legend.labels.vp)
	for(i in 1:length(leg.labels)){
		grid.text(leg.labels[i],x=0.1, y=(i/length(leg.labels)-0.5/length(leg.labels)), just="left", gp=gpar(cex=cex.label*1.5))}
	popViewport()
	
	subtitle.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "top"))
	pushViewport(subtitle.vp)
	grid.text(paste(nrow(eset),"probes - Cluster method = ",atr$method,"with",dist.method, "distance"))
	popViewport()
	
	}else{
		
	grid.newpage()
	
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.9,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("centre", "center"))
	pushViewport(title.vp)
	grid.text(main, just="top")
	popViewport()
	
	phenobar.height <- (1.6/ncol(eset))*ncol(phenotypes)
	if(ncol(phenotypes) == 1){phenobar.height <- phenobar.height/2}
	
	phenobar.vp <- viewport(x = unit(0.25, "npc"),
	y = unit((0.7-phenobar.height),  "npc"),
	width = unit(0.6, "npc"),
	height = unit(phenobar.height, "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,nrow(pheno.matrix)+0.5),
	yscale=c(0.5,ncol(pheno.matrix)+0.5))
	pushViewport(phenobar.vp)
	do.call("panel.levelplot", trellis.panelArgs(pheno.mat.grob, 1))
	popViewport()
	
	phenobar.labels.vp <- viewport(x = unit(0.85, "npc"),
	y = unit((0.7-phenobar.height),  "npc"),
	width = unit(0.15, "npc"),
	height = unit(phenobar.height, "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,nrow(pheno.matrix)+0.5),
	yscale=c(0.5,ncol(pheno.matrix)+0.5))
	pushViewport(phenobar.labels.vp)
	if(ncol(phenotypes) == 1) {grid.text(names(phenotypes)[1], x=0.05, y=0.5, just=c("left","centre"))}else{
		for(i in 1:ncol(phenotypes)){
		grid.text(names(phenotypes)[i],x=0.05, y=(i/ncol(phenotypes)-0.5/ncol(phenotypes)), just="left", gp=gpar(cex=cex.label*1.5))}
		}
	popViewport()
	
	sample.vp <- viewport(x = unit(0.25, "npc"),
	y = unit(0.7-phenobar.height,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "top"))
	pushViewport(sample.vp)
	for(i in 1:ncol(eset)){
		grid.text(sampleNames(eset)[i],x=(i/ncol(eset)-0.5/ncol(eset)),
		y=0.95, just=c("right","centre"), gp=gpar(cex=cex.label), rot=90)}
	popViewport()
	
	legend.height <- (1.6/ncol(eset))*length(leg.labels)
	
	legend.vp <- viewport(x = unit(0.05, "npc"),
	y = unit(0.7-phenobar.height,  "npc"),
	width = unit(legend.height/(2*length(leg.labels)), "npc"),
	height = unit(min(legend.height,0.8), "npc"),
	just = c("left", "centre"),
	xscale = leg.grob$x.limits,
    yscale = leg.grob$y.limits)
	pushViewport(legend.vp)
	grid.rect()
	do.call("panel.levelplot", trellis.panelArgs(leg.grob, 1))
	popViewport()
	
	legend.labels.vp <- viewport(x = unit((legend.height/(2*length(leg.labels)))+0.05, "npc"),
	y = unit(0.7-phenobar.height,  "npc"),
	width = unit(0.19-(legend.height/(2*length(leg.labels))), "npc"),
	height = unit(min(legend.height,0.8), "npc"),
	just = c("left", "centre"),
	xscale = leg.grob$x.limits,
    yscale = leg.grob$y.limits)
	pushViewport(legend.labels.vp)
	for(i in 1:length(leg.labels)){
		grid.text(leg.labels[i],x=0.1, y=(i/length(leg.labels)-0.5/length(leg.labels)), just="left", gp=gpar(cex=cex.label*1.5))}
	popViewport()	
	
	subtitle.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "top"))
	pushViewport(subtitle.vp)
	grid.text(paste(nrow(eset),"probes - Cluster method = ",atr$method,"with",dist.method, "distance"))
	popViewport()
		
		}
	
	}
	
	if(length(device) != 0 && device != "quartz") dev.off()
	
	
	if(return.atr == T) return(atr)
}

centerGenes <- function(eset, center="median")
#	Centre expression set values by row(probe)
{
	require(Biobase)
	if(center=="mean"){
		means <- (apply(t(exprs(eset)),2,mean,na.rm=T))
		exprs.cent <- scale(t(exprs(eset)), center=means, scale = F)
		}
	if(center=="median"){
		medians <- (apply(t(exprs(eset)),2,median,na.rm=T))
		exprs.cent <- scale(t(exprs(eset)), center=medians, scale = F)
		}
	
	new.exprs <- t(exprs.cent)
	row.names(new.exprs) <- featureNames(eset)
	colnames(new.exprs) <- sampleNames(eset)
	assayDataElement(eset, "exprs") <- new.exprs

	return(eset)
}

centerArrays <- function(eset, center="median", quant=0.5)
#	Centre expression set values by column(array)
{
	require(Biobase)
	medians <- apply(exprs(eset),2,function(x) quantile(x, quant,na.rm=T))
	means <- apply(exprs(eset),2,mean,na.rm=T)
	if(center=="mean"){exprs(eset) <- scale(exprs(eset), center=means, scale = F)}else{
		exprs(eset) <- scale(exprs(eset), center=medians, scale = F)}
		
	if(!is.null(assayData(eset)$smo)){
	if(center=="mean"){new.smo <- scale(assayData(eset)$smo, center=means, scale = F)}else{
		new.smo <- scale(assayData(eset)$smo, center=medians, scale = F)}
		row.names(new.smo) <- featureNames(eset)
		colnames(new.smo) <- sampleNames(eset)
		assayDataElement(eset, "smo") <- new.smo
	}	
		
	return(eset)
}

repProbeMean <- function(eset, fData.id="probeID")
{
	require(Biobase)
	rep.list <- fData(eset)[,fData.id]
	eset.uni <- eset[match(unique(rep.list), rep.list),]
	cat("Combining",nrow(eset.uni),"replicates\n")
	
	for(j in 1:length(assayDataElementNames(eset)))
	{
		AD.j <- assayData(eset)[[assayDataElementNames(eset)[j]]]
		AD.uni.j <- assayData(eset.uni)[[assayDataElementNames(eset.uni)[j]]]
			for (i in 1:nrow(AD.uni.j))
			{
			AD.i <- AD.j[which(is.element(rep.list, unique(rep.list)[i])), ,drop=F]
			if(nrow(AD.i) > 1) AD.uni.j[i,] <- colMeans(AD.i, na.rm=T)
			}
			
		assayDataElement(eset.uni, assayDataElementNames(eset.uni)[j]) <- AD.uni.j
	}
	
	if(!is.null(assayData(eset.uni)$flags)) {
		flags.uni <- assayData(eset.uni)$flags
		flags.uni[flags.uni != 1] <- 0
		assayDataElement(eset.uni, "flags") <- flags.uni
		}
		
	eset.uni	
	
}

removeProbesByFlags <- function(eset, flag.limit)
#	remove probes from MA or cgh object
#	which have more than a specified number of flagged values - flag.limit
{
	require(Biobase)
	cat("Removing probes with",flag.limit,"flags or more\n")
	if(is.null(assayData(eset)$flags)) assayDataElement(eset,"flags") <- is.na(exprs(eset))
	flag.count <- apply(assayData(eset)$flags, 1,sum)
	eset <- eset[which(flag.count < flag.limit),]
	cat("",nrow(eset),"probes remaining\n","Done\n")
	eset
}

removeUnmappedProbes <- function(eset, verbose=T)
{	
	require(Biobase)
	if(is.null(fData(eset)$chrom)|is.null(fData(eset)$start)|is.null(fData(eset)$end)) stop("Expression set featureData must be annotated with \"chrom\", \"start\" and \"end\"")
	if(verbose==T) cat("Removing unmapped probes\n")	
	fData(eset)$chrom[fData(eset)$chrom == "X"] <- 23
	fData(eset)$chrom[fData(eset)$chrom == "Y"] <- 24

	fData(eset)$chrom <- as.integer(fData(eset)$chrom)
	fData(eset)$start <- as.integer(fData(eset)$start)
	fData(eset)$end <- as.integer(fData(eset)$end)
	
	eset <- eset[!is.na(fData(eset)$chrom),]
	eset <- eset[!is.na(fData(eset)$start),]
	eset <- eset[!is.na(fData(eset)$end),]
	eset <- eset[which(is.element(fData(eset)$chrom,c(1:24))),]
		
	if(verbose==T) cat("",nrow(eset),"mapped probes remaining\n","Done\n")
	eset	
}


calculateFoldChange <- function(eset, pheno)
#	calculate the fold change for one phenotype in an expression set vs the rest of the expression set
#	returns the fold change as a vector
{
	require(Biobase)
	positive <- eset[,which(!is.na(pheno))]
	eset.control <- eset[,!is.element(sampleNames(eset), sampleNames(positive))]
	fold <- 2^(as.vector(esApply(positive, 1, mean, na.rm=T)) - as.vector(esApply(eset.control, 1, mean, na.rm=T)))
	fold
}
	

toPAM <- function(eset,pheno,eset.id="probeID", eset.geneid="symbol")
{
	require(Biobase)
	geneID <- fData(eset)[,which(fvarLabels(eset) == eset.id)]

	genenames <- fData(eset)[,which(fvarLabels(eset) == eset.geneid)]
	
	pamo <- list(geneid=geneID,genenames=genenames,x=exprs(eset),y=pheno,samplelabels=sampleNames(eset))
   return(pamo)
}


runmedByChroms<- function(M,chrom, k)
{
	running.med <- M
	M.chroms <- split(M[is.finite(M)], chrom[is.finite(M)])
	runmed.chroms <- M.chroms
	for(i in 1:length(M.chroms))
	{
		bw <- min(k, length(M.chroms[[i]]))
		if(bw%%2 == 0) bw <- bw-1
		runmed.chroms[[i]]<- runmed(M.chroms[[i]], k=bw)	}	
	running.med[is.finite(M)] <- unlist(runmed.chroms)
	running.med
}

calculateMAD = function(eset, chrom=fData(eset)$chrom, bw=151)
#	Function to calculate median absolute deviation from running median in aCGH
#	returns eset 
{	
	require(Biobase)
	rmed <- matrix(NA, nrow(eset), ncol(eset))
	for(i in 1:ncol(eset))
	{
		rmed[,i] <- runmedByChroms(exprs(eset)[,i], chrom, k=bw)
	}
	rmed.diffs <- (abs(exprs(eset)-rmed))
	MAD <- 1.4826*apply(rmed.diffs,2,median, na.rm=T)
	pData(eset)$MAD <- round(MAD,3)
	varMetadata(eset)$labelDescription[varLabels(eset) == "MAD"] <- "Genome MAD - Genome-wide Median Absolute Deviation from running median"
	
	eset
	
}


calculateMB <- function(cgh)
#	calculate mid-point MB position for plotting genomes
{
	require(Biobase)
	cgh <- cgh[order(fData(cgh)$chrom, fData(cgh)$start),]
	chroms <- unique(fData(cgh)$chrom)
	chrom.1 <- fData(cgh)[which(fData(cgh)$chrom == chroms[1]),]
	MB <- (chrom.1$start+chrom.1$end)/2000000
	for(i in 2:length(chroms))
	{
		ichrom <- fData(cgh)[which(fData(cgh)$chrom==chroms[i]),]
		start <- MB[length(MB)] +20
		mids <- ((ichrom$start + ichrom$end)/2000000) +start
		MB <- c(MB,mids)
	}
	MB
}

removeOutliersByMAD <- function(eset, MAD=1, bw=3)
{
	require(Biobase)
	cat(paste("Removing probes with a",MAD,"MAD deviation from the running median of",bw,"\n"))
#	commands <- c(cgh$commands,deparse(match.call()))
	if(is.null(pData(eset)$MAD)) eset <- calculateMAD(eset)
	new.exprs <- exprs(eset)
	for(i in 1:ncol(eset))
		{
		rmed <- runmedByChroms(M=exprs(eset)[,i],fData(eset)$chrom, k = bw)
		rmed.diff <- abs(exprs(eset)[,i]-rmed)/exp(abs(exprs(eset)[,i]))
		new.exprs[rmed.diff>pData(eset)$MAD[i]*MAD,i] <- NA
		cat(i,sampleNames(eset)[i],"\n")
		}
	
	assayDataElement(eset, "exprs") <- new.exprs
	cat("Recalculating MADs...")
	eset <- calculateMAD(eset)
#	cgh$commands <- commands	
	cat("Done\n")
	eset
}


imputeCGH <- function(eset, M=exprs(eset), span=2)
{
#	commands <- c(cgh$commands,deparse(match.call()))
	require(Biobase)
	cat("Imputing missing values\n")
	
	for(j in 1:ncol(eset))
	{
		for(i in unique(fData(eset)$chrom))
		{
			iM <- M[(fData(eset)$chrom==i),j]
			M[(fData(eset)$chrom==i),j] <- imputeChrom(iM, span=span)
			cat(i, "\t")
		}
	cat(sampleNames(eset)[j], "\n")	
	}
#	cgh$commands <- commands
	cat("Done\n")
	M
}


imputeChrom <- function(M, span)
#	Impute missing values into cgh object
#	Use the median of the running median
#	of the values either side of the missing value
#	Width of the values to use is set by span
{
	require(Biobase)
	start.mean <- rep(NA,length(M))
	start.mean <- mean(M[!is.na(M)][1:span])
	start.means <- rep(start.mean,span)
		
	end.mean <- rep(NA,length(M))
	end.mean <- mean(M[!is.na(M)][(length(M[!is.na(M)])):(length(M[!is.na(M)])-span)])
	end.means <- rep(end.mean,span)

	rmed <- c(start.means,M,end.means)
	nas <- (1:length(rmed))[is.na(rmed)]
	rmed.out <- rmed
	
	for(j in nas)
	{
	below <- (rmed[1:j])[!is.na(rmed[1:j])]
	below <- below[length(below):(length(below)-span+1)]
	above <- (rmed[j:length(rmed)])[!is.na(rmed[j:length(rmed)])]
	above <- above[1:span]	
	rmed.out[j] <- mean(c(above,below))
	}

	rmed.out <- rmed.out[(span+1):(length(rmed.out)-span)]
	rmed.out
}

cbsCGH <- function(cgh, cbs.alpha=NULL, ...)
{
	require(Biobase)
	require(DNAcopy)
	options(warn=-1)
#	commands <- c(cgh$commands,deparse(match.call()))
	cat("Smoothing data using cbs\n")
	
	MB <- calculateMB(cgh)
	if(sum(is.na(exprs(cgh))) > 0) cgh <- imputeCGH(cgh, span=2)
	
	cbs.id <- paste("cbs", c(1:length(sampleNames(cgh))), sep=".")
	cnaObject <- CNA(exprs(cgh), fData(cgh)$chrom, MB, data.type = "logratio",sampleid = cbs.id)

	DC.smo <- smooth.CNA(cnaObject)
	if(is.null(cbs.alpha)) cbs.alpha=0.01
	DC.seg <- segment(DC.smo, verbose = 2, alpha=cbs.alpha, ...)
	
	cbs <- exprs(cgh)
	for (i in 1:length(cbs.id))
	{
	output.i <- DC.seg$output[which(DC.seg$output$ID == cbs.id[i]),]
	ok <- is.finite(cbs[,i])
	cbs[ok,i] <- rep(output.i$seg.mean,output.i$num.mark)
	}
	
	assayDataElement(cgh, "smo") <- cbs
	
	cat("Done")
#	cgh$commands <- commands
	cgh
}

kernelSmoothCGH <- function(cgh, bandwidth=NULL)
{
	require(Biobase)
	require(aroma.core)

#	commands <- c(cgh$commands,deparse(match.call()))
	cat("Kernel Smoothing data with a bandwidth of",bandwidth, "\n")
	
	cgh <- cgh[order(fData(cgh)$chrom, fData(cgh)$start),]
	MB <- calculateMB(cgh)
	if(sum(is.na(exprs(cgh))) > 0) assayDataElement(cgh, "exprs") <- imputeCGH(cgh, span=2)
	
	kern.exprs <- exprs(cgh)
	chroms <- unique(fData(cgh)$chrom)
	for (i in 1:length(chroms))
	{
		
	kern.exprs[fData(cgh)$chrom == chroms[i],] <- colKernelSmoothing(Y=exprs(cgh)[fData(cgh)$chrom == chroms[i],], x=1:sum(fData(cgh)$chrom == chroms[i]), kernel="gaussian", h=bandwidth, na.rm=T)
	cat("Chromosome\t")
	cat(unique(fData(cgh)$chrom)[i], "\t")
	}
	
	row.names(kern.exprs) <- featureNames(cgh)
	colnames(kern.exprs) <- sampleNames(cgh)
	
	assayDataElement(cgh, "exprs") <- kern.exprs
	
	cat("\nDone")
#	cgh$commands <- commands
	cgh
}


awsCGH <- function(cgh)
{
	require(Biobase)
	require(aws)
#	commands <- c(cgh$commands,deparse(match.call()))
	cat("Smoothing data using aws\n")
	MB <- calculateMB(cgh)
	smooth <- exprs(cgh)
	for(i in 1:ncol(cgh))
		{
		ok <- is.finite(exprs(cgh)[,i])
		M.ok <- exprs(cgh)[ok,i]
		chroms.ok <- fData(cgh)$chrom[ok]
		MB.ok <- MB[ok]

		chroms <- split(M.ok, chroms.ok)
		MB.chroms <- split(MB.ok, chroms.ok)
		
		smo <- chroms
		for (j in 1:length(chroms))
			{
			MB.start <- MB.chroms[[j]][1]
			MB.end <- MB.chroms[[j]][length(MB.chroms[[j]])]
			ismo <- aws(y=chroms[[j]], hmax=length(chroms[[j]]))
			smo[[j]] <- ismo$theta
			cat(j,"\t")
			}
		cat(sampleNames(cgh)[i],"\n")
		
		smoothed.data <- unlist(smo)
		smooth[ok,i] <- smoothed.data
		}
	assayDataElement(cgh, "smo") <- smooth
	cat("Done")
#	cgh$commands <- commands
	cgh
}

thresholdCGH <- function(cgh, gainthresh, lossthresh=NULL, ampthresh=NULL, delthresh=NULL, verbose=T)
{
	require(Biobase)
	if(verbose == T) cat("thresholding aCGH data\n")
	GL <- matrix(0,nrow(cgh),ncol(cgh))
	
	smo <- assayData(cgh)$smo

	GL[smo  > gainthresh] <- 1
	if(is.null(lossthresh))
		{lossthresh <- gainthresh*(-1)}
	GL[smo  < lossthresh] <- -1
	
	if(!is.null(ampthresh))
		{
		if(is.null(delthresh))
			{delthresh <- ampthresh*(-1)}
		GL[smo  > ampthresh] <- 2		
		GL[smo  < delthresh] <- -2
		}
	
	rownames(GL) <- featureNames(cgh)
	colnames(GL) <- sampleNames(cgh)
	thresholds <- list(gainthresh, lossthresh, ampthresh, delthresh)
	names(thresholds) <- c("gainthresh","lossthresh", "ampthresh","delthresh")
	result <- list(GL, thresholds)
	names(result) <- c("GL", "thresholds")
	result
}



callCGHStatesThreshold <- function(cgh, gainthresh, lossthresh=NULL, ampthresh=NULL, delthresh=NULL,contig=3, verbose=T)
#	categorize smoothed aCGH data based upon thresholds
#	currently uses two thresholds
#	gainthresh for "gain/loss" and ampthresh for "amp/del"
#	records gain/loss as 1/-1 and amp/del as 2/-2 in a GL matrix
#	if ampthresh is not specified only gain/loss is recorded
#	currently symmetrical although aCGH data is not!
#	filters GL table using contigStates for contiguous regions >contig
#	returns GL matrix of integral data
#	counts "gains/losses" and "amps/dels" in a GLAD matrix
{
#	commands <- c(cgh@commands,deparse(match.call()))
	require(Biobase)
	smo <- assayData(cgh)$smo
	
	GL.result <- thresholdCGH(cgh=cgh, gainthresh=gainthresh, lossthresh=lossthresh, ampthresh=ampthresh, delthresh=delthresh, verbose=verbose)
	GL.raw <- GL.result$GL
	cgh@thresholds <- GL.result$thresholds
		
	GL <- GL.raw
	for (i in 1:ncol(GL.raw))
	{
		GLi <- as.vector(GL.raw[,i])
		GL.chroms <- split(GLi, fData(cgh)$chrom)
		for (j in 1:length(GL.chroms))
			{
			GL.chroms[[j]] <- contigStates(GL.chroms[[j]], contig=contig)
			cat(j,"\t")
			}
		cat(sampleNames(cgh)[i],"\n")
		GL[,i] <- unlist(GL.chroms)
	}
	smo[GL!=GL.raw] <- NA
	
	if(verbose == T) cat("Filtering contiguous states...")
	assayDataElement(cgh, "smo") <- smo
	if(verbose == T) cat("Done\n")
	assayDataElement(cgh, "smo") <- imputeCGH(cgh, M=assayData(cgh)$smo)
	
#	cgh$contig <- (contig)
#	names(cgh$contig) <- "contig"
	assayDataElement(cgh, "GL") <- GL

#	cgh@commands <- commands
	cgh
	
}


thresholdMADCGH <- function(cgh, gainMADS, lossMADS=NULL, ampMADS=NULL, delMADS=NULL, verbose=T)
{
	require(Biobase)
	if(verbose==T) cat("thresholding aCGH data\n")
	
	if(is.null(pData(cgh)$MAD)) stop("Please calculate genome MADs first")
	GL <- matrix(0, nrow(cgh), ncol(cgh))
	
	smo <- assayData(cgh)$smo

	for(i in 1:ncol(cgh))
	{
		MAD <- pData(cgh)$MAD[i]
		
		gainthresh <- MAD*gainMADS 
		GL[smo[,i]  >= gainthresh,i] <- 1

		if(is.null(lossMADS)) {lossMADS <- gainMADS}
		lossthresh <- (MAD*lossMADS*-1)
		GL[smo[,i]  <= lossthresh,i] <- -1
	
		if(!is.null(ampMADS)) {ampthresh <- MAD*ampMADS}
		GL[smo[,i]  >= ampthresh,i] <- 2

		if(is.null(delMADS)) {delMADS <- ampMADS}
		delthresh <- (MAD*delMADS*-1)
		GL[smo[,i]  <= delthresh,i] <- -2
		
		cat(i, sampleNames(cgh)[i], "\n")
		
	}

	rownames(GL) <- featureNames(cgh)
	colnames(GL) <- sampleNames(cgh)
	MADS <- list(gainMADS, lossMADS, ampMADS, delMADS)
	names(MADS) <- c("gainMADS","lossMADS", "ampMADS","delMADS")
	result <- list(GL, MADS)
	names(result) <- c("GL", "MADS")
	result
}


callCGHStatesMAD <- function(cgh, gainMADS, lossMADS=NULL, ampMADS=NULL, delMADS=NULL, contig=3)
#	categorize smoothed aCGH data based upon thresholds
#	currently uses two thresholds
#	gainthresh for "gain/loss" and ampthresh for "amp/del"
#	records gain/loss as 1/-1 and amp/del as 2/-2 in a GL matrix
#	if ampthresh is not specified only gain/loss is recorded
#	currently symmetrical although aCGH data is not!
#	filters GL table using contigStates for contiguous regions >contig
#	returns GL matrix of integral data
#	counts "gains/losses" and "amps/dels" in a GLAD matrix
{
	require(Biobase)
	smo <- assayData(cgh)$smo	
	
	GL.result <- thresholdMADCGH(cgh=cgh, gainMADS=gainMADS, lossMADS=lossMADS, ampMADS=ampMADS, delMADS=delMADS, verbose=F)
	GL.raw <- GL.result$GL
	cgh@MADS <- GL.result$MADS
	
	GL <- GL.raw
	for (i in 1:ncol(GL.raw))
	{
		GLi <- as.vector(GL.raw[,i])
		GL.chroms <- split(GLi, fData(cgh)$chrom)
		for (j in 1:length(GL.chroms))
			{
			GL.chroms[[j]] <- contigStates(GL.chroms[[j]], contig=contig)
			cat(j,"\t")
			}
		cat(sampleNames(cgh)[i],"\n")
		GL[,i] <- unlist(GL.chroms)
	}
	smo[GL!=GL.raw] <- NA
	
	cat("Filtering contiguous states...")
	assayDataElement(cgh, "smo") <- smo
	
	assayDataElement(cgh, "smo") <- imputeCGH(cgh, M=assayData(cgh)$smo)
	
#	cgh$contig <- (contig)
#	names(cgh$contig) <- "contig"
		
	assayDataElement(cgh, "GL") <- GL
	cgh
	
}


contigStates <- function(states, contig=3)
{
#	a function to take integral data as a vector of states (-2,-1,0,1,2)
#	record the positions at which the input states change
#	calculate the lengths of the contiguous regions
#	replace the values in short regions (below contig)
#	currently set to identify maximum amplicons
#	with values from the neighbouring large regions
	require(Biobase)
	len <- length(states)
	start <- states[1]
	end <- states[len]
	index <- 1:len
	plus1 <- c(states[2:len],end)
	min1 <- c(start, states[1:len-1])
	offset <- rbind(index,min1,states,plus1)
	offset <- as.matrix(offset)
	
# 	offset is a matrix of
#	index
#	min1
#	states
#	plus1
# 	contains all input values
#	use offsets to calculate positive and negative lengths
#	record index positions in which to replace values
	
	p.breaks <- offset[,which(offset[3,] != offset[2,])]
	if(length(p.breaks) > 4)
	{
	pbc <- ncol(p.breaks)
	p.lengths <- as.vector(p.breaks[1,2:pbc])-as.vector(p.breaks[1,1:(pbc-1)])
	p.lengths <- c(p.lengths,(len-p.breaks[1,pbc]))
	p.breaks.table <- rbind(p.breaks,p.lengths)	
	n.breaks <- offset[,which(offset[3,] != offset[4,])]
	nbc <- ncol(n.breaks)
	n.lengths <- as.vector(n.breaks[1,2:nbc])-as.vector(n.breaks[1,1:(nbc-1)])
	n.lengths <- c(n.breaks[1,1],n.lengths)
	n.breaks.table <- rbind(n.breaks,n.lengths)
	start.breaks <- as.matrix(c(offset[,1],p.breaks.table[1,1]-1))
	end.breaks <- as.matrix(c(offset[,len],len-p.breaks.table[1,pbc]))
	breaks.table <- cbind(start.breaks,p.breaks.table,n.breaks.table,end.breaks)
	dimnames(breaks.table)[1] <- list(c("index","minus1","states","plus1","lengths"))
	breaks.table <- breaks.table[,order(breaks.table[1,])]

#	breaks.table is a matrix of
#	index
#	min1
#	states
#	plus1
#	lengths
#	for every point where input is different to input+1 and input-1
#	p denoted the plus lengths and n denotes the negative lengths in the input
#	start and end states are added as placeholders
	
#	make a matrix of all the breaks bigger than contig called big.breaks
	ncb <-ncol(breaks.table)
#	breaks just lists the breaks without the start and end states
	breaks <- breaks.table[,2:(ncb-1)]

#	make a matrix of all the breaks bigger than contig called big.breaks
	big.breaks <- breaks[,which(breaks[5,]>=contig)]

#	add placeholders back for start and end
#	this is the matrix in which to look for the correct states - rep.states
	breaks.index <- cbind(breaks.table[,1],big.breaks,breaks.table[,ncb])

#	make a matrix of all the breaks smaller than contig called small.breaks
#	these are the values we'd like to replace from breaks.index
	small.breaks <- as.matrix(breaks[,which(breaks[5,]<contig)])
	
#	loop through the small.breaks table
# 	look up index positions and states
# 	for previous and next positions in small.breaks table
#	and replace the states in the original offset table
#	with the previous or the next state from the big.breaks table
#	whichever you want
	
	if(ncol(small.breaks) >0)
	{
		for(i in 1:ncol(small.breaks))
		{
			current.index <- small.breaks[1,i]

			bigger <- breaks.index[1,which(breaks.index[1,]>current.index)]
			if(length(bigger)==0){next.index <-length(states)}
			if(length(bigger)>0){next.index <-min(bigger)}

			next.state <- states[next.index]
			
			smaller <- breaks.index[1,which(breaks.index[1,]<current.index)]
			if(length(smaller) == 0){prev.index <- 1}
			if(length(smaller) > 0){prev.index <-max(smaller)}

			prev.state <- states[prev.index]
			
			av.state <- sum(next.state, prev.state)
			
#	pick which of the neighbouring states to use - biggest or smallest
#	if you want minimum amplicons set these to max and then min
#	if you want maximum amplicons set these to min and then max

			if(av.state<0){rep.state <- max(next.state, prev.state)}
			if(av.state>0){rep.state <- min(next.state, prev.state)}
			if(av.state==0){rep.state <- 0}
			states[(prev.index+1):current.index] <- rep.state
		}	
	}
	}
	states
}
	
RLEStates <- function(states, contig=3, amps="minimum")
{
#	a function to take integral data as a vector of states (-2,-1,0,1,2)
#	record the positions at which the input states change
#	calculate the lengths of the contiguous regions
#	replace the values in short regions (below contig)
#	currently set to identify maximum amplicons
#	with values from the neighbouring large regions
#	Use run length encoding to identify coniguous regions

	require(Biobase)
	rle <- rle(states)
	index <- cumsum(rle$lengths)
	values <- rle$values
	lengths <- rle$lengths
	rle <- rbind(index,values,lengths)
		
#	rle.result is a matrix of
#	index
#	state up to that point
#	lengths

#	make a matrix of all the breaks bigger than contig called big.breaks
	big.breaks <- rle[,which(rle[3,]>=contig)]

#	add placeholders back for start and end
#	this is the matrix in which to look for the correct states - rep.states
	breaks.index <- cbind(c(1,states[1],1),big.breaks)

#	make a matrix of all the breaks smaller than contig called small.breaks
#	these are the values we'd like to replace from breaks.index
	small.breaks <- rle[,which(rle[3,]<contig)]

#	check there is something to replace!
	
#	loop through the small.breaks table
# 	look up index positions and states
# 	for previous and next positions in small.breaks table
#	and replace the states in the original offset table
#	with the previous or the next state from the big.breaks table
#	whichever you want
	
	if(ncol(small.breaks) >0)
	{
		for(i in 1:ncol(small.breaks))
		{
			current.index <- small.breaks[1,i]
			current.length <- small.breaks[2,i]

			bigger <- breaks.index[1,which(breaks.index[1,]>current.index)]
			if(length(bigger)==0){next.index <-length(states)}
			if(length(bigger)>0){next.index <-min(bigger)}

			next.state <- states[next.index]
			
			smaller <- breaks.index[1,which(breaks.index[1,]<current.index)]
			if(length(smaller) == 0){prev.index <- 1}
			if(length(smaller) > 0){prev.index <-max(smaller)}

			prev.state <- as.numeric(states[prev.index])
			
			av.state <- sum(next.state, prev.state)
			
#	pick which of the neighbouring states to use - biggest or smallest
#	if you want minimum amplicons set these to max and then min
#	if you want maximum amplicons set these to min and then max
			if(amps == "minimum")
			{
			if(av.state<0){rep.state <- max(next.state, prev.state)}
			if(av.state>0){rep.state <- min(next.state, prev.state)}
			if(av.state==0){rep.state <- 0}
			}
			if(amps == "maximum")
			{
			if(av.state>0){rep.state <- max(next.state, prev.state)}
			if(av.state<0){rep.state <- min(next.state, prev.state)}
			if(av.state==0){rep.state <- 0}
			}
			states[(prev.index+1):current.index] <- rep.state
		}	
	}
	states
}


addCytobands <- function(cgh, cytoband.table="cytobands.txt", verbose=T)
{
	cytobands.table <- read.delim(cytoband.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	
	cytobands <- rep(NA, nrow(cgh))
	for(i in 1:nrow(cgh))
	{
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == fData(cgh)$chrom[i]),,drop=F]
		ibands.max <- ibands.chrom[which(ibands.chrom$end > fData(cgh)$start[i]),,drop=F]
		ibands <- ibands.max[which(ibands.max$start <= fData(cgh)$end[i]),,drop=F]

		if(nrow(ibands)>1) {cytobands[i] <- paste(ibands$cytoband[1],ibands$cytoband[nrow(ibands)],sep="-")}
		if(nrow(ibands) == 1) {cytobands[i] <- ibands$cytoband}
		if(verbose==T & i %% 100 == 0) cat(i," ")
	}
		fData(cgh)$cytoband <- cytobands
		cgh	
}

armChangesCGH <- function(cgh, WAG.thresh, PG.contig=3, chrom.arms=NULL, project="aCGH")
{

require(Biobase)
if(is.null(chrom.arms)){
if(is.null(fData(cgh)$cytoband)) stop("cgh object must be annotated with cytoband")
chrom.arms <- gsub("-", ".",fData(cgh)$cytoband)
chrom.arms <- gsub("q\\w*\\.*\\w*", "q", chrom.arms, perl=T)
chrom.arms <- gsub("q.*", "q", chrom.arms, perl=T)
chrom.arms <- gsub("p\\w*\\.*\\w*", "p", chrom.arms, perl=T)
chrom.arms <- gsub("p.*", "p", chrom.arms, perl=T)
chrom.arms <- paste(fData(cgh)$chrom, chrom.arms, sep="")
}

chr.arms <- unique(chrom.arms)
cat("Chromosome arms\n", chr.arms, "\n")

chr.arm.lengths <- as.vector(table(chrom.arms))

chr.arm.L <- chr.arm.G  <- chr.arm.prop.gain <- chr.arm.prop.loss <- matrix("NC", ncol(cgh), length(chr.arms))

cat("Counting losses and gains in chromosome arms\n")
for(i in 1:length(chr.arms)){

	GL.i <- assayData(cgh)$GL[chrom.arms == chr.arms[i],]
	GL.i.gain <- GL.i >= 1
	
	G.prop <- chr.arm.prop.gain[,i] <- signif((apply(GL.i.gain,2,sum)/nrow(GL.i))*100,3)
	chr.arm.G[which(apply(GL.i.gain,2,sum) > PG.contig),i] <- "PG"
	chr.arm.G[which(G.prop > WAG.thresh),i] <- "WAG"

	GL.i.loss <- GL.i <= -1
	L.prop <- chr.arm.prop.loss[,i] <- signif((apply(GL.i.loss,2,sum)/nrow(GL.i))*100,3)

	chr.arm.L[which(apply(GL.i.gain,2,sum) > PG.contig),i] <- "PL"
	chr.arm.L[which(L.prop > WAG.thresh),i] <- "WAL"

	cat(i, "\t")
}

arm.table.prop.loss <- cbind.data.frame(pData(cgh), chr.arm.prop.loss)
arm.table.prop.gain <- cbind.data.frame(pData(cgh), chr.arm.prop.gain)

names(arm.table.prop.loss) <- c(names(pData(cgh)), paste(chr.arms, "%.loss", sep="."))
names(arm.table.prop.gain) <- c(names(pData(cgh)), paste(chr.arms, "%.gain", sep="."))

write.table(arm.table.prop.gain, paste(project,"prop.gained.xls", sep="."), sep="\t", row.names=F)
write.table(arm.table.prop.loss, paste(project,"prop.lost.xls", sep="."), sep="\t", row.names=F)

arm.table.L <- cbind.data.frame(pData(cgh), chr.arm.L)
arm.table.G <- cbind.data.frame(pData(cgh), chr.arm.G)

names(arm.table.L) <- c(names(pData(cgh)), paste(chr.arms, "loss", sep="."))
names(arm.table.G) <- c(names(pData(cgh)), paste(chr.arms, "gain", sep="."))

write.table(arm.table.G, paste(project,"arm.calls.gains.xls", sep="."), sep="\t", row.names=F)
write.table(arm.table.L, paste(project,"arm.calls.losses.xls", sep="."), sep="\t", row.names=F)

chr.arm.prop.loss.df <- data.frame(chr.arm.prop.loss) 
names(chr.arm.prop.loss.df) <- paste("chr",chr.arms,sep=".")
chr.arm.prop.gain.df <- data.frame(chr.arm.prop.gain) 
names(chr.arm.prop.gain.df) <- paste("chr",chr.arms,sep=".")

chr.arm.L.df <- data.frame(chr.arm.L) 
names(chr.arm.L.df) <- paste("chr",chr.arms,sep=".")
chr.arm.G.df <- data.frame(chr.arm.G) 
names(chr.arm.G.df) <- paste("chr",chr.arms,sep=".")

return(list(loss.props=chr.arm.prop.loss.df, loss.calls=chr.arm.L.df, gain.props=chr.arm.prop.gain.df, gain.calls=chr.arm.G.df))

}

HicksScore <- function(cgh, project)
#	a script to take integral data  from a cgh object as a vector of states (-2,-1,0,1,2)
#	record the positions at which the input states change
#	calculate the lengths of the contiguous regions
# 	Calculate Hicks F score for genomic rearrangements
{
#	commands <- c(cgh$commands,deparse(match.call()))
	require(Biobase)
	cat("Calculating F scores...")

	F.score.gain <- numeric(ncol(cgh))
	
	gain.table <- assayData(cgh)$GL
	
	for(i in 1:ncol(cgh))
		{
		chroms <- split(gain.table[,i], fData(cgh)$chrom)
		breaks <- rle(chroms[[1]])
		lengths <- breaks$lengths
		lr <- lengths[-1]
		ll <- lengths[-length(lengths)]
		F.table <- cbind(ll, lr)
		F.score <- 2/apply(F.table, 1,sum)
		F.table <- cbind(F.table, F.score)

		for(j in 2:22)
			{
			breaks <- rle(chroms[[j]])
			lengths <- breaks$lengths
			lr <- lengths[-1]
			ll <- lengths[-length(lengths)]
			F.table.j <- cbind(ll, lr)
			F.score <- 2/apply(F.table.j, 1,sum)
			F.table.j <- cbind(F.table.j, F.score)
			F.table <- rbind(F.table, F.table.j)
			}
	
		F.score.gain[i] <- sum(F.table[,3])
		}
	
	F.score.amp <- numeric(ncol(cgh))

	amp.table <- matrix(0, nrow(cgh), ncol(cgh))
	amp.table[assayData(cgh)$GL == 2] <- 2
	amp.table[assayData(cgh)$GL == -2] <- -2
	
	for(i in 1:ncol(assayData(cgh)$GL))
		{
		chroms <- split(amp.table [,i], fData(cgh)$chrom)
		breaks <- rle(chroms[[1]])
		lengths <- breaks$lengths
		lr <- lengths[-1]
		ll <- lengths[-length(lengths)]
		F.table <- cbind(ll, lr)
		F.score <- 2/apply(F.table, 1,sum)
		F.table <- cbind(F.table, F.score)

		for(j in 2:22)
			{
			breaks <- rle(chroms[[j]])
			lengths <- breaks$lengths
			lr <- lengths[-1]
			ll <- lengths[-length(lengths)]
			F.table.j <- cbind(ll, lr)
			F.score <- 2/apply(F.table.j, 1,sum)
			F.table.j <- cbind(F.table.j, F.score)
			F.table <- rbind(F.table, F.table.j)
			}
	
		F.score.amp[i] <- sum(F.table[,3])
		}
	
	pData(cgh)$F.score.amp <- round(F.score.amp,2)
	pData(cgh)$F.score.gain <- round(F.score.gain,2)
	Hicks.result <- data.frame(sampleNames=sampleNames(cgh), F.score.amp=pData(cgh)$F.score.amp, F.score.gain=pData(cgh)$F.score.gain)
	write.table(Hicks.result, file=paste(project,"Hicks.scores.xls", sep="."), sep="\t",row.names=F)
#	cgh$commands <- commands
	cat("Done\n")
	cgh
}



countEventsCGH <- function(cgh, project, return.CGH=F){
	require(Biobase)
	total.counts <- gains <- losses <- amps <- dels <- total.prop <- gain.prop <- loss.prop <- amp.prop <- del.prop <- rep(NA, ncol(cgh))

	rle.states <- apply(assayData(cgh)$GL, 2, rle)
	for(i in 1:ncol(assayData(cgh)$GL)){
		if(length(which(rle.states[[i]]$values != 0)) > 0){
		total.counts[i] <- length(which(rle.states[[i]]$values != 0))
		gains[i] <- length(which(rle.states[[i]]$values == 1))
		losses[i] <- length(which(rle.states[[i]]$values == -1))
		amps[i] <- length(which(rle.states[[i]]$values == 2))
		dels[i] <- length(which(rle.states[[i]]$values == -2))
		
		total.prop[i] <- signif(sum(assayData(cgh)$GL[,i] != 0)/nrow(cgh),3)
		gain.prop[i] <- signif(sum(assayData(cgh)$GL[,i] >= 1)/nrow(cgh),3)
		loss.prop[i] <- signif(sum(assayData(cgh)$GL[,i] <= -1)/nrow(cgh),3)
		amp.prop[i] <- signif(sum(assayData(cgh)$GL[,i] >= 2)/nrow(cgh),3)
		del.prop[i] <- signif(sum(assayData(cgh)$GL[,i] <= -2)/nrow(cgh),3)
		
		}}
	
	aCGH.counts <- data.frame(total.counts, gains, losses, amps, dels, total.prop, gain.prop, loss.prop, amp.prop, del.prop)
	aCGH.counts.pheno <- data.frame(sampleNames=sampleNames(cgh), aCGH.counts)
	
	if(is.null(pData(cgh))){
		pData(cgh) <- data.frame(sampleNames=sampleNames(cgh), aCGH.counts)
		}else{
			pData(cgh) <- data.frame(pData(cgh), aCGH.counts)
			}
	
	
		write.table(aCGH.counts.pheno, paste(project, "aCGH.events.xls", sep="."), sep="\t", row.names=F, na="")
		
	if(identical(return.CGH,T)){
		return(cgh)
		}
		
}

genomePlot <- function(cgh, main=NULL, chroms=NULL, case, thresh=NULL, yAxis=NULL, colourGL=T, plot.smo=T, point.cex=0.5, point.pch=16)
{
	require(Biobase)
	if(!is.numeric(case)) {
		case <- match(case, sampleNames(cgh))
		if(is.na(case)) stop("SampleName is not found in the ExpressionSet")
		}
	if(!is.element(case,1:ncol(cgh))) stop("Case is not found in the ExpressionSet")	
	if(!is.null(chroms))
	{
	cgh <- cgh[which(is.element(fData(cgh)$chrom,chroms)),]
	}else{chroms <- unique(fData(cgh)$chrom)}

	if(is.null(main))
	{
	main<- sampleNames(cgh)[case]
	}
	ok <- is.finite(exprs(cgh)[,case])
	
	if(length(chroms)==1){MB <- fData(cgh)$MB <- ((fData(cgh)$start +fData(cgh)$end)/2000000)}else{
		MB <- fData(cgh)$MB <- calculateMB(cgh)}
	
	chrom.maps <- split(fData(cgh)[,],fData(cgh)$chrom)
	ends <- sapply(chrom.maps,function(x) max(x$MB, na.rm=T))
	ends <- ends+10
	begs <- sapply(chrom.maps,function(x) min(x$MB, na.rm=T))
	begs <- begs-10
	mids <- (begs + ends)/2
	
	M <- exprs(cgh)[,case]
	GL <- assayData(cgh)$GL[,case]
	
	if(is.null(yAxis))
		{
		yAxis <- c(min(round(min(M, na.rm=T))-0.5,-1), max(round(max(M, na.rm=T))+0.5,1))
		}
		
	chlabs <- as.character(unique(fData(cgh)$chrom))
	chlabs[which(chlabs==23)] <- "X"
	chlabs[which(chlabs==24)] <- "Y"
	yscale <- seq(-30,30,1)
	yscale <- yscale[yscale >= yAxis[1] & yscale <= yAxis[2]]
		
	gains.MB <- losses.MB <- amps.MB <- dels.MB <- MB[ok]
	gains <- losses <- amps <- dels <- M[ok]
	
	gains.MB[GL<=0] <- NA
	gains[GL<=0] <- NA
	losses.MB[GL>=0] <- NA
	losses[GL>=0] <- NA
		
	amps.MB[GL<=1] <- NA
	amps[GL<=1] <- NA
	dels.MB[GL>= -1] <- NA
	dels[GL>= -1] <- NA
	
	ylab <- "log ratio"
	
	if(length(chroms)==1) {
		xlab=""
#		xlab <- paste("Chromosome", chroms, "- MBp")
		plot(MB,M,ylim=yAxis,axes=T,xlab=xlab,					main=main,ylab=ylab, pch=16, cex=0.5, las=1, type="n")
		abline(v=seq(from=min(MB),to=max(MB), by=10), col="grey", lty=8)
		abline(h=seq(from=min(yscale),to=max(yscale), by=1), col="grey")
		abline(h=seq(from=0.5,to=max(yscale), by=1),lty=8, col="grey")
		abline(h=seq(from=-0.5,to=min(yscale), by=-1),lty=8, col="grey")
		abline(h=0, col="grey")
		if(!is.null(thresh)){abline(h=c(-thresh,thresh),lty=8, lwd=2, col="grey")}
		points(MB,M, pch=point.pch, cex=point.cex)
		}else{
			xlab <- ""
		#	xlab <- "Chromosome"
			plot(MB,M,ylim=yAxis,axes=F,xlab=xlab,					main=main,ylab=ylab,cex=1, type="n")
			abline(v=(begs)[-1], col="grey")
			abline(h=seq(from=min(yscale),to=max(yscale), by=1), col="grey")
			abline(h=seq(from=0.5,to=max(yscale), by=1),lty=8, col="grey")
			abline(h=seq(from=-0.5,to=min(yscale), by=-1),lty=8, col="grey")
			abline(h=0, col="grey")
			if(!is.null(thresh)){abline(h=c(-thresh,thresh),lty=8, lwd=2, col="grey")}
			points(MB,M, pch=point.pch, cex=point.cex)
			axis(side=1, at=mids,labels=chlabs,tick=F)
			axis(side=2,at=yscale,labels=yscale, tick=T, adj=1, cex=1, las=1)
			box(lwd=1)
			}
	
		
	if(!is.null(assayData(cgh)$GL)){
		if(colourGL=="GLAD"|colourGL==T)
		{
		points(gains.MB, gains, pch=point.pch, cex=point.cex, col="green4")
		points(losses.MB, losses, pch=point.pch, cex=point.cex, col="red4")
		points(amps.MB, amps, pch=point.pch, cex=point.cex, col="green")
		points(dels.MB, dels, pch=point.pch, cex=point.cex, col="red")
		}
		
		if(colourGL=="GL")
		{
		points(gains.MB, gains, pch=point.pch, cex=point.cex, col="green4")
		points(losses.MB, losses, pch=point.pch, cex=point.cex, col="red4")
		points(amps.MB, amps, pch=point.pch, cex=point.cex, col="green4")
		points(dels.MB, dels, pch=point.pch, cex=point.cex, col="red4")
		}
		
		if(colourGL=="AD")
		{
		points(amps.MB, amps, pch=point.pch, cex=point.cex, col="green")
		points(dels.MB, dels, pch=point.pch, cex=point.cex, col="red")
		}
	}
	
	if(plot.smo == TRUE && !is.null(assayData(cgh)$smo))
		{
			smo.split <- split(assayData(cgh)$smo[,case],fData(cgh)$chrom)
			smo.split <- c(smo.split,NA)
			MB.split <- split(MB,fData(cgh)$chrom)
			MB.split <- c(MB.split,NA)
			for(i in 1:(length(smo.split))) lines(MB.split[[i]],smo.split[[i]],col="purple", lwd=2)
		}
				
		
}

multiGenomePlots <- function(cgh, project, chroms=NULL, yAxis=NULL, device="quartz", thresh=NULL, colourGL=T, plot.smo=T, point.cex=0.5, point.pch=16)
{
	require(Biobase)
	if(!is.element("genome.plots", list.files())) dir.create("genome.plots")
	for(i in 1:ncol(cgh))
	{
		cat(sampleNames(cgh)[i], "\n")
		if(device == "quartz")
			{
			quartz(title= sampleNames(cgh)[i],width=6, height=4, dpi=75)
			genomePlot(cgh=cgh, case=i, thresh=thresh, chroms=chroms, yAxis=yAxis, colourGL=colourGL, plot.smo=plot.smo, point.pch=point.pch, point.cex=point.cex)
			}
		
		if(device == "PDF")
			{
			plot.file <- paste(project,i,sampleNames(cgh)[i], "pdf", sep=".")
			pdf(file=paste("genome.plots/",plot.file,sep=""),width=9, height=6)
			genomePlot(cgh=cgh, case=i, thresh=thresh, chroms=chroms, yAxis=yAxis, colourGL=colourGL, plot.smo=plot.smo, point.pch=point.pch, point.cex=point.cex)
			dev.off()
			}
		
		if(device == "PNG")
			{
			plot.file <- paste(project,i,sampleNames(cgh)[i], "png", sep=".")
			png(file=paste("genome.plots/",plot.file,sep=""),width=1200, height=800, pointsize=12)
			genomePlot(cgh=cgh, case=i, thresh=thresh, chroms=chroms, yAxis=yAxis, colourGL=colourGL, plot.smo=plot.smo, point.pch=point.pch, point.cex=point.cex)
			dev.off()
			}
		
		if(device == "JPEG")
			{
			plot.file <- paste(project,i,sampleNames(cgh)[i], "jpg", sep=".")
			jpeg(file=paste("genome.plots/",plot.file,sep=""),width=1200, height=800, pointsize=12, quality=100)
			genomePlot(cgh=cgh, case=i, thresh=thresh, chroms=chroms, yAxis=yAxis, colourGL=colourGL, plot.smo=plot.smo, point.pch=point.pch, point.cex=point.cex)
			dev.off()
			}
			
			if(device == "PS")
			{
			plot.file <- paste(project,i,sampleNames(cgh)[i], "ps", sep=".")
			postscript(file=paste("genome.plots/",plot.file,sep=""),width=9, height=6, pointsize=12)
			genomePlot(cgh=cgh, case=i, thresh=thresh, chroms=chroms, yAxis=yAxis, colourGL=colourGL,plot.smo=plot.smo,  point.pch=point.pch, point.cex=point.cex)
			dev.off()
			}
	}

}

genomePlotGrid <- function(cgh, main=NULL, smooth=NULL, bw=15, chroms=NULL, case, thresh=NULL, yAxis=NULL, colourGL=T)
{
	require(Biobase)
	require(grid)
	if(!is.numeric(case)) {
		case <- match(case, sampleNames(cgh))
		if(is.na(case)) stop("SampleName is not found in the ExpressionSet")
		}
	if(!is.element(case,1:ncol(cgh))) stop("Case is not found in the ExpressionSet")	
	if(!is.null(chroms))
	{
	cgh <- cgh[which(is.element(fData(cgh)$chrom,chroms)),]
	}else{chroms <- unique(fData(cgh)$chrom)}

	if(is.null(main))
	{
	main<- sampleNames(cgh)[case]
	}
	ok <- is.finite(exprs(cgh)[,case])
	
	if(length(chroms)==1){MB <- fData(cgh)$MB <- ((fData(cgh)$start +fData(cgh)$end)/2000000)
		
		}else{
		MB <- fData(cgh)$MB <- calculateMB(cgh)
		
		chrom.maps <- split(fData(cgh)[,],fData(cgh)$chrom)
		ends <- sapply(chrom.maps,function(x) max(x$MB, na.rm=T))
		ends <- ends+10
		begs <- sapply(chrom.maps,function(x) min(x$MB, na.rm=T))
		begs <- begs-10
		mids <- (begs + ends)/2
	
		mid.coords <- mids/max(MB)
		end.coords <- ends/max(MB)
		beg.coords <- begs/max(MB)
		break.coords <- (ends[1:(length(ends)-1)] + begs[2:length(begs)])/2/max(MB)
	
		if(is.element("smo", assayDataElementNames(cgh))) smo.split <- split(assayData(cgh)$smo[,case],fData(cgh)$chrom)
		
		for(i in 1:length(smo.split)) smo.split[[i]] <- c(smo.split[[i]],NA)
		smo.split <- as.vector(unlist(smo.split))

		MB.split <- split(MB,fData(cgh)$chrom)
		for(i in 1:length(MB.split)) MB.split[[i]] <- c(MB.split[[i]],NA)
		MB.split <- as.vector(unlist(MB.split))	
		
		}
	
	M <- exprs(cgh)[,case]
	GL <- assayData(cgh)$GL[,case]
	
	if(is.null(yAxis))
		{
		yAxis <- c(min(range(M)[1]-0.2,-1), max(range(M)[2]+0.2,1))
		}
		
	chlabs <- as.character(unique(fData(cgh)$chrom))
	chlabs[which(chlabs==23)] <- "X"
	chlabs[which(chlabs==24)] <- "Y"
	
	
	yscale.min <- seq(-30.5,30.5,1)
	yscale.min <- yscale.min[yscale.min > yAxis[1] & yscale.min < yAxis[2]]
	yscale.max <- seq(-30,30,1)
	yscale.max <- yscale.max[yscale.max > yAxis[1] & yscale.max < yAxis[2]]
	
	gains.MB <- losses.MB <- amps.MB <- dels.MB <- MB[ok]
	gains <- losses <- amps <- dels <- M[ok]
	
	gains.MB[GL<=0] <- NA
	gains[GL<=0] <- NA
	losses.MB[GL>=0] <- NA
	losses[GL>=0] <- NA
		
	amps.MB[GL<=1] <- NA
	amps[GL<=1] <- NA
	dels.MB[GL>= -1] <- NA
	dels[GL>= -1] <- NA
	
	ylab <- "log ratio"
	
	axis.vp <- viewport(x = unit(0.17, "npc"),
	y = unit(0.2,  "npc"),
	width = unit(0.77, "npc"),
	height = unit(0.7, "npc"),
	just = c("left", "bottom"),
	xscale = range(MB),
	yscale = yAxis)
	
	
	grid.newpage()
	pushViewport(axis.vp)
	grid.rect(height=unit(1, "npc"))
	grid.text("log ratio", x = unit(-0.15, "npc"), gp = gpar(fontsize = 14),rot = 90) 
	grid.yaxis()
	
	grid.lines(range(MB),0, default.units="native", gp=gpar(col="grey"))	
	for(i in 1:length(yscale.min)){
		grid.lines(range(MB),yscale.min[i], gp=gpar(lty=8, col="grey"), default.units="native")
		}
	for(i in 1:length(yscale.max)){
		grid.lines(range(MB),yscale.max[i], gp=gpar(col="grey"), default.units="native")
		}

	
	popViewport()
		
	xaxis.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.2, "npc"),
	just = c("left", "bottom"),
	xscale = range(MB),
	yscale = yAxis, clip="off")
	
	pushViewport(xaxis.vp)
	
	if(length(chroms)==1){grid.text(paste("Chromosome", chroms), y = unit(0.4, "npc"), gp = gpar(cex=1))}else{grid.text("Chromosome", y = unit(0.4, "npc"), gp = gpar(cex=1))}
	
	popViewport()
	
	data.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.2,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.7, "npc"),
	just = c("left", "bottom"),
	xscale = range(MB),
	yscale = yAxis, clip="off")
	
	pushViewport(data.vp)
		
	if(length(chroms) ==1) {
		grid.xaxis()	
		j <- 20
	while(j < max(MB)){
		grid.lines(j,yAxis, gp=gpar(lty=8, col="grey"), default.units="native")
		j <- j+20
		}
		}else{
		
		for(i in 1:length(break.coords))
			{
			grid.lines(break.coords[i],c(0,1),gp=gpar(col="black", lwd=1), default.units="npc")
				for(i in 1:length(chroms))
				{
				grid.text(chlabs[i], x=unit(mid.coords[i],"npc"),y = unit(-0.05, "npc"), just=c("center", "bottom"),gp = gpar(cex=1))}
				}
			}
		
	grid.text(main, y = unit(1, "npc") + unit(1, "lines"), gp = gpar(fontsize = 16))
	
	grid.points(MB, M, pch=16, gp=gpar(cex=0.5))
	
	if(!is.null(assayData(cgh)$GL)){
		if(colourGL=="GLAD"|colourGL==T)
		{
		grid.points(gains.MB, gains, pch=16, gp=gpar(cex=0.5, col="green4"))
		grid.points(losses.MB, losses, pch=16, gp=gpar(cex=0.5, col="red4"))
		grid.points(amps.MB, amps, pch=16, gp=gpar(cex=0.5, col="green"))
		grid.points(dels.MB, dels, pch=16, gp=gpar(cex=0.5, col="red"))
		}
		
		if(colourGL=="GL")
		{
		grid.points(gains.MB, gains, pch=16, gp=gpar(cex=0.5, col="green4"))
		grid.points(losses.MB, losses, pch=16, gp=gpar(cex=0.5, col="red4"))
		grid.points(amps.MB, amps, pch=16, gp=gpar(cex=0.5, col="green4"))
		grid.points(dels.MB, dels, pch=16, gp=gpar(cex=0.5, col="red4"))
		}
		
		if(colourGL=="AD")
		{
		grid.points(amps.MB, amps, pch=16, gp=gpar(cex=0.5, col="green"))
		grid.points(dels.MB, dels, pch=16, gp=gpar(cex=0.5, col="red"))
		}
	}
	
	if(length(chroms)==1){
		if(is.element("smo", assayDataElementNames(cgh))) grid.lines(MB,assayData(cgh)$smo[,case],gp=gpar(col="blue", lwd=2), default.units="native")
	}else{
		if(is.element("smo", assayDataElementNames(cgh))) grid.lines(MB.split,smo.split,gp=gpar(col="blue", lwd=2), default.units="native")
		
		}
	
	popViewport()
		
}

genomePlotHTML <- function(cgh, project, thresh=NULL, yAxis=NULL, colourGL=T, probe.id="bac.id", MB.window=0.1)
{
	require(grid)
	require(Biobase)

	project.folder <- gsub(" ", ".", project)
	
	if(!is.element(project.folder, list.files())) dir.create(project.folder)
	
	html.index <- paste(project.folder, paste(project, "index.html", sep="."), sep="/")
	cat("<body>\n", file=html.index, append=F)
	cat("<p> Interactive genomes - ",project,"project </p>\n", file=html.index, append=T)


	for (i in 1:ncol(cgh))
		{
			cat(paste("<p><a href=",
			paste(paste(paste(sampleNames(cgh)[i],sampleNames(cgh)[i],sep="\\"),"html",sep="."),sampleNames(cgh)[i],sep=">")
			,"</p>", sep=""),"\n", file=html.index, append=T)
		}
	cat("<body>", file=html.index, append=T)
	cat("<html>", file=html.index, append=T)

	for(case in 1:ncol(cgh)){

	case.folder <- paste(project.folder, sampleNames(cgh)[case], sep="/")
	if(!is.element(sampleNames(cgh)[case], list.files(project.folder))) dir.create(case.folder)
	
	MB <- fData(cgh)$MB <- calculateMB(cgh)
	chroms <- unique(fData(cgh)$chrom)
	
	chrom.maps <- split(fData(cgh)[,],fData(cgh)$chrom)
	ends <- sapply(chrom.maps,function(x) max(x$MB, na.rm=T))
	ends <- ends+10
	begs <- sapply(chrom.maps,function(x) min(x$MB, na.rm=T))
	begs <- begs-10
	mids <- (begs + ends)/2
	
		mid.coords <- mids/max(MB)
		end.coords <- ends/max(MB)
		beg.coords <- begs/max(MB)
		break.coords <- c(0,(ends[1:(length(ends)-1)] + begs[2:length(begs)])/2/max(MB),1)

	image.width <- 1200
	image.height <- 600
	
	genome.image.file <- paste(sampleNames(cgh)[case], "png", sep=".")
	
	genome.image <- paste(case.folder,genome.image.file, sep="/")
	genome.page <- paste(case.folder,paste(sampleNames(cgh)[case], "html", sep="."), sep="/")
	genome.map <- paste(case.folder,paste(sampleNames(cgh)[case], "map", sep="."), sep="/")
	
	png(genome.image, width=image.width, height=image.height)
	genomePlotGrid(cgh, case=case, colourGL=colourGL)
	dev.off()
	
	chlabs <- as.character(unique(fData(cgh)$chrom))
		chlabs[which(chlabs==23)] <- "X"
		chlabs[which(chlabs==24)] <- "Y"
		
		
	chroms.alt <- paste("alt= chr",chlabs," title= chr",chlabs, " />", sep="")
	chrom.links <- paste("href=chr",chlabs,"html", sep=".")
		
	cat("<body>", file=genome.page, append=F)
	cat("<p><a \"Interactive genome Plot\" </a> </p>\n", file=genome.page, append=T)
	cat("<p><a <img height=",image.height, paste("src= \"",genome.image.file,"\" width=", sep=""),image.width, " usemap=#map> </a> </p>\n", file=genome.page, append=T)
	cat("<map name=map>\n", file=genome.page, append=T)
	chrom.coords <- (break.coords*(image.width*0.7))+(image.width*0.2)
	for(i in 2:length(chrom.coords))
	{	
	coords.i <- paste(chrom.coords[i],(image.height*0.1),chrom.coords[i-1],(image.height*0.9), sep=",")
	cat("<area shape=rect coords=",coords.i,chrom.links[i-1],chroms.alt[i-1], "</p>\n", file=genome.page, append=T)
	}
	cat("</map>", file=genome.page, append=T)
	cat("</body>", file=genome.page, append=T)
	
	cat("<map name=map>\n", file=genome.map, append=T)
	
	chrom.coords <- (break.coords*(image.width*0.8))+(image.width*0.1)
	for(i in 2:length(chrom.coords))
	{	
	coords.i <- paste(chrom.coords[i],(image.height*0.1),chrom.coords[i-1],(image.height*0.9), sep=",")
	cat("<area shape=rect coords=",coords.i,chrom.links[i-1],chroms.alt[i-1], "</p>\n", file=genome.map, append=T)
	}
	cat("</map>", file=genome.map, append=T)	
	
	for(j in 1:length(chroms))
	{
		chr.image.file <- paste("chr",chlabs[j],"png", sep=".")
		chr.image <- paste(case.folder,paste("chr",chlabs[j],"png", sep="."),sep="/")
		png(chr.image, width=image.width, height=image.height)
		
		cgh.chrom <- cgh[is.element(fData(cgh)$chrom,j),]
		MB <- (fData(cgh.chrom)$start+fData(cgh.chrom)$end)/2000000
		
		yAxis.chrom <- c(min(range(exprs(cgh.chrom))[1]-0.2,-1), max(range(exprs(cgh.chrom))[2]+0.2,1))
		genomePlotGrid(cgh, chroms=chroms[j], case=case, yAxis=yAxis.chrom, colourGL=colourGL)
		dev.off()

	probe.ids <- fData(cgh.chrom)[,match(probe.id, fvarLabels(cgh))]
	window.starts <- fData(cgh.chrom)$start - (MB.window*1000000)
	window.ends <- fData(cgh.chrom)$end + (MB.window*1000000)
#	clone.ids.alt <- paste("alt=",probe.ids," title=",probe.ids, " />", sep="")
	clone.maps <- paste(chlabs[j],(paste(fData(cgh.chrom)$start,fData(cgh.chrom)$end,sep="-")), sep=":")
	clone.maps.window <- paste(chlabs[j],(paste(window.starts,window.ends,sep="-")), sep=":")
	ensembl.root <- "href=http://www.ensembl.org/Homo_sapiens/Location/View?r="
	ensembl.links <- paste(ensembl.root, clone.maps, sep="")
	clone.ids.alt <- paste("alt=",paste(probe.ids,clone.maps,sep=" - ")," title=",paste(probe.ids,clone.maps,sep=" - "), " />", sep="")
	
	chr.html <- paste(case.folder,paste("chr",chlabs[j],"html", sep="."),sep="/")
	cat("<body>", file=chr.html, append=F)
	cat("<p><a \"Interactive genome Plot\" </a> </p>\n", file=chr.html, append=T)
	cat("<p><a <img height=",image.height, paste("src= \"",chr.image.file,"\" width=", sep=""),image.width, "usemap=#map> </a> </p>\n", file=chr.html, append=T)
	cat("<map name=map>\n", file=chr.html, append=T)
	MB.coords <- ((MB/MB[length(MB)])*(image.width*0.7)) +(image.width*0.2)
	
	for( i in 1:length(MB))
	{	
	if(i == 1) {coords.i <- paste(MB.coords[1],(image.height*0.2),MB.coords[1],(image.height*0.9), sep=",")}else{
		if (i == length(MB)) {coords.i <- paste(MB.coords[i],(image.height*0.2),MB.coords[i],(image.height*0.9), sep=",")}else
			{coords.i <- paste((MB.coords[i-1]+MB.coords[i])/2,(image.height*0.2),(MB.coords[i]+MB.coords[i+1])/2,(image.height*0.9), sep=",")}}
	cat("<area shape=rect coords=",coords.i,ensembl.links[i],clone.ids.alt[i], "</p>\n", file=chr.html, append=T)
	}
	
	cat("</map>", file=chr.html, append=T)
	cat("</body>", file=chr.html, append=T)
	
	chr.map <- paste(case.folder,paste("chr",chlabs[j],"map", sep="."),sep="/")
	
	cat("<map name=map>\n", file=chr.map, append=T)
	MB.coords <- ((MB/MB[length(MB)])*(image.width*0.7)) +(image.width*0.2)
	for( i in 1:length(MB))
	{	
	if(i == 1) {coords.i <- paste(MB.coords[1],(image.height*0.2),MB.coords[1],(image.height*0.9), sep=",")}else{
		if (i == length(MB)) {coords.i <- paste(MB.coords[i],(image.height*0.2),MB.coords[i],(image.height*0.9), sep=",")}else
			{coords.i <- paste((MB.coords[i-1]+MB.coords[i])/2,(image.height*0.2),(MB.coords[i]+MB.coords[i+1])/2,(image.height*0.9), sep=",")}}
	cat("<area shape=rect coords=",coords.i,ensembl.links[i],clone.ids.alt[i], "</p>\n", file=chr.map, append=T)
	}
	cat("</map>", file=chr.map, append=T)
	
	cat(chroms[j], "\t")
	}
	cat(sampleNames(cgh)[case], "\n")
	
	}
}


genomePlotCGHExpression <- function(cgh, main=NULL, bw=15, chrom, cgh.case, thresh=NULL, yAxis=NULL, expression)
{
	
	require(Biobase)
	if(!is.numeric(cgh.case)) {
		cgh.case <- match(cgh.case, sampleNames(cgh))
		if(is.na(cgh.case)) stop("SampleName is not found in the CGH ExpressionSet")
		}
	if(!is.element(cgh.case,1:ncol(cgh))) stop("Case is not found in the CGH ExpressionSet")
	
	cgh <- cgh[which(is.element(fData(cgh)$chrom,chrom)),]

	if(is.null(main))
	{
	main<- sampleNames(cgh)[cgh.case]
	}
	
	cgh.MB <- ((fData(cgh)$start +fData(cgh)$end)/2000000)
	cgh.M <- exprs(cgh)[,cgh.case]
	
	GL <- assayData(cgh)$GL[,cgh.case]
	
	exp.case <- match(sampleNames(cgh)[cgh.case], sampleNames(expression))
	if(is.na(exp.case)) stop("Case is not found in the Expression ExpressionSet")
	expression <- removeUnmappedProbes(expression, verbose=F)
	expression <- expression[order(fData(expression)$chrom, fData(expression)$start),]

	expression <- expression[which(is.element(fData(expression)$chrom,chrom)),]

	exp.MB <- ((fData(expression)$start +fData(expression)$end)/2000000)
	exp.M <- exprs(expression)[,exp.case]
	exp.M <- exp.M*(max(cgh.M)/(max(exp.M)))
	
	if(is.null(yAxis))
		{
		yAxis <- c(
		min(
		round(min(exp.M, na.rm=T))-0.5,
		round(min(cgh.M, na.rm=T))-0.5,
		-1), 
		max(
		round(max(exp.M, na.rm=T))+0.5,
		round(max(cgh.M, na.rm=T))+0.5,
		1))
		}

	yscale <- seq(-30,30,1)
	yscale <- yscale[yscale >= yAxis[1] & yscale <= yAxis[2]]
		
	gains.MB <- losses.MB <- amps.MB <- dels.MB <- cgh.MB[is.finite(cgh.M)]
	gains <- losses <- amps <- dels <- cgh.M[is.finite(cgh.M)]
	
	gains.MB[GL<=0] <- NA
	gains[GL<=0] <- NA
	losses.MB[GL>=0] <- NA
	losses[GL>=0] <- NA
		
	amps.MB[GL<=1] <- NA
	amps[GL<=1] <- NA
	dels.MB[GL>= -1] <- NA
	dels[GL>= -1] <- NA
	
	ylab <- "log ratio"

	xlab <- paste("chromosome", chrom, "- MBp")
	plot(cgh.MB,cgh.M,ylim=yAxis,axes=T,xlab=xlab,					main=main,ylab=ylab, pch=16, cex=0.5, las=1, type="n")
	abline(v=seq(from=min(cgh.MB),to=max(cgh.MB), by=10), col="grey", lty=8)
	abline(h=seq(from=min(yscale),to=max(yscale), by=1), col="grey")
	abline(h=seq(from=0.5,to=max(yscale), by=1),lty=8, col="grey")
	abline(h=seq(from=-0.5,to=min(yscale), by=-1),lty=8, col="grey")
	abline(h=0, col="grey")
	if(!is.null(thresh)){abline(h=c(-thresh,thresh),col="darkgrey", lty=8, lwd=2)}
	points(cgh.MB,cgh.M, pch=16, cex=0.5)
	points(exp.MB,exp.M, pch=16, cex=0.5, col="blue")
	
	if(is.element("smo", assayDataElementNames(cgh)))
		{
		lines(cgh.MB,assayData(cgh)$smo[,cgh.case],col="purple", lwd=2)
		par(xpd=T)
		}

}

plotFrequency <- function(cgh, chroms=NULL, yAxis.GL=NULL, yAxis.AD=NULL, yAxis.amp=NULL, device=NULL, main=NULL, project)
{
	require(Biobase)
	if(!is.null(chroms))
		{
		cgh <- cgh[which(is.element(fData(cgh)$chrom,chroms)),]
		}

	fData(cgh)$MB <- calculateMB(cgh)
		
	chlabs <- as.character(c(1:22,"X","Y"))
	chrom.MB <- (fData(cgh)$end+fData(cgh)$start)/2000000
	
	gain.f <- apply(assayData(cgh)$GL,1,function(x) sum(x >= 1))/ncol((cgh))
	loss.f <- -1*(apply(assayData(cgh)$GL,1,function(x) sum(x <= -1))/ncol((cgh)))
	amp.f <- apply(assayData(cgh)$GL,1,function(x) sum(x == 2))/ncol((cgh))
	del.f <- -1*(apply(assayData(cgh)$GL,1,function(x) sum(x == -2))/ncol((cgh)))
	
	if(is.null(yAxis.GL)) yAxis.GL <- c(min(round(range(loss.f),1),-0.1),max(round(range(gain.f)+0.1, 1),0.1))
	if(is.null(yAxis.AD)) yAxis.AD <- c(min(round(range(del.f),1),-0.1),max(round(range(amp.f), 1)+0.1,0.1))
	if(is.null(yAxis.amp)) yAxis.amp <- c(0,max(round(range(amp.f)+0.1, 1),0.1))
	
	gain.f[gain.f ==0] <- NA
	loss.f[loss.f ==0] <- NA
	amp.f[amp.f ==0] <- NA
	del.f[del.f ==0] <- NA
	
	chr.nms <- names(table(fData(cgh)$chrom))
	chr.nms[match("23",chr.nms,nomatch=0)] <- "X"
	chr.nms[match("24",chr.nms,nomatch=0)] <- "Y"
	ends <- sapply(split(fData(cgh)$MB,fData(cgh)$chrom),max)
	begs <- sapply(split(fData(cgh)$MB,fData(cgh)$chrom),min)
	mids <- (begs+ends)/2
	nums <- seq(-1,1,0.2)
	MB.labs <- seq(0, max(chrom.MB),10)
	MB.lines <- seq(0, max(chrom.MB),5)
	
	if(!is.null(chroms))
	{
		project <- paste(project,"chr",sep=".")
		for(i in 1:length(chroms))
		{project <- paste(project,chroms[i],sep=".")}
	}
	
	gain.title <- paste(main," Frequency of Gains and Losses")
	amp.title <- paste(main," Frequency of Amplifications and Deletions")
	amp.only.title <- paste(main," Frequency of Amplifications")
	
	if(length(device) !=0){
	if(!is.element("Frequency.plots", list.files())) dir.create("Frequency.plots")
			
	if(device == "quartz") {quartz(title=paste(project,gain.title),width=10, height=6)}
	if(device == "PNG")
		{GL.file <- paste(project,"GL.frequency","png", sep=".")
		png(file=paste("Frequency.plots/",GL.file,sep=""),width=2400, height=1200, pointsize=24, type="cairo")}
		
	if(device == "JPEG")
		{GL.file <- paste(project,"GL.frequency","jpeg", sep=".")
		jpeg(file=paste("Frequency.plots/",GL.file,sep=""),width=2400, height=1200, pointsize=24, quality=100)}
		
		
	if(device == "PDF")
		{GL.file <- paste(project,"GL.frequency","pdf", sep=".")
		pdf(file=paste("Frequency.plots/",GL.file,sep=""),width=12, height=6, pointsize=12)}
	if(device == "PS")
		{GL.file <- paste(project,"GL.frequency","ps", sep=".")
		postscript(file=paste("Frequency.plots/",GL.file,sep=""),width=10, height=5, pointsize=8)}}
	
	if(length(chroms)==1){
		plot(chrom.MB,gain.f,type="n",ylab="Frequency",axes=F,
		xlab=paste("Chromosome",chroms,"- MB"),ylim=yAxis.GL,main=gain.title)
		abline(h=seq(from=0.1,to=max(yAxis.GL), by=0.1),lty=8, col="grey")
		abline(h=seq(from=-0.1,to=min(yAxis.GL), by=-0.1),lty=8, col="grey")
		points(chrom.MB,gain.f,type="h",col="dark green", lwd=2)
		points(chrom.MB,loss.f,type="h",col="red", lwd=2)
		axis(side=1,at=MB.labs,tick=T,cex=0.6)
		axis(side=2,at=nums,labels=nums,tick=T,cex=0.6, las=1)
		abline(v=MB.lines,col=16)
		abline(h=0) 
		box()}else{
		plot(fData(cgh)$MB,gain.f,type="n",ylab="Frequency",
		axes=F,xlab="Chromosome",ylim=yAxis.GL,main=gain.title)
		abline(h=seq(from=0.1,to=max(yAxis.GL), by=0.1),lty=8, col="grey")
		abline(h=seq(from=-0.1,to=min(yAxis.GL), by=-0.1),lty=8, col="grey")
		points(fData(cgh)$MB,gain.f,type="h",col="dark green")
		points(fData(cgh)$MB,loss.f,type="h",col="red")
		axis(side=1,at=mids,labels=chr.nms,tick=F,cex=0.6)
		axis(side=2,at=nums,labels=nums,tick=T,cex=0.6, las=1)
		abline(v=c(1,ends),col=16)
		abline(h=0) 
		box()}
		
	if(length(device) != 0 && device != "quartz") dev.off()
	
	if(length(device) !=0){	
	if(device == "quartz") {quartz(title=paste(project,amp.title),width=10, height=6)}
	if(device == "PNG")
		{AD.file <- paste(project,"AD.frequency","png", sep=".")
		png(file=paste("Frequency.plots/",AD.file,sep=""),width=2400, height=1200, pointsize=24, type="cairo")}
		
	if(device == "JPEG")
		{AD.file <- paste(project,"AD.frequency","jpeg", sep=".")
		jpeg(file=paste("Frequency.plots/",AD.file,sep=""),width=2400, height=1200, pointsize=24, quality=100)}
		
		
	if(device == "PDF")
		{AD.file <- paste(project,"AD.frequency","pdf", sep=".")
		pdf(file=paste("Frequency.plots/",AD.file,sep=""),width=12, height=6, pointsize=12)}
	if(device == "PS")
		{AD.file <- paste(project,"AD.frequency","ps", sep=".")
		postscript(file=paste("Frequency.plots/",AD.file,sep=""),width=10, height=5, pointsize=8)}}

	if(length(chroms)==1){
		plot(chrom.MB,amp.f,type="n",ylab="Frequency",axes=F,
		xlab=paste("Chromosome",chroms,"- MB"),ylim=yAxis.GL,main=gain.title)
		abline(h=seq(from=0.1,to=max(yAxis.AD), by=0.1),lty=8, col="grey")
		abline(h=seq(from=0.1,to=min(yAxis.AD), by=-0.1),lty=8, col="grey")
		points(chrom.MB,amp.f,type="h",col="dark green")
		points(chrom.MB,del.f,type="h",col="red")
		axis(side=1,at=MB.labs,tick=T,cex=0.6)
		axis(side=2,at=nums,labels=nums,tick=T,cex=0.6, las=1)
		abline(v=MB.lines,col=16)
		abline(h=0) 
		box()}else{
		plot(fData(cgh)$MB,amp.f,type="n",ylab="Frequency",
		axes=F,xlab="Chromosome",ylim=yAxis.AD,main=amp.title)
		abline(h=seq(from=0.1,to=max(yAxis.AD), by=0.1),lty=8, col="grey")
		abline(h=seq(from=-0.1,to=min(yAxis.AD), by=-0.1),lty=8, col="grey")
		points(fData(cgh)$MB,amp.f,type="h",col="dark green")
		points(fData(cgh)$MB,del.f,type="h",col="red")
		axis(side=1,at=mids,labels=chr.nms,tick=F,cex=0.6)
		axis(side=2,at=nums,labels=nums,tick=T,cex=0.6, las=1)
		abline(v=c(1,ends),col=16)
		abline(h=0) 
		box()}
	
	if(length(device) != 0 && device != "quartz") dev.off()
	
	
	if(length(device) !=0){	
	if(device == "quartz") {quartz(title=paste(project,amp.only.title),width=10, height=6)}
	if(device == "PNG")
		{A.file <- paste(project,"Amp.frequency","png", sep=".")
		png(file=paste("Frequency.plots/",A.file,sep=""),width=2400, height=1200, pointsize=12, type="cairo")}
		
		
	if(device == "JPEG")
		{A.file <- paste(project,"Amp.frequency","jpeg", sep=".")
		jpeg(file=paste("Frequency.plots/",A.file,sep=""),width=2400, height=1200, pointsize=24, quality=100)}
	
	if(device == "PDF")
		{A.file <- paste(project,"Amp.frequency","pdf", sep=".")
		pdf(file=paste("Frequency.plots/",A.file,sep=""),width=12, height=6, pointsize=12)}
	if(device == "PS")
		{A.file <- paste(project,"Amp.frequency","ps", sep=".")
		postscript(file=paste("Frequency.plots/",A.file,sep=""),width=10, height=5, pointsize=8)}}

	if(length(chroms)==1){
		plot(chrom.MB,amp.f,type="n",ylab="Frequency",axes=F,
		xlab=paste("Chromosome",chroms,"- MB"),ylim=yAxis.amp,main=amp.only.title)
		abline(h=seq(from=0.1,to=max(yAxis.amp), by=0.1),lty=8, col="grey")
		points(chrom.MB,amp.f,type="h",col="dark green")
		axis(side=1,at=MB.labs,tick=T,cex=0.6)
		axis(side=2,at=nums,labels=nums,tick=T,cex=0.6, las=1)
		abline(v=MB.lines,col=16)
		abline(h=0) 
		box()}else{
		plot(fData(cgh)$MB,amp.f,type="n",ylab="Frequency",
		axes=F,xlab="Chromosome",ylim=yAxis.amp,main=amp.only.title)
		abline(h=seq(from=0.1,to=max(yAxis.amp), by=0.1),lty=8, col="grey")
		points(fData(cgh)$MB,amp.f,type="h",col="dark green")
		axis(side=1,at=mids,labels=chr.nms,tick=F,cex=0.6)
		axis(side=2,at=nums,labels=nums,tick=T,cex=0.6, las=1)
		abline(v=c(1,ends),col=16)
		abline(h=0) 
		box()}
	
	if(length(device) != 0 && device != "quartz") dev.off()

}

plotFrequencyLattice <- function(cgh, pheno, chroms=NULL, project=NULL, device=NULL)
{
	require(lattice)
	require(Biobase)
	
	cgh <- cgh[,!is.na(pheno)]
	pheno <- factor(pheno[!is.na(pheno)])
	
	if(!is.null(chroms))
		{
		cgh <- cgh[which(is.element(fData(cgh)$chrom,chroms)),]
		}
	fData(cgh)$MB <- calculateMB(cgh)
	chrom.MB <- (fData(cgh)$end+fData(cgh)$start)/2000000
	
	cgh.i <- gain.f <- loss.f <- amp.f <- del.f <- list(NULL)
	for(i in 1:length(levels(pheno))){
		cgh.i[[i]] <- cgh[,which(pheno == levels(pheno)[i])]
		
		gain.f[[i]] <- apply(assayData(cgh.i[[i]])$GL,1, function(x) sum(x >= 1))/ncol((cgh.i[[i]]))
		loss.f[[i]] <- -1*(apply(assayData(cgh.i[[i]])$GL,1, function(x) sum(x <= -1))/ncol((cgh.i[[i]])))
		amp.f[[i]] <- apply(assayData(cgh.i[[i]])$GL,1, function(x) sum(x >= 2))/ncol((cgh.i[[i]]))
		del.f[[i]] <- -1*(apply(assayData(cgh.i[[i]])$GL,1, function(x) sum(x <= -2))/ncol((cgh.i[[i]])))
		}
	
	chr.nms <- names(table(fData(cgh)$chrom))
	chr.nms[match("23",chr.nms,nomatch=0)] <- "X"
	chr.nms[match("24",chr.nms,nomatch=0)] <- "Y"
	ends <- sapply(split(fData(cgh)$MB,fData(cgh)$chrom),max)
	begs <- sapply(split(fData(cgh)$MB,fData(cgh)$chrom),min)
	mids <- (begs+ends)/2
	nums <- seq(-1,1,0.2)
	MB.labs <- seq(0, max(chrom.MB),10)
	MB.lines <- seq(0, max(chrom.MB),5)
	
	breaks <- c(ends,max(chrom.MB))-c(0,begs)
	
	if(!is.element("Frequency.plots", list.files())) dir.create("Frequency.plots", showWarnings=F)
	
	freq.matrix.GL <- data.frame(GL=rep("Gains and Losses",nrow(cgh)*length(levels(pheno))),gains=unlist(gain.f), losses=unlist(loss.f), MB=fData(cgh)$MB, pheno=rep(levels(pheno),each=nrow(cgh)))
	freq.matrix.GL[freq.matrix.GL == 0] <- NA
######## CN: fixed "darkgreen" instead of "dark green"

	GL.plot <- xyplot(freq.matrix.GL$gains+freq.matrix.GL$losses~freq.matrix.GL$MB|freq.matrix.GL$pheno, type="h", layout=c(1,length(levels(pheno))),
	col=c("darkgreen","red"), main=paste(project,"- Frequency of gains and losses"), xlab="Chromosome", ylab="Frequency",
	scales=list(x=list(at=mids,labels=chr.nms)),
	panel = function(...)
	{panel.abline(h = 0, lty = 1)
	panel.abline(h=seq(from=-1,to=1, by=0.1),lty=8, col="grey")
	panel.abline(v=c(-10,ends+10),col="grey")
	panel.xyplot(...)})
	
	if(length(device) !=0){
	if(device == "quartz") {quartz(title=paste(project,gain.title),width=6, height=10)}
	if(device == "PNG")
		{GL.file <- paste(project,"GL.Lattice.Frequency","png", sep=".")
		trellis.device("png", color=T, file=paste("Frequency.plots/",GL.file,sep=""),width=2400, height=1200, pointsize=24, type="cairo")}
	if(device == "JPEG")
		{GL.file <- paste(project,"GL.Lattice.Frequency","jpeg", sep=".")
		trellis.device("jpeg", color=T, file=paste("Frequency.plots/",GL.file,sep=""),width=2400, height=1200, pointsize=24, quality=100)}
	if(device == "PDF")
		{GL.file <- paste(project,"GL.Lattice.Frequency","pdf", sep=".")
		trellis.device("pdf", color=T, file=paste("Frequency.plots/",GL.file,sep=""),width=8, height=6, pointsize=8)}
	if(device == "PS")
		{GL.file <- paste(project,"GL.Lattice.Frequency","ps", sep=".")
		trellis.device("postscript", color=T, file=paste("Frequency.plots/",GL.file,sep=""),width=10, height=5, pointsize=8)}
	}
	plot(GL.plot)
	if(length(device) != 0 && device != "quartz") dev.off()
		
	freq.matrix.AD <- data.frame(GL=rep("Gains and Losses",nrow(cgh)*length(levels(pheno))),amps=unlist(amp.f), dels=unlist(del.f), MB=fData(cgh)$MB, pheno=rep(levels(pheno),each=nrow(cgh)))
	freq.matrix.AD[freq.matrix.AD == 0] <- NA
	AD.plot <- xyplot(freq.matrix.AD$amps+freq.matrix.AD$dels~freq.matrix.AD$MB|freq.matrix.AD$pheno, type="h", layout=c(1,length(levels(pheno))),
	col=c("dark green","red"), main=paste(project,"- Frequency of amplifications and deletions"), xlab="Chromosome", ylab="Frequency",
	scales=list(x=list(at=mids,labels=chr.nms)),
	panel = function(...)
	{panel.abline(h = 0, lty = 1)
	panel.abline(h=seq(from=-1,to=1, by=0.1),lty=8, col="grey")
	panel.abline(v=c(-10,ends+10),col="grey")
	panel.xyplot(...)})
	
	if(length(device) !=0){
	if(device == "quartz") {quartz(title=paste(project,gain.title),width=6, height=10)}
	if(device == "PNG")
		{AD.file <- paste(project,"AD.Lattice.Frequency","png", sep=".")
		trellis.device("png", color=T, file=paste("Frequency.plots/",AD.file,sep=""),width=2400, height=1200, pointsize=12, type="cairo")}
	if(device == "JPEG")
		{AD.file <- paste(project,"AD.Lattice.Frequency","jpeg", sep=".")
		trellis.device("jpeg", color=T, file=paste("Frequency.plots/",AD.file,sep=""),width=2400, height=1200, pointsize=12, quality=100)}
	if(device == "PDF")
		{AD.file <- paste(project,"AD.Lattice.Frequency","pdf", sep=".")
		trellis.device("pdf", color=T, file=paste("Frequency.plots/",AD.file,sep=""),width=8, height=6, pointsize=8)}
	if(device == "PS")
		{AD.file <- paste(project,"AD.Lattice.Frequency","ps", sep=".")
		trellis.device("postscript", color=T, file=paste("Frequency.plots/",AD.file,sep=""),width=10, height=5, pointsize=8)}
	}
	
	plot(AD.plot)
	if(length(device) != 0 && device != "quartz") dev.off()
}


frequencyShingles <- function(cgh, main, project, device=NULL)
{

require(Biobase)
if(is.null(assayData(cgh)$GL)) stop("Call aCGH gains and losses first\n")
MB <- calculateMB(cgh)
gain.mids <- loss.mids <- amp.mids <- del.mids <- as.list(rep(NA,ncol(cgh)))
gain.widths <- loss.widths <- amp.widths <- del.widths <- as.list(rep(NA,ncol(cgh)))

cat("Listing aCGH states\n")

for(i in 1:ncol(cgh)){
iGL <- assayData(cgh)$GL[,i]
gain.iGL <- as.numeric(iGL >= 1)
amp.iGL <- as.numeric(iGL == 2)
loss.iGL <- as.numeric(iGL <= -1)
del.iGL <- as.numeric(iGL == -2)

gains.table <- listBreaks(cgh=cgh, GL.column=gain.iGL, contig=3)
gain.breaks <- gains.table[,which(gains.table[2,] == 1), drop=F]
amp.table <- listBreaks(cgh=cgh, GL.column=amp.iGL, contig=3)
amp.breaks <- amp.table[,which(amp.table[2,] == 1), drop=F]
loss.table <- listBreaks(cgh=cgh, GL.column=loss.iGL, contig=3)
loss.breaks <- loss.table[,which(loss.table[2,] == 1), drop=F]
del.table <- listBreaks(cgh=cgh, GL.column=del.iGL, contig=3)
del.breaks <- del.table[,which(del.table[2,] == 1), drop=F]

if(max(gains.table[2,]) != 0){
gain.ends <- gain.breaks[1,]
gain.starts <- gain.breaks[1,]-gain.breaks[3,]+1
gain.mids[[i]]<- ((MB[gain.ends]+MB[gain.starts])/2)/max(MB)
gain.widths[[i]] <- (MB[gain.ends]-MB[gain.starts])/max(MB)
}
if(max(loss.table[2,]) != 0){
loss.ends <- loss.breaks[1,]
loss.starts <- loss.breaks[1,]-loss.breaks[3,]+1
loss.mids[[i]]<- ((MB[loss.ends]+MB[loss.starts])/2)/max(MB)
loss.widths[[i]] <- (MB[loss.ends]-MB[loss.starts])/max(MB)
}
if(max(amp.table[2,]) != 0){
amp.ends <- amp.breaks[1,]
amp.starts <- amp.breaks[1,]-amp.breaks[3,]+1
amp.mids[[i]]<- ((MB[amp.ends]+MB[amp.starts])/2)/max(MB)
amp.widths[[i]] <- (MB[amp.ends]-MB[amp.starts])/max(MB)
}
if(max(del.table[2,]) != 0){
del.ends <- del.breaks[1,]
del.starts <- del.breaks[1,]-del.breaks[3,]+1
del.mids[[i]]<- ((MB[del.ends]+MB[del.starts])/2)/max(MB)
del.widths[[i]] <- (MB[del.ends]-MB[del.starts])/max(MB)
}
cat(i, "\t")
}


	chrom.maps <- split(MB,fData(cgh)$chrom)
		ends <- sapply(chrom.maps,function(x) max(x))
		begs <- sapply(chrom.maps,function(x) min(x))
		midMB <- (begs + ends)/2
		mid.coords <- midMB/max(MB)
		end.coords <- ends/max(MB)
		beg.coords <- begs/max(MB)
		break.coords <- (ends[1:(length(ends)-1)] + begs[2:length(begs)])/2/max(MB)
	
	chlabs <- as.character(unique(fData(cgh)$chrom))
		chlabs[which(chlabs==23)] <- "X"
		chlabs[which(chlabs==24)] <- "Y"	
	
	library(grid)
	
	if(length(device) !=0){
	if(!is.element("frequency.shingles", list.files())) dir.create("frequency.shingles")
			
	if(device == "quartz") {quartz(title=paste(project,"Gains and losses"),width=10, height=6)}
	if(device == "PNG")
		{GL.file <- paste(project,"GL.frequency.shingle","png", sep=".")
		png(file=paste("frequency.shingles/",GL.file,sep=""),width=1600, height=800, pointsize=12)}
		
	if(device == "JPEG")
		{GL.file <- paste(project,"GL.frequency.shingle","jpeg", sep=".")
		jpeg(file=paste("frequency.shingles/",GL.file,sep=""),width=800, height=400, pointsize=12)}
		
		
	if(device == "PDF")
		{GL.file <- paste(project,"GL.frequency.shingle","pdf", sep=".")
		pdf(file=paste("frequency.shingles/",GL.file,sep=""),width=12, height=6, pointsize=12)}
	if(device == "PS")
		{GL.file <- paste(project,"GL.frequency.shingle","ps", sep=".")
		postscript(file=paste("frequency.shingles/",GL.file,sep=""),width=10, height=5, pointsize=8)}}
	
	grid.newpage()
	
	gain.vp <- viewport(x = unit(0.1, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(gain.vp)
	grid.rect()
	for(i in 1:ncol(cgh)){
	if(!is.na(gain.mids)[i]) grid.rect(x=gain.mids[[i]], y=i/(ncol(cgh)+1), width=gain.widths[[i]],
	height=1/(ncol(cgh)+1), gp = gpar(fill = "dark green", col="dark green"))
	}
	
	for(i in 1:length(break.coords))
		{
		grid.lines(break.coords[i],c(0,1),gp=gpar(col="black", lwd=1), default.units="npc")
		}
		
	popViewport()
	
	#	viewport for labels
	gain.samples.vp <- viewport(x = unit(0.8, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(0.2, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(gain.samples.vp)
	#grid.rect()
	for(i in 1:ncol(cgh)){
		grid.text(sampleNames(cgh)[i], x=0.05, y=i/(ncol(cgh)+1), hjust=0, gp=gpar(cex=0.5))}
	popViewport()
		
	chlabs.vp <- viewport(x = unit(0.1, "npc"),
	y = unit(0.45,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.05, "npc"),
	just = c("left", "bottom"))
	pushViewport(chlabs.vp)
	#grid.rect()
	for(i in 1:length(unique(fData(cgh)$chrom)))
		{
		grid.text(chlabs[i], x=unit(mid.coords[i],"npc"),y = unit(0.5, "npc"), gp = gpar(cex=0.7))
		}
	popViewport()
	
	loss.vp <- viewport(x = unit(0.1, "npc"),
	y = unit(0.05,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(loss.vp)
	grid.rect()
	for(i in 1:ncol(cgh)){
	if(!is.na(loss.mids)[i]) grid.rect(x=loss.mids[[i]], y=i/(ncol(cgh)+1), width=loss.widths[[i]],
	height=1/(ncol(cgh)+1), gp = gpar(fill = "dark red", col="dark red"))
	}
	for(i in 1:length(break.coords))
		{
		grid.lines(break.coords[i],c(0,1),gp=gpar(col="black", lwd=1), default.units="npc")
		}
	popViewport()
	
	loss.samples.vp <- viewport(x = unit(0.8, "npc"),
	y = unit(0.05,  "npc"),
	width = unit(0.2, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(loss.samples.vp)
	#grid.rect()
	for(i in 1:ncol(cgh)){
		grid.text(sampleNames(cgh)[i], x=0.05, y=i/(ncol(cgh)+1), hjust=0, gp=gpar(cex=0.5))}
	popViewport()
	
	gain.axis.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(0.1, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(gain.axis.vp)
	#grid.rect()
	grid.text("Gains", rot=90)
	popViewport()
	
	loss.axis.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.05,  "npc"),
	width = unit(0.1, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(loss.axis.vp)
	#grid.rect()
	grid.text("Losses", rot=90)
	popViewport()
	
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.9,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "bottom"))
	pushViewport(title.vp)
	#grid.rect()
	grid.text(main)
	popViewport()
	
if(length(device) != 0 && device != "quartz") dev.off()

if(length(device) !=0){
	if(!is.element("frequency.shingles", list.files())) dir.create("frequency.shingles")
			
	if(device == "quartz") {quartz(title=paste(project,"Amps and Dels"),width=10, height=6)}
	if(device == "PNG")
		{AD.file <- paste(project,"AD.frequency.shingle","png", sep=".")
		png(file=paste("frequency.shingles/",AD.file,sep=""),width=1600, height=800, pointsize=12)}
		
	if(device == "JPEG")
		{AD.file <- paste(project,"AD.frequency.shingle","jpeg", sep=".")
		jpeg(file=paste("frequency.shingles/",AD.file,sep=""),width=800, height=400, pointsize=12)}
		
		
	if(device == "PDF")
		{AD.file <- paste(project,"AD.frequency.shingle","pdf", sep=".")
		pdf(file=paste("frequency.shingles/",AD.file,sep=""),width=12, height=6, pointsize=12)}
	if(device == "PS")
		{AD.file <- paste(project,"AD.frequency.shingle","ps", sep=".")
		postscript(file=paste("frequency.shingles/",AD.file,sep=""),width=10, height=5, pointsize=8)}}
		
	grid.newpage()
	
	amp.vp <- viewport(x = unit(0.1, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(amp.vp)
	grid.rect()
	for(i in 1:ncol(cgh)){
	if(!is.na(amp.mids)[i]) grid.rect(x=amp.mids[[i]], y=i/(ncol(cgh)+1), width=amp.widths[[i]],
	height=1/(ncol(cgh)+1), gp = gpar(fill = "dark green", col="dark green"))
	}
	
	for(i in 1:length(break.coords))
		{
		grid.lines(break.coords[i],c(0,1),gp=gpar(col="black", lwd=1), default.units="npc")
		}
		
	popViewport()
	
	#	viewport for labels
	amp.samples.vp <- viewport(x = unit(0.8, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(0.2, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(amp.samples.vp)
	#grid.rect()
	for(i in 1:ncol(cgh)){
		grid.text(sampleNames(cgh)[i], x=0.05, y=i/(ncol(cgh)+1), hjust=0, gp=gpar(cex=0.5))}
	popViewport()
		
	chlabs.vp <- viewport(x = unit(0.1, "npc"),
	y = unit(0.45,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.05, "npc"),
	just = c("left", "bottom"))
	pushViewport(chlabs.vp)
	#grid.rect()
	for(i in 1:length(unique(fData(cgh)$chrom)))
		{
		grid.text(chlabs[i], x=unit(mid.coords[i],"npc"),y = unit(0.5, "npc"), gp = gpar(cex=0.7))
		}
	popViewport()
	
	del.vp <- viewport(x = unit(0.1, "npc"),
	y = unit(0.05,  "npc"),
	width = unit(0.7, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(del.vp)
	grid.rect()
	for(i in 1:ncol(cgh)){
	if(!is.na(del.mids)[i]) grid.rect(x=del.mids[[i]], y=i/(ncol(cgh)+1), width=del.widths[[i]],
	height=1/(ncol(cgh)+1), gp = gpar(fill = "dark red", col="dark red"))
	}
	for(i in 1:length(break.coords))
		{
		grid.lines(break.coords[i],c(0,1),gp=gpar(col="black", lwd=1), default.units="npc")
		}
	popViewport()
	
	del.samples.vp <- viewport(x = unit(0.8, "npc"),
	y = unit(0.05,  "npc"),
	width = unit(0.2, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(del.samples.vp)
	#grid.rect()
	for(i in 1:ncol(cgh)){
		grid.text(sampleNames(cgh)[i], x=0.05, y=i/(ncol(cgh)+1), hjust=0, gp=gpar(cex=0.5))}
	popViewport()
	
	amp.axis.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.5,  "npc"),
	width = unit(0.1, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(amp.axis.vp)
	#grid.rect()
	grid.text("Amplifications", rot=90)
	popViewport()
	
	del.axis.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.05,  "npc"),
	width = unit(0.1, "npc"),
	height = unit(0.4, "npc"),
	just = c("left", "bottom"))
	pushViewport(del.axis.vp)
	#grid.rect()
	grid.text("Deletions", rot=90)
	popViewport()
	
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.9,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "bottom"))
	pushViewport(title.vp)
	#grid.rect()
	grid.text("Frequency Shingles")
	popViewport()
	
	if(length(device) != 0 && device != "quartz") dev.off()
	
}


readBeads <- function(beadstudio, ann.library=NULL, norm.method="rsn", flags.on=F, pval.thresh=0.01)
{
	require(Biobase)
	require(lumi)
	require(lumiHumanAll.db)
	require(annotate)
	
	cat("\nReading in BeadStudio output...\n")
	lumi.raw <- lumiR(beadstudio)
	cat("Done\n")
	cat("\nVariance stabilisation normalisation...\n")
	lumi.T <- lumiT(lumi.raw)
	cat("Done\n")
	cat("\nNormalisation - method = ", norm.method, "\n")
	lumi.N <- lumiN(lumi.T, method=norm.method)
	cat("Done\n")
	cat("\nLumibatch summary\n")
	summary(lumi.N)
	cat("QC summary\n")
	summary(lumi.N, "QC")
	
	BS.file <- read.delim(beadstudio, sep="\t", stringsAsFactors=F, header=T, na.strings=c("", " ", "NA", "N/A", "#N/A"))
	
	if(!is.null(BS.file$PROBE_SEQUENCE)) {
		cat(sum(duplicated(BS.file$PROBE_SEQUENCE)), "Duplicated sequences\n")
		BS.file <- BS.file[which(!duplicated(BS.file$PROBE_SEQUENCE)),]
		lumi.N <- lumi.N[which(!duplicated(BS.file$PROBE_SEQUENCE)),]
		}
		
	if(!is.null(BS.file$TRANSCRIPT)) {cat("Retrieving TRANSCRIPT IDs...")
		fData(lumi.N)$TRANSCRIPT <- BS.file$TRANSCRIPT
		cat("Done\n")}
	
	if(!is.null(BS.file$PROBE_ID)) {cat("Retrieving PROBE IDs...")
		fData(lumi.N)$probeID <- BS.file$PROBE_ID
		cat("Done\n")}
	
	if(is.null(ann.library))
	{
	if(!is.null(BS.file$CHROMOSOME)) {cat("Retrieving mapping positions...")
		chrom <- BS.file$CHROMOSOME
		chrom[chrom == "X"] <- "23"
		chrom[chrom == "Y"] <- "24"
		fData(lumi.N)$chrom <- as.numeric(chrom)}
		
	if(!is.null(BS.file$PROBE_COORDINATES)) {chrloc <- BS.file$PROBE_COORDINATES
	start <- numeric(length(chrloc))
	end <- numeric(length(chrloc))
	for(i in 1:length(chrloc))
		{
		coords <- unlist(strsplit((unlist(strsplit(chrloc[i],":"))),"-"))
			if(!is.null(coords))
			{
			start[i] <- min(as.numeric(coords))
			end[i] <- max(as.numeric(coords))	
			}else{
			start[i] <- NA
			end[i] <- NA	
			}	
		}
	fData(lumi.N)$start <- start
	fData(lumi.N)$end <- end
	cat("Done\n")}
	
	if(!is.null(BS.file$ACCESSION)) {cat("Retrieving accession numbers...")
		fData(lumi.N)$accession <- BS.file$ACCESSION
		cat("Done\n")}
	
	
	if(!is.null(BS.file$SYMBOL)) {cat("Retrieving gene symbols...")
		fData(lumi.N)$symbol <- BS.file$SYMBOL
		cat("Done\n")}
	
	
	if(!is.null(BS.file$DEFINITION)) {cat("Retrieving genenames...")
		fData(lumi.N)$description <- BS.file$DEFINITION
		cat("Done\n")}

	if(!is.null(BS.file$REFSEQ_ID)) {cat("Retrieving refseq...")
		fData(lumi.N)$refseq <- BS.file$REFSEQ_ID
		cat("Done\n")}
	
	if(!is.null(BS.file$UNIGENE_ID)) {cat("Retrieving unigene...")
		fData(lumi.N)$unigene <- BS.file$UNIGENE_ID
		cat("Done\n")}
	
	if(!is.null(BS.file$ENTREZ_GENE_ID)) {cat("Retrieving entrez...")
		fData(lumi.N)$entrez <- BS.file$ENTREZ_GENE_ID
		cat("Done\n")}
	
	if(!is.null(BS.file$PROBE_SEQUENCE_ID)) {cat("Retrieving sequences...")
		fData(lumi.N)$sequence <- BS.file$PROBE_SEQUENCE_ID
		cat("Done\n")}
	
	if(!is.null(BS.file$ARRAY_ADDRESS_ID)) {cat("Retrieving array address...")
		fData(lumi.N)$address <- BS.file$ARRAY_ADDRESS_ID
		cat("Done\n")}
	
	if(!is.null(BS.file$CYTOBAND)) {cat("Retrieving cytoband...")
		fData(lumi.N)$cytoband <- BS.file$CYTOBAND
		cat("Done\n")}
	
	
}else{
	cat("Retrieving mapping positions...")
	chrom.ann <- as.character(lookUp(featureNames(lumi.N),ann.library,'CHR'))
	chrom.ann[chrom.ann == "X"] <- "23"
	chrom.ann[chrom.ann == "Y"] <- "24"
	fData(lumi.N)$chrom <- as.numeric(chrom.ann)
	chrloc.ann <- lookUp(featureNames(lumi.N),ann.library,'CHRLOC')
	start.ann <- numeric(length(chrloc.ann))
	end.ann <- numeric(length(chrloc.ann))
	for(i in 1:length(chrloc.ann))
		{
		start.ann[i] <- min(abs(unlist(chrloc.ann[i])))
		end.ann[i] <- max(abs(unlist(chrloc.ann[i])))		}
	fData(lumi.N)$start <- start.ann
	fData(lumi.N)$end <- end.ann
	cat("Done\n")
	cat("Retrieving genenames...")
	fData(lumi.N)$description <- as.character(lookUp(featureNames(lumi.N),ann.library,'GENENAME'))
	cat("Done\n")
	cat("Retrieving gene symbols...")
	fData(lumi.N)$symbol <- as.character(lookUp(featureNames(lumi.N),ann.library,'SYMBOL'))
	cat("Done\n")
	cat("Retrieving entrezIDs...")
	fData(lumi.N)$entrez <- as.character(lookUp(featureNames(lumi.N),ann.library,'ENTREZID'))
	cat("Done\n")
	cat("Retrieving accession numbers...")
	fData(lumi.N)$refseq <- as.character(lookUp(featureNames(lumi.N),ann.library,'ACCNUM'))
	cat("Done\n")
	cat("Retrieving Ensembl IDs...")
	fData(lumi.N)$ensg <- as.character(lookUp(featureNames(lumi.N),ann.library,'ENSEMBL'))
	cat("Done\n")
	cat("Retrieving Unigene IDs...")
	unigene <- lookUp(featureNames(lumi.N),ann.library,'UNIGENE')
	unigene.one <- character(length(unigene))
	for(i in 1:length(unigene))
		{
		unigene.one[i] <- unlist(unigene[i])[1]
		}
	fData(lumi.N)$unigene <- unigene.one
	cat("Done\n")
	cat("Retrieving Cytobands...")
	fData(lumi.N)$cytoband <- as.character(lookUp(featureNames(lumi.N),ann.library,'MAP'))
	cat("Done\n")
}

#	make a matrix of present absent calls based upon p.value detection
flags <- detectionCall(lumi.N, type="matrix", Th=pval.thresh)
flags <- flags=="A"

#	If flags.on=T remove intensity values from the M matrix and replace with NA

if(flags.on) is.na(exprs(lumi.N)) <- flags

lumi.out <- as(lumi.N, "ExpressionSet")
assayDataElement(lumi.out, "flags") <- flags
lumi.out <- as(lumi.out, "BACE.exp")
lumi.out
}

readBeadsMethyLumi <- function(Beadstudio, pval.thresh=0.05, flags.on=F, beta.cuts = c(0.2, 0.8), mapfun = c("atan", "ratio"))
{
	require(Biobase)
	require(lumi)
	require(methylumi)
	
	
	cat("Reading BeadStudio methylation data\n")
	
	metLum <- methylumiR(Beadstudio)
	
	pdata <- data.frame(sampleNames=sampleNames(metLum))
	row.names(pdata) <- sampleNames(metLum)
	phenoData(metLum) <- new("AnnotatedDataFrame", pdata)
	
	flags <- matrix(F, nrow(metLum), ncol(metLum))
	flags[assayData(metLum)$pvals >  pval.thresh] <- T
	
	names(fData(metLum))[names(fData(metLum))=="SYMBOL"] <- "symbol"
	fData(metLum)$chrom <- fData(metLum)$CHROMOSOME
	fData(metLum)$chrom <- as.character(fData(metLum)$chrom)
	
	fData(metLum)$chrom[fData(metLum)$chrom == "X"] <- 23
	fData(metLum)$chrom[fData(metLum)$chrom == "Y"] <- 24
	fData(metLum)$chrom[!is.element(fData(metLum)$chrom,1:24)] <- NA
	fData(metLum)$chrom <- as.numeric(fData(metLum)$chrom)
	fData(metLum)$start <- fData(metLum)$CPG_COORDINATE
	fData(metLum)$end <- fData(metLum)$start + unlist(lapply(fData(metLum)$INPUT_SEQUENCE,nchar))
	
	if(flags.on) is.na(betas(metLum)) <- flags
	assayDataElement(metLum, "flags") <- flags
	
	return(metLum)

}




normalizeMethyLumiBACE <- function (eset, beta.cuts = c(0.2, 0.8), mapfun = c("atan", "ratio")) 
{
    require(Biobase)
    mapfun = match.arg(mapfun)
    
    good <- rep(TRUE, ncol(eset))
    cy3 <- assayData(eset)$unmethylated
    cy3[cy3 < 0] <- NA
    cy5 <- assayData(eset)$methylated
    cy5[cy5 < 0] <- NA
    for (i in 1:ncol(cy5)) {
        cy3inc <- (!is.na(exprs(eset)[, i]) & !is.na(cy3[, i]))
        cy5inc <- (!is.na(exprs(eset)[, i]) & !is.na(cy5[, i]))
        cy3vec <- cy3[cy3inc, i]
        cy5vec <- cy5[cy5inc, i]
        cy3h <- median(cy3vec[exprs(eset)[cy3inc, i] < beta.cuts[1]])
        cy3l <- median(cy3vec[exprs(eset)[cy3inc, i] > beta.cuts[2]])
        cy5l <- median(cy5vec[exprs(eset)[cy5inc, i] < beta.cuts[1]])
        cy5h <- median(cy5vec[exprs(eset)[cy5inc, i] > beta.cuts[2]])
        corfactor <- (cy3h - cy3l)/(cy5h - cy5l)
        cy5[, i] <- cy5[, i] * (corfactor)
        cy5vec <- cy5[cy5inc, i]
        newcy5l <- median(cy5vec[exprs(eset)[cy5inc, i] < beta.cuts[1]])
        if (newcy5l < cy3l) {
            cy5[, i] <- cy5[, i] + (cy3l - newcy5l)
        }
        else {
            cy3[, i] <- cy3[, i] + (newcy5l - cy3l)
        }
        if (corfactor < 0) {
            good[i] <- FALSE
            warning(sprintf("Sample %d has medians that do not make sense for a normal sample\n(cy3l=%f ,cy5l=%f ,cy3h=%f ,cy5h=%f)\nRemoving sample!  Check quality control.", 
                i, cy3l, cy5l, cy3h, cy5h))
            cy5[, i] <- NA
            cy3[, i] <- NA
        }
    }
    newbeta <- 0
    if (mapfun == "atan") {
        newbeta <- atan((cy5)/(cy3))/(pi/2)
    }
    else {
        newbeta <- cy5/(cy5 + cy3 + 100)
    }
    
    assayDataElement(eset, "unmethylated") <- cy3
    assayDataElement(eset, "methylated") <-  cy5
    assayDataElement(eset, "exprs") <- newbeta
    
    return(eset)
}



callStatesMeth <- function(meth, thresholds=c(0.2,0.8))
{
	
	require(Biobase)
	mstates <- matrix(0, nrow(meth), ncol(meth))
	mstates[exprs(meth) < thresholds[1]] <- -1
	mstates[exprs(meth) > thresholds[2]] <- 1
	
	assayDataElement(meth,"mstates") <- mstates
	
	fData(meth)$UM.count <- apply(assayData(meth)$mstates, 1, function(x) sum(x == -1))
	UM.samples <- apply(assayData(meth)$mstates, 1, function(x) which(x == -1))
	fData(meth)$UM.sample.ids <- rep(NA, nrow(assayData(meth)$mstates))
	for(i in 1:length(UM.samples)){fData(meth)$unmethylated.samples[i] <- paste(sampleNames(meth)[unlist(UM.samples[i])],collapse=", ")}
	
	fData(meth)$PM.count <- apply(assayData(meth)$mstates, 1, function(x) sum(x == 0))
	PM.samples <- apply(assayData(meth)$mstates, 1, function(x) which(x == 0))
	fData(meth)$PM.sample.ids <- rep(NA, nrow(assayData(meth)$mstates))
	for(i in 1:length(PM.samples)){fData(meth)$partially.methylated.samples[i] <- paste(sampleNames(meth)[unlist(PM.samples[i])],collapse=", ")}
	
	fData(meth)$M.count <- apply(assayData(meth)$mstates, 1, function(x) sum(x == 1))
	M.samples <- apply(assayData(meth)$mstates, 1, function(x) which(x == 1))
	fData(meth)$M.sample.ids <- rep(NA, nrow(assayData(meth)$mstates))
	for(i in 1:length(M.samples)){fData(meth)$methylated.samples[i] <- paste(sampleNames(meth)[unlist(M.samples[i])],collapse=", ")}
	
	meth
}


expressionMethOverlay <- function(eset.expression, eset.meth, eset.id="symbol", meth.id="symbol")
#	Function to combine meth score matrix and probne info with an expression object
{
	#	index to put meth in the right order
	require(Biobase)
	cat("Combining expression and mehtylation ExpressionSets\n")
	expression.index <- which(fvarLabels(eset.expression) == eset.id)
	if(length(expression.index) != 1) stop("eset.id not matched")
	meth.index <- which(fvarLabels(eset.meth) == meth.id)
	if(length(meth.index) != 1) stop("meth.id not matched")
	
	#	probes in the meth which are in the expression
	expression.probes <- fData(eset.expression)[,expression.index]
	meth.probes <- fData(eset.meth)[,meth.index]
	meth.probes.in.exp <- meth.probes[which(is.element(meth.probes,expression.probes))]
	
	#	expand if there are replicates
	exp.probe.index <- pmatch(meth.probes.in.exp,expression.probes, duplicates.ok=T, nomatch=0)
	
	meth.samples.in.exp <- sampleNames(eset.meth)[which(is.element(sampleNames(eset.meth), sampleNames(eset.expression)))]
	exp.samples.index <- pmatch(meth.samples.in.exp, sampleNames(eset.expression), duplicates.ok=T, nomatch=0)
	
	#	subset meth to those present in the expression
	eset.meth <- eset.meth[meth.probes %in% expression.probes, sampleNames(eset.meth) %in% sampleNames(eset.expression)]
	expression.featurenames <- featureNames(eset.expression)
	#	subset expression to those present in the meth and order the hybs
	eset.expression <- eset.expression[exp.probe.index,exp.samples.index]
	featureNames(eset.expression) <- make.names(expression.featurenames[exp.probe.index], unique=T)

	#	cat the overlay
	if(all(sampleNames(eset.expression) == sampleNames(eset.meth)) ) {cat(ncol(eset.expression),"matching sample.ids\n")}else{stop("sample.ids do not match")}
	if(all(fData(eset.expression)[,expression.index] == fData(eset.meth)[,meth.index])) {cat(nrow(eset.expression),"probe.ids in common\n")}else{stop("probe.ids do not match")}
	
	fvarLabels(eset.meth) <- paste("meth", fvarLabels(eset.meth), sep=".")
	fData(eset.expression) <- cbind.data.frame(fData(eset.expression), fData(eset.meth))
	validObject(eset.expression)
	meth.exprs <- exprs(eset.meth)
	row.names(meth.exprs) <- featureNames(eset.expression)
	assayDataElement(eset.expression, "betas") <- meth.exprs
	meth.flags <- assayData(eset.meth)$flags
	row.names(meth.flags) <- featureNames(eset.expression)
	assayDataElement(eset.expression, "meth.flags") <- meth.flags
	
	meth.methylated <- assayData(eset.meth)$methylated
	row.names(meth.methylated) <- featureNames(eset.expression)
	assayDataElement(eset.expression, "methylated") <- meth.methylated
	meth.unmethylated <- assayData(eset.meth)$unmethylated
	row.names(meth.unmethylated) <- featureNames(eset.expression)
	assayDataElement(eset.expression, "unmethylated") <- meth.unmethylated
	
	
	if(is.null(assayData(eset.meth)$mstates)){eset.meth <- callStatesMeth(eset.meth)}else{
		meth.mstates <- assayData(eset.meth)$mstates
		row.names(meth.mstates) <- featureNames(eset.expression)
		assayDataElement(eset.expression, "mstates") <- meth.mstates}
	
	eset.expression
}

expressionMethCorrelation <- function(exp.meth, correlation.method="pearson", adj.method="BH", pval.thresh=0.05)
#	exp.meth is a list object with matched data.frames for expression and methylation
#	exp.meth is a result of medianMeth?
#	
{
	options(warn=-1)

#	commands <- c(exp.meth$commands,deparse(match.call()))
	require(Biobase)
	cat("Correlation of gene expression with methylation values\n")
		
	meth.pearson.p <- meth.pearson.cor <- meth.spearman.p <- meth.spearman.cor <- meth.glm.p <- meth.glm.cor <- rep(NA,nrow(exp.meth))
	
	for (i in 1:nrow(exp.meth))
	
		{
			if((sum(is.finite(exprs(exp.meth)[i,])) > 3) && (sum(is.finite(assayData(exp.meth)$beta[i,])) > 3))
			{
#				
				if(correlation.method=="pearson"){
				pearson.score <- cor.test(exprs(exp.meth)[i,],assayData(exp.meth)$betas[i,], method="pearson", na.rm=T)
				meth.pearson.p[i] <- pearson.score$p.value
				meth.pearson.cor[i] <- pearson.score$estimate}
				if(correlation.method=="spearman"){
				spearman.score <- cor.test(exprs(exp.meth)[i,],assayData(exp.meth)$betas[i,], method="spearman", na.rm=T)
				meth.spearman.p[i]<- spearman.score$p.value
				meth.spearman.cor[i] <- spearman.score$estimate}
				if(correlation.method=="glm"){
				glm.summary <- summary(glm(assayData(exp.meth)$betas[i,]~exprs(exp.meth)[i,], family=quasibinomial(link="logit")), mustart=c(0,1))
				meth.glm.p[i] <- glm.summary$coefficients[nrow(glm.summary$coefficients),4]
				meth.glm.cor[i] <- glm.summary$coefficients[nrow(glm.summary$coefficients),1]
				}
#				cat(i," ")
#				if(i %% 100 == 0) cat(i," ")
			}
		
		}
		
	cat("\nDone\n")
	
	if(correlation.method=="pearson"){
	fData(exp.meth)$meth.pearson.p <- meth.pearson.p
	fData(exp.meth)$meth.pearson.adjp <- p.adjust(meth.pearson.p, method=adj.method)
	fData(exp.meth)$meth.pearson.cor <- meth.pearson.cor
#	sort.order <- fData(exp.meth)$meth.pearson.p
	cat(length(which(fData(exp.meth)$meth.pearson.adjp < pval.thresh)),"significantly correlated genes\n")
	}
	if(correlation.method=="spearman"){
	fData(exp.meth)$meth.spearman.p <- meth.spearman.p
	fData(exp.meth)$meth.spearman.adjp <- p.adjust(meth.spearman.p, method=adj.method)
	fData(exp.meth)$meth.spearman.cor <- meth.spearman.cor
#	sort.order <- fData(exp.meth)$meth.spearman.p
	cat(length(which(fData(exp.meth)$meth.spearman.adjp < pval.thresh)),"significantly correlated genes\n")
	}
	if(correlation.method=="glm"){
	fData(exp.meth)$meth.glm.p <- meth.glm.p
	fData(exp.meth)$meth.glm.adjp <- p.adjust(meth.glm.p, method=adj.method)
	fData(exp.meth)$meth.glm.cor <- meth.glm.cor
#	sort.order <- fData(exp.meth)$meth.spearman.p
	cat(length(which(fData(exp.meth)$meth.glm.adjp < pval.thresh)),"significantly correlated genes\n")
	}
#	exp.meth$commands <- commands
	cat("\nDone\n")

	exp.meth
}


wilcoxTestExpMeth <- function(exp.meth, count=2, pval.thresh=0.05, adj.method="BH")
#	Function to run Wilcoxon rank sum tests 
#	to correlate gene expression with methylation states
#	Runs on the product of expressionMethOverlay
{

#	commands <- c(exp.meth$commands,deparse(match.call()))
	require(Biobase)
	cat("Wilcoxon rank sum tests of expression based on methylation states\n")

	M.wilcox.p <- UM.wilcox.p <- M.fold <- UM.fold <- rep(NA,nrow(exp.meth))
	
	for (i in 1:nrow(exp.meth))
		{

		meth <- which(assayData(exp.meth)$mstates[i,] == 1)
		if((sum(!is.na(exprs(exp.meth)[i,meth])) >= count) && (sum(!is.na(exprs(exp.meth)[i,-meth])) >= count)){
			M.wilcox.p[i] <- wilcox.test(exprs(exp.meth)[i,meth],exprs(exp.meth)[i,-meth])$p.value}
		
		unmeth <- which(assayData(exp.meth)$mstates[i,] == -1)
		if((sum(!is.na(exprs(exp.meth)[i,unmeth])) >= count) && (sum(!is.na(exprs(exp.meth)[i,-unmeth])) >= count)){
			UM.wilcox.p[i] <- wilcox.test(exprs(exp.meth)[i,unmeth],exprs(exp.meth)[i,-unmeth])$p.value}
		
		M.fold[i] <- 2^(mean(exprs(exp.meth)[i,meth],na.rm=T))/2^(mean(exprs(exp.meth)[i,-meth],na.rm=T))
		UM.fold[i] <- 2^(mean(exprs(exp.meth)[i,unmeth],na.rm=T))/2^(mean(exprs(exp.meth)[i,-unmeth],na.rm=T))

		if(i %% 100 == 0) cat(i," ")
		
		}

	cat("\nDone\n")
	
	M.wilcox.p <- M.wilcox.p
	
	fData(exp.meth)$M.fold <- M.fold
	fData(exp.meth)$UM.fold <- UM.fold
	
	fData(exp.meth)$M.wilcox.p <- M.wilcox.p
	fData(exp.meth)$M.wilcox.adjp <- p.adjust(M.wilcox.p, method=adj.method)
	fData(exp.meth)$UM.wilcox.p <- UM.wilcox.p
	fData(exp.meth)$UM.wilcox.adjp <- p.adjust(UM.wilcox.p, method=adj.method)	
#	exp.meth$commands <- commands
	cat("\nDone\n")
	exp.meth
}


wilcoxTestCGHMeth <- function(cgh.meth, threshold.data=T, count=2, pval.thresh=0.05)
#	Function to run Wilcoxon rank sum tests 
#	to correlate gene expression with aCGH states
#	Runs on the product of medianCGH
#	converts GL matrix into 0/1 counts for gain loss amplification and deletion
#	using cgh$thresholds
#	runs wilcox test on each row of M and GL matrices
#	for GL.gain, GL.loss, GL.amp, GL.del
{
	require(Biobase)
	require(multtest)
#	commands <- c(cgh.meth$commands,deparse(match.call()))
	cat("Wilcoxon rank sum tests of expression based on copy number changes\n")
#	cgh.meth <- subsetMA(cgh.meth, probes = order(fData(cgh.meth)$chrom, fData(cgh.meth)$start))
	wilcox.p <- matrix(NA,nrow(cgh.meth),4)
	gain.fold <- loss.fold <- amp.fold <- del.fold <- rep(NA,nrow(cgh.meth))
	
	if(threshold.data == T)
	{
			
	for (i in 1:nrow(cgh.meth))
		{
		
		gains <- which(assayData(cgh.meth)$GL[i,] >= 1)
		if((sum(!is.na(assayData(cgh.meth)$beta[i,gains])) >= count) && (sum(!is.na(assayData(cgh.meth)$betas[i,-gains])) >= count)){
			wilcox.p[i,1] <- wilcox.test(assayData(cgh.meth)$betas[i,gains],assayData(cgh.meth)$betas[i,-gains])$p.value}
		
		losses <- which(assayData(cgh.meth)$GL[i,] <= -1)
		if((sum(!is.na(assayData(cgh.meth)$betas[i,losses])) >= count) && (sum(!is.na(assayData(cgh.meth)$betas[i,-losses])) >= count)){
			wilcox.p[i,2] <- wilcox.test(assayData(cgh.meth)$betas[i,losses],assayData(cgh.meth)$betas[i,-losses])$p.value}
		
		amps <- which(assayData(cgh.meth)$GL[i,]== 2)
		if((sum(!is.na(assayData(cgh.meth)$betas[i,amps])) >= count) && (sum(!is.na(assayData(cgh.meth)$betas[i,-amps])) >= count)){
			wilcox.p[i,3] <- wilcox.test(assayData(cgh.meth)$betas[i,amps],assayData(cgh.meth)$betas[i,-amps])$p.value}
		
		dels <- which(assayData(cgh.meth)$GL[i,]== -2)
		if((sum(!is.na(assayData(cgh.meth)$betas[i,dels])) > count) && (sum(!is.na(assayData(cgh.meth)$betas[i,-dels])) > count)){
			wilcox.p[i,4] <- wilcox.test(assayData(cgh.meth)$beta[i,dels],assayData(cgh.meth)$beta[i,-dels])$p.value}
		
		gain.fold[i] <- 2^(mean(assayData(cgh.meth)$betas[i,gains],na.rm=T))/2^(mean(assayData(cgh.meth)$betas[i,-gains],na.rm=T))
		loss.fold[i] <- 2^(mean(assayData(cgh.meth)$betas[i,losses],na.rm=T))/2^(mean(assayData(cgh.meth)$betas[i,-losses],na.rm=T))
		amp.fold[i] <- 2^(mean(assayData(cgh.meth)$betas[i,amps],na.rm=T))/2^(mean(assayData(cgh.meth)$betas[i,-amps],na.rm=T))
		del.fold[i] <- 2^(mean(assayData(cgh.meth)$betas[i,dels],na.rm=T))/2^(mean(assayData(cgh.meth)$betas[i,-dels],na.rm=T))
		
		if(i %% 100 == 0) cat(i," ")
		
		}
	}
	
	if(threshold.data == F)
	{
	for (i in 1:nrow(cgh.meth))
		{

		if((sum(!is.na(assayData(cgh.meth)$betas[i,assayData(cgh.meth)$GL[i,] == 1])) >= count) && (sum(!is.na(assayData(cgh.meth)$betas[i,assayData(cgh.meth)$GL[i,] == 0])) >= count))
		{wilcox.p[i,1] <- wilcox.test(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 1)],assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 0)])$p.value
			}
		if((sum(!is.na(assayData(cgh.meth)$betas[i,assayData(cgh.meth)$GL[i,] == -1])) >= count) && (sum(!is.na(assayData(cgh.meth)$betas[i,assayData(cgh.meth)$GL[i,] == 0])) >= count))
		{wilcox.p[i,2] <- wilcox.test(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== -1)],assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 0)])$p.value
			}
		if((sum(!is.na(assayData(cgh.meth)$betas[i,assayData(cgh.meth)$GL[i,] == 2])) >= count) && (sum(!is.na(assayData(cgh.meth)$betas[i,assayData(cgh.meth)$GL[i,] == 0])) >= count))
		{wilcox.p[i,3] <- wilcox.test(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 2)],assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 0)])$p.value
			}
		if((sum(!is.na(assayData(cgh.meth)$betas[i,assayData(cgh.meth)$GL[i,] == -2])) >= count) && (sum(!is.na(assayData(cgh.meth)$betas[i,assayData(cgh.meth)$GL[i,] == 0])) >= count))
		{wilcox.p[i,4] <- wilcox.test(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== -2)],assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 0)])$p.value
			}
		
		gain.fold[i] <- 2^(mean(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 1)],na.rm=T))/2^(mean(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 0)],na.rm=T))
		loss.fold[i] <- 2^(mean(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== -1)],na.rm=T))/2^(mean(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 0)],na.rm=T))
		amp.fold[i] <- 2^(mean(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 2)],na.rm=T))/2^(mean(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 0)],na.rm=T))
		del.fold[i] <- 2^(mean(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== -2)],na.rm=T))/2^(mean(assayData(cgh.meth)$betas[i,which(assayData(cgh.meth)$GL[i,]== 0)],na.rm=T))
		
		if(i %% 100 == 0) cat(i," ")
		}
	}
	
	cat("\nDone\n")
	
	wilcox.adjp <- wilcox.p
	
	glad <- countGLAD(cgh.meth)[,1:4]
	
	
	for(i in 1:4) glad[which(assayData(cgh.meth)$GLAD[,i] >= count),i] <- 1
	
	gain.breaks <- listBreaks(cgh.meth, GL.column=glad[,1], contig=1)
	gain.breaks <- gain.breaks[,which(gain.breaks[2,] > 0), drop=F]
	loss.breaks <- listBreaks(cgh=cgh.meth, GL.column=glad[,2], contig=1)
	loss.breaks <- loss.breaks[,which(loss.breaks[2,] > 0) ,drop=F]
	amp.breaks <- listBreaks(cgh=cgh.meth, GL.column=glad[,3], contig=1)
	amp.breaks <- amp.breaks[,which(amp.breaks[2,] > 0),drop=F]
	del.breaks <- listBreaks(cgh=cgh.meth, GL.column=glad[,4], contig=1)
	del.breaks <- del.breaks[,which(del.breaks[2,] > 0),drop=F]
	
	if(ncol(gain.breaks) > 0) for (i in 1:ncol(gain.breaks))
		{		
		end.index <- gain.breaks[1,i]
		start.index <- gain.breaks[1,i]-(gain.breaks[3,i]-1)
		pvals.i <- wilcox.p[start.index:end.index,1]
		if(sum(!is.na(pvals.i)) >1)
		{wilcox.adjp[start.index:end.index,1] <- p.adjust(pvals.i, method="BH")}
		}
	if(ncol(loss.breaks) > 0) for (i in 1:ncol(loss.breaks))
		{		
		end.index <- loss.breaks[1,i]
		start.index <- loss.breaks[1,i]-(loss.breaks[3,i]-1)
		pvals.i <- wilcox.p[start.index:end.index,2]
		if(sum(!is.na(pvals.i)) >1)
		{wilcox.adjp[start.index:end.index,2] <- p.adjust(pvals.i, method="BH")}
		}
	if(ncol(amp.breaks) > 0) for (i in 1:ncol(amp.breaks))
		{		
		end.index <- amp.breaks[1,i]
		start.index <- amp.breaks[1,i]-(amp.breaks[3,i]-1)
		pvals.i <- wilcox.p[start.index:end.index,3]
		if(sum(!is.na(pvals.i)) >1)
		{wilcox.adjp[start.index:end.index,3] <- p.adjust(pvals.i, method="BH")}
		}
	if(ncol(del.breaks) > 0) for (i in 1:ncol(del.breaks))
		{		
		end.index <- del.breaks[1,i]
		start.index <- del.breaks[1,i]-(del.breaks[3,i]-1)
		pvals.i <- wilcox.p[start.index:end.index,4]
		if(sum(!is.na(pvals.i)) >1)
		{wilcox.adjp[start.index:end.index,4] <- p.adjust(pvals.i, method="BH")}
		}
	
	fData(cgh.meth)$beta.gain.fold <- gain.fold
	fData(cgh.meth)$beta.loss.fold <- loss.fold
	fData(cgh.meth)$beta.amp.fold <- amp.fold
	fData(cgh.meth)$beta.del.fold <- del.fold
	
	if(threshold.data == T)
	{
		fData(cgh.meth)$Wilcox.beta.p.gain <- wilcox.p[,1]
		fData(cgh.meth)$Wilcox.beta.adjp.gain <- wilcox.adjp[,1]
		fData(cgh.meth)$Wilcox.beta.p.loss <- wilcox.p[,2]
		fData(cgh.meth)$Wilcox.beta.adjp.loss <- wilcox.adjp[,2]
		fData(cgh.meth)$Wilcox.beta.p.amp <- wilcox.p[,3]
		fData(cgh.meth)$Wilcox.beta.adjp.amp <- wilcox.adjp[,3]
		fData(cgh.meth)$Wilcox.beta.p.del <- wilcox.p[,4]
		fData(cgh.meth)$Wilcox.beta.adjp.del <- wilcox.adjp[,4]

	}
	
	if(threshold.data == F)
	{
		fData(cgh.meth)$Wilcox.beta.p.gain.vs.NC <- wilcox.p[,1]
		fData(cgh.meth)$Wilcox.beta.adjp.gain.vs.NC <- wilcox.adjp[,1]
		fData(cgh.meth)$Wilcox.beta.p.loss.vs.NC <- wilcox.p[,2]
		fData(cgh.meth)$Wilcox.beta.adjp.loss.vs.NC <- wilcox.adjp[,2]
		fData(cgh.meth)$Wilcox.beta.p.amp.vs.NC <- wilcox.p[,3]
		fData(cgh.meth)$Wilcox.beta.adjp.amp.vs.NC <- wilcox.adjp[,3]
		fData(cgh.meth)$Wilcox.beta.p.del.vs.NC <- wilcox.p[,4]
		fData(cgh.meth)$Wilcox.beta.adjp.del.vs.NC <- wilcox.adjp[,4]
		
	}

#	cgh.meth$commands <- commands

	cgh.meth
	
}

listBreaks <- function(cgh, GL.column, contig)
{
	#	a script to take integral data  from a cgh object as a vector of states (-2,-1,0,1,2)
	#	record the positions at which the input states change
	#	calculate the lengths of the contiguous regions
	require(Biobase)
	
	if(length(unique(fData(cgh)$chrom)) >1){
		GL.chroms <- split(GL.column, fData(cgh)$chrom)
		
	breaks <- rle(GL.chroms[[1]])
	index <- cumsum(breaks$lengths)
	values <- breaks$values
	lengths <- breaks$lengths
	chrom <- unique(fData(cgh)$chrom)[1]
	breaks.table <- rbind(index,values,lengths,chrom)
	start.state <- c(1,GL.chroms[[1]][1],1,fData(cgh)$chrom[index[1]])
	breaks.table <- cbind(start.state,breaks.table)
	rownames(breaks.table) <-c("index", "values", "lengths", "chrom")
	
	for (i in 2:length(GL.chroms))
			{
			start.index <- breaks.table[1,ncol(breaks.table)]
			breaks <- rle(GL.chroms[[i]])
			index <- cumsum(breaks$lengths)+start.index
			values <- breaks$values
			lengths <- breaks$lengths
			chrom <- unique(fData(cgh)$chrom)[i]
			ibreaks.table <- rbind(index,values,lengths,chrom)
			start.state <- c(start.index+1,GL.chroms[[i]][1],1,fData(cgh)$chrom[index[1]])
			ibreaks.table <- cbind(start.state,ibreaks.table)
			rownames(ibreaks.table) <-c("index", "values", "lengths", "chrom")
			breaks.table <- cbind(breaks.table,ibreaks.table)
			}
	breaks.table <- breaks.table[,which(breaks.table[3,] >= contig),drop=F]
	
	}else{
	breaks <- rle(GL.column)
	index <- cumsum(breaks$lengths)
	values <- breaks$values
	lengths <- breaks$lengths
	chrom <- unique(fData(cgh)$chrom)[1]
	breaks.table <- rbind(index,values,lengths,chrom)
	start.state <- c(1,GL.column[1],1,fData(cgh)$chrom[index[1]])
	breaks.table <- cbind(start.state,breaks.table)
	rownames(breaks.table) <-c("index", "values", "lengths", "chrom")
		
	}
			
	breaks.table

}


listBreaksGL <- function(cgh, contig=3, acgh.cnv.table="acgh.cnvs.txt", full.cnv.table="full.cnvs.txt", genes.table="all.genes.txt", cytoband.table="cytobands.txt", mirnas.table="mirnas.txt", project)
#	Function to list aCGH copy number changes in a CGH object
#	Lists gains losses amplifications and deletions in individual cases
#	Lists contiguous probes longer than a specified region (contig)
#	Uses spitGALTables to output regions into results files for each case
{
	require(Biobase)
	if(!is.element("sample.GALS", list.files())) dir.create("sample.GALS")
	for(i in 1:ncol(cgh))
	{
		
	gain.column <- loss.column <- amp.column <- del.column <- rep(0, nrow(cgh))
	gain.column[assayData(cgh)$GL[,i] >= 1] <- 1
	loss.column[assayData(cgh)$GL[,i] <= -1] <- 1
	amp.column[assayData(cgh)$GL[,i] == 2] <- 1
	del.column[assayData(cgh)$GL[,i] == -2] <- 1
	
	gain.breaks.table <- listBreaks(cgh=cgh[,i], GL.column=gain.column, contig=contig)
	gain.breaks <- gain.breaks.table[,which(gain.breaks.table[2,] == 1),drop=F]
	loss.breaks.table <- listBreaks(cgh=cgh[,i], GL.column=loss.column, contig=contig)
	loss.breaks <- loss.breaks.table[,which(loss.breaks.table[2,] == 1),drop=F]
	
	amp.breaks.table <- listBreaks(cgh=cgh[,i], GL.column=amp.column, contig=contig)
	amp.breaks <- amp.breaks.table[,which(amp.breaks.table[2,] == 1), drop=F]
	del.breaks.table <- listBreaks(cgh=cgh[,i], GL.column=del.column, contig=contig)
	del.breaks <- del.breaks.table[,which(del.breaks.table[2,] == 1), drop=F]
	breaks <- list(gain.breaks, loss.breaks, amp.breaks, del.breaks)
	names(breaks) <- c("gain.breaks", "loss.breaks", "amp.breaks", "del.breaks")
	
	filename <- paste(project, i, sampleNames(cgh)[i], "GALS", "xls", sep=".")
	spitGALTables(cgh=cgh[,i], gain.count=0, loss.count=0, amp.count=0, del.count=0, contig=contig, breaks=breaks, file=paste("sample.GALS/",filename,sep=""), acgh.cnv.table=acgh.cnv.table, full.cnv.table=full.cnv.table, genes.table=genes.table, cytoband.table=cytoband.table, mirnas.table=mirnas.table)
	cat(sampleNames(cgh)[i], "\n")
	}
}


listBreaksGLAD <- function(cgh, contig=3, gain.count, loss.count=NULL, amp.count=NULL, del.count=NULL, acgh.cnv.table="acgh.cnvs.txt", full.cnv.table="full.cnvs.txt", genes.table="all.genes.txt", cytoband.table="cytobands.txt", mirnas.table="mirnas.txt", project)
#	Function to list aCGH copy number changes in a CGH object
{
	require(Biobase)
	if(is.null(loss.count)) loss.count <- gain.count
	if(is.null(amp.count)) amp.count <- gain.count
	if(is.null(del.count)) del.count <- amp.count
	
	gain.counts <- as.vector(apply(assayData(cgh)$GL >= 1, 1, sum))
	loss.counts <- as.vector(apply(assayData(cgh)$GL <= -1, 1, sum))
	amp.counts <- as.vector(apply(assayData(cgh)$GL == 2, 1, sum))
	del.counts <- as.vector(apply(assayData(cgh)$GL == -2, 1, sum))
	
	gain.counts[gain.counts < gain.count] <- 0
	gain.counts[gain.counts >= gain.count] <- 1
	loss.counts[loss.counts < loss.count] <- 0
	loss.counts[loss.counts >= loss.count] <- 1
	amp.counts[amp.counts < amp.count] <- 0
	amp.counts[amp.counts >= amp.count] <- 1
	del.counts[del.counts < del.count] <- 0
	del.counts[del.counts >= del.count] <- 1
	
	gain.breaks <- listBreaks(cgh=cgh, GL.column=gain.counts, contig=contig)
	gain.breaks <- gain.breaks[,which(gain.breaks[2,] > 0), drop=F]
	loss.breaks <- listBreaks(cgh=cgh, GL.column=loss.counts, contig=contig)
	loss.breaks <- loss.breaks[,which(loss.breaks[2,] > 0) ,drop=F]
	amp.breaks <- listBreaks(cgh=cgh, GL.column=amp.counts, contig=contig)
	amp.breaks <- amp.breaks[,which(amp.breaks[2,] > 0),drop=F]
	del.breaks <- listBreaks(cgh=cgh, GL.column=del.counts, contig=contig)
	del.breaks <- del.breaks[,which(del.breaks[2,] > 0),drop=F]

	breaks <- list(gain.breaks, loss.breaks,amp.breaks, del.breaks)
	
	names(breaks) <- c("gain.breaks", "loss.breaks", "amp.breaks", "del.breaks")
	filename <- paste(project, "GLAD","counts", "xls", sep=".")
	spitGALTables(cgh=cgh, gain.count=gain.count, loss.count=loss.count, amp.count=amp.count, del.count=del.count, contig=contig, breaks=breaks, file=filename, acgh.cnv.table=acgh.cnv.table, full.cnv.table=full.cnv.table, genes.table=genes.table, cytoband.table=cytoband.table, mirnas.table=mirnas.table)
}



spitGALTables <- function(cgh, gain.count, loss.count=NULL, amp.count=NULL, del.count=NULL, contig=NULL, breaks, file, acgh.cnv.table="acgh.cnvs.txt", full.cnv.table="full.cnvs.txt", genes.table="all.genes.txt", cytoband.table="cytobands.txt", mirnas.table="mirnas.txt", probeID=featureNames(cgh))
#	function to spit out start, stop probes and positions from thresholded data
#	and to assign genes in the region and anything else you can think of
{
	require(Biobase)
	cgh.name <- deparse(match.call()$cgh)

	genes.table <- read.delim(genes.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	mirnas.table <- read.delim(mirnas.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	cytobands.table <- read.delim(cytoband.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	acgh.cnv.table <- read.delim(acgh.cnv.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	full.cnv.table <- read.delim(full.cnv.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	
	if(length(cgh@thresholds) > 0)
	{
		gainthresh <- cgh@thresholds$gainthresh
		lossthresh <- cgh@thresholds$lossthresh
		ampthresh <- cgh@thresholds$ampthresh
		delthresh <- cgh@thresholds$delthresh
	}
	if(length(cgh@MADS) > 0)
	{
		gainthresh <- paste(cgh@MADS$gainMADS, "MADS")
		lossthresh <- paste(cgh@MADS$lossMADS, "MADS")
		ampthresh <- paste(cgh@MADS$ampMADS, "MADS")
		delthresh <- paste(cgh@MADS$delMADS, "MADS")
	}

	if(ncol(cgh) >1)
		{
		cat(file=file,c("CGH object\t",cgh.name, "\n"), append=F, sep="")
		cat(file=file,c("Gains in ",gain.count," or more cases\n"), append=T, sep="")
		cat(file=file,"Gain threshold = ",gainthresh,"\n", append=T, sep="")
		cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","probes","length MB", "maxM","max.overlap","cases", "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
		}
	if(ncol(cgh) ==1)
		{
		cat(file=file,c(sampleNames(cgh)[1], "\n"), append=F, sep="")
		cat(file=file,c("Gains\n"), append=T, sep="")
		cat(file=file,"Gain threshold = ",gainthresh,"\n", append=T, sep="")
		cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","probes","length MB", "maxM", "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
		}
	
	
	if(ncol(breaks$gain.breaks) > 0) for (i in 1:ncol(breaks$gain.breaks))
		{
		end.index <- breaks$gain.breaks[1,i]
		start.index <- breaks$gain.breaks[1,i]-(breaks$gain.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]	
		iprobes <- breaks$gain.breaks[3,i]
		
		start <- fData(cgh)$start[start.index]	
		end <- fData(cgh)$end[end.index]
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.probe <- probeID[start.index]
		end.probe <- probeID[end.index]
		iM <- exprs(cgh)[start.index:end.index,]
		
		maxM <- max(iM[!is.na(iM)])
		if(ncol(cgh) > 1)
		{
		istates <- apply(assayData(cgh)$GL[start.index:end.index,] >= 1,1,sum)
		icases <- max(istates)

		GL <- assayData(cgh)$GL[start.index:end.index,,drop=F]
		iGL <- matrix(0, nrow(GL), ncol(GL))
		iGL[GL >= 1] <- 1
		gains <- apply(iGL,2,max)
		gain.cases <- sampleNames(cgh)[gains==1]
		}
		
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end >= start),]
		igenes <- igenes.max[which(igenes.max$start <= end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end >= start),]
		imirnas <- imirnas.max[which(imirnas.max$start <= end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end >= start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start <= end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end >= start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start <= end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end >= start),]
		ibands <- ibands.max[which(ibands.max$start <= end),]
		ibands <- ibands$cytoband
		
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.probe, end.probe,iprobes, ilength,maxM)
		if(ncol(cgh)>1) row.result <- c(row.result,icases)
		
		cat(file=file, row.result, append=T, sep="\t")
		cat(file=file, "\t",append=T, sep="")
		if(ncol(cgh)>1)
		{cat(file=file, gain.cases,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")}
		cat(file=file, igenes,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, imirnas,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.acgh.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.full.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, ibands,append=T, sep="-")
		cat(file=file, "\n",append=T, sep="\t")

		}
	
	
	if(ncol(cgh) >1)
		{
			cat(file=file,c("\n\nLosses in ",loss.count," or more cases\n"), append=T, sep="")
			cat(file=file,"Loss threshold = ",lossthresh,"\n", append=T, sep="")
			cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","probes","length MB", "minM","max.overlap","cases", "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")}
	if(ncol(cgh) ==1)
		{
			cat(file=file,c("\n\nLosses\n"), append=T, sep="")
			cat(file=file,"Loss threshold = ",lossthresh,"\n", append=T, sep="")
			cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","probes","length MB", "minM","genes","mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
			}
	
	
	
	if(ncol(breaks$loss.breaks) > 0) for (i in 1:ncol(breaks$loss.breaks))
		{
		end.index <- breaks$loss.breaks[1,i]
		start.index <- breaks$loss.breaks[1,i]-(breaks$loss.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]	
		iprobes <- breaks$loss.breaks[3,i]
		
		start <- fData(cgh)$start[start.index]
		end <- fData(cgh)$end[end.index]
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.probe <- probeID[start.index]
		end.probe <- probeID[end.index]
		iM <- exprs(cgh)[start.index:end.index,]
		
		minM <- min(iM[!is.na(iM)])
		if(ncol(cgh)>1)
		{
		istates <- apply(assayData(cgh)$GL[start.index:end.index,] <= -1,1,sum)
		icases <- max(istates)
		
		GL <- assayData(cgh)$GL[start.index:end.index,,drop=F]
		iGL <- matrix(0, nrow(GL), ncol(GL))
		iGL[GL <= -1] <- 1
		losses <- apply(iGL,2,max)
		loss.cases <- sampleNames(cgh)[losses==1]
		}		
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end >= start),]
		igenes <- igenes.max[which(igenes.max$start <= end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end >= start),]
		imirnas <- imirnas.max[which(imirnas.max$start <= end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end >= start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start <= end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end >= start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start <= end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end >= start),]
		ibands <- ibands.max[which(ibands.max$start <= end),]
		ibands <- ibands$cytoband
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.probe, end.probe, iprobes, ilength,minM)
		if(ncol(cgh)>1) row.result <- c(row.result,icases)
		
		cat(file=file, row.result, append=T, sep="\t")
		cat(file=file, "\t",append=T, sep="")
		if(ncol(cgh)>1)
		{cat(file=file, loss.cases,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")}
		cat(file=file, igenes,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, imirnas,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.acgh.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.full.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, ibands,append=T, sep="-")
		cat(file=file, "\n",append=T, sep="\t")

		}
	
	if(ncol(cgh)>1)
	{
		cat(file=file,c("\n\nAmplifications in ",amp.count," or more cases\n"), append=T, sep="")
		cat(file=file,"Amplification threshold = ",ampthresh,"\n", append=T, sep="")
		cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","probes","length MB", "maxM","max.overlap","cases", "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
		}
	if(ncol(cgh)==1)
	{
		cat(file=file,c("\n\nAmplifications\n"), append=T, sep="")
		cat(file=file,"Amplification threshold = ",ampthresh,"\n", append=T, sep="")
		cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","probes","length MB", "maxM","genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
		}
	
	
	
	if(ncol(breaks$amp.breaks) > 0) for (i in 1:ncol(breaks$amp.breaks))
		{
		end.index <- breaks$amp.breaks[1,i]
		start.index <- breaks$amp.breaks[1,i]-(breaks$amp.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]	
		iprobes <- breaks$amp.breaks[3,i]
		
		start <- fData(cgh)$start[start.index]
		end <- fData(cgh)$end[end.index]
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.probe <- probeID[start.index]
		end.probe <- probeID[end.index]
		iM <- exprs(cgh)[start.index:end.index,]
		
		maxM <- max(iM[!is.na(iM)])
		if(ncol(cgh)>1)
		{
		istates <- apply(assayData(cgh)$GL[start.index:end.index,] > 1,1,sum)
		icases <- max(istates)
		
		GL <- assayData(cgh)$GL[start.index:end.index,,drop=F]
		iGL <- matrix(0, nrow(GL), ncol(GL))
		iGL[GL == 2] <- 1
		amps <- apply(iGL,2,max)
		amp.cases <- sampleNames(cgh)[amps==1]
		}
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end >= start),]
		igenes <- igenes.max[which(igenes.max$start <= end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end >= start),]
		imirnas <- imirnas.max[which(imirnas.max$start <= end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end >= start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start <= end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end >= start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start <= end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end >= start),]
		ibands <- ibands.max[which(ibands.max$start <= end),]
		ibands <- ibands$cytoband
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.probe, end.probe, iprobes, ilength,maxM)
		if(ncol(cgh)>1) row.result <- c(row.result,icases)
		
		cat(file=file, row.result, append=T, sep="\t")
		cat(file=file, "\t",append=T, sep="")
		if(ncol(cgh)>1) 
		{cat(file=file, amp.cases,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")}
		cat(file=file, igenes,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, imirnas,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.acgh.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.full.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, ibands,append=T, sep="-")
		cat(file=file, "\n",append=T, sep="\t")

		}
	
	if(del.count >0)
	{
		cat(file=file,c("\n\nDeletions in ",del.count," or more cases\n"), append=T, sep="")
		cat(file=file,"Deletion threshold = ",delthresh, "\n",append=T, sep="")
		cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","probes","length MB", "minM","max.overlap","cases", "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
		}
	if(del.count ==0)
	{
		cat(file=file,c("\n\nDeletions\n"), append=T, sep="")
		cat(file=file,"Deletion threshold = ",delthresh, "\n",append=T, sep="")
		cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","probes","length MB", "minM","genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
		}
	
	
	
	if(ncol(breaks$del.breaks) > 0) for (i in 1:ncol(breaks$del.breaks))
		{
		end.index <- breaks$del.breaks[1,i]
		start.index <- breaks$del.breaks[1,i]-(breaks$del.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]	
		iprobes <- breaks$del.breaks[3,i]
		
		start <- fData(cgh)$start[start.index]
		end <- fData(cgh)$end[end.index]
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.probe <- probeID[start.index]
		end.probe <- probeID[end.index]
		iM <- exprs(cgh)[start.index:end.index,]
		minM <- min(iM[!is.na(iM)])
		if(ncol(cgh)>1)
		{
		istates <- apply(assayData(cgh)$GL[start.index:end.index,] < -1,1,sum)
		icases <- max(istates)
				
		GL <- assayData(cgh)$GL[start.index:end.index,,drop=F]
		iGL <- matrix(0, nrow(GL), ncol(GL))
		iGL[GL == -2] <- 1
		dels <- apply(iGL,2,max)
		del.cases <- sampleNames(cgh)[dels==1]
		}		
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end >= start),]
		igenes <- igenes.max[which(igenes.max$start <= end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end >= start),]
		imirnas <- imirnas.max[which(imirnas.max$start <= end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end >= start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start <= end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end >= start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start <= end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end >= start),]
		ibands <- ibands.max[which(ibands.max$start <= end),]
		ibands <- ibands$cytoband
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.probe, end.probe, iprobes, ilength, minM)
		if(ncol(cgh)>1) row.result <- c(row.result,icases)
		
		cat(file=file, row.result, append=T, sep="\t")
		cat(file=file, "\t",append=T, sep="")
		if(ncol(cgh)>1) 
		{cat(file=file, del.cases,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")}
		cat(file=file, igenes,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, imirnas,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.acgh.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.full.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, ibands,append=T, sep="-")
		cat(file=file, "\n",append=T, sep="\t")

		}
}



spitGLADTable <-function(cgh, filename, pheno=NULL)
{
	require(Biobase)
	if(!is.null(pheno)){
		
	all.GLAD <- countGLAD(cgh)	
	names(all.GLAD) <- paste("Total", names(all.GLAD), sep=".")
	
	cgh <- cgh[,!is.na(pheno)]
	pheno <- factor(pheno[!is.na(pheno)])
	
	cgh.list <- GLAD.list <- list(NULL)
	for(i in 1:length(unique(pheno))){
		cgh.list[[i]] <- cgh[,which(pheno==unique(pheno)[i])]
		}
	for(i in 1:length(unique(pheno))){
	GLAD.list[[i]] <- countGLAD(cgh.list[[i]])
	names(GLAD.list[[i]]) <- paste(unique(pheno)[i], names(GLAD.list[[i]]), sep=".")
	}
	

	GLAD.table <- do.call(cbind, GLAD.list)
	GLAD.table <- cbind(all.GLAD, GLAD.table)
	}else{
		
		GLAD.table <- countGLAD(cgh)
		
		}
	cat("Writing GLAD table...")
	
	GLAD.table <- cbind.data.frame(fData(cgh), GLAD.table)

	write.table(GLAD.table, file=filename, sep="\t",row.names=F, na="")
	cat(paste(nrow(cgh), "probes\n"))
	cat("Done\n")
}


countGLAD <- function(cgh){
	
	require(Biobase)
	if(!is.element("GL", assayDataElementNames(cgh))) stop("No GL assayData - call aCGH states")
	
	gains <- apply(assayData(cgh)$GL, 1, function(x) sum(x >= 1))
	if(sum(gains) == 0){gain.cases <- rep("", nrow(cgh))}else{
	gain.cases <- apply(assayData(cgh)$GL, 1, function(x) sampleNames(cgh)[which(x >= 1)])
	gain.cases <- as.vector(unlist((lapply(gain.cases, function(x) paste(unlist(x), collapse=", ")))))
	}
	
	losses <- apply(assayData(cgh)$GL, 1, function(x) sum(x <= -1))
	if(sum(losses) == 0){loss.cases <- rep("", nrow(cgh))}else{
	loss.cases <- apply(assayData(cgh)$GL, 1, function(x) sampleNames(cgh)[which(x <= -1)])
	loss.cases <- as.vector(unlist((lapply(loss.cases, function(x) paste(unlist(x), collapse=", ")))))
	}
	
	amps <- apply(assayData(cgh)$GL, 1, function(x) sum(x == 2))
	if(sum(amps) == 0){amp.cases <- rep("", nrow(cgh))}else{
	amp.cases <- apply(assayData(cgh)$GL, 1, function(x) sampleNames(cgh)[which(x == 2)])
	amp.cases <- as.vector(unlist((lapply(amp.cases, function(x) paste(unlist(x), collapse=", ")))))	}
	
	dels <- apply(assayData(cgh)$GL, 1, function(x) sum(x == -2))
	if(sum(dels) == 0){del.cases <- rep("", nrow(cgh))}else{
	del.cases <- apply(assayData(cgh)$GL, 1, function(x) sampleNames(cgh)[which(x == -2)])
	del.cases <- as.vector(unlist((lapply(del.cases, function(x) paste(unlist(x), collapse=", ")))))
	}
	
	GLAD <- cbind.data.frame(gains, losses, amps,  dels, gain.cases, loss.cases, amp.cases, del.cases)
	GLAD
}


fisherTestCGH <- function(cgh,pheno, project=NULL, p.val.adjustment="fdr", pval.thresh=0.05, LOOCV=F, return.eset=T)
{
	require(Biobase)
	options(warn=-1)
	
	if(is.null(project)) stop("Please spcefiy \"project\"")
	cgh <- cgh[,!is.na(pheno)]
	pheno <- factor(pheno[!is.na(pheno)])
	
	phenos <- unique(as.character(pheno))
	cat("Performing Fishers Exact tests for gains, losses, amps and dels in",
		levels(pheno)[1:(length(levels(pheno))-1)],"and",levels(pheno)[length(levels(pheno))],"groups\n")
	
	gain <- (assayData(cgh)$GL >= 1)
	loss <-  (assayData(cgh)$GL <= -1)
	amp <- (assayData(cgh)$GL >= 2)
	del <-  (assayData(cgh)$GL <= -2)	
	
	pval.gain <- pval.loss <- pval.amp <- pval.del <- rep(NA, nrow(cgh))
	
	
	if(identical(LOOCV,T)){
		LOOCV.p.gain <- LOOCV.p.loss <- LOOCV.p.amp <- LOOCV.p.del <- matrix(NA, nrow(cgh), ncol(cgh))
		LOOCV.adjp.gain <- LOOCV.adjp.loss <- LOOCV.adjp.amp <- LOOCV.adjp.del <- matrix(NA, nrow(cgh), ncol(cgh))}
		
	gain.count <- loss.count <- amp.count <- del.count <- matrix(0,nrow=nrow(cgh), ncol=length(phenos))
		
	for(i in 1:nrow(cgh))
	{
	
	if(ncol(cgh) > sum(gain[i,]) && sum(gain[i,]) > 0) pval.gain[i] <- fisher.test(table(gain[i,], pheno))$p.value
	if(ncol(cgh) > sum(loss[i,]) && sum(loss[i,]) > 0) pval.loss[i] <-  fisher.test(table(loss[i,], pheno))$p.value
	if(ncol(cgh) > sum(amp[i,]) && sum(amp[i,]) > 0) pval.amp[i] <-  fisher.test(table(amp[i,], pheno))$p.value
	if(ncol(cgh) > sum(del[i,]) && sum(del[i,]) > 0) pval.del[i] <-  fisher.test(table(del[i,], pheno))$p.value
		
	g.count <- l.count <- a.count <- d.count <- rep(0, length(phenos))
				
	for (p in 1:length(phenos))
		{
		g.count[p] <- sum(assayData(cgh)$GL[i,phenos[p]==pheno]>= 1)
		l.count[p] <- sum(assayData(cgh)$GL[i,phenos[p]==pheno]<= -1)
		a.count[p] <- sum(assayData(cgh)$GL[i,phenos[p]==pheno]== 2)
		d.count[p] <- sum(assayData(cgh)$GL[i,phenos[p]==pheno]== -2)
		}
		
	gain.count[i,] <- g.count
	loss.count[i,] <- l.count
	amp.count[i,] <- a.count
	del.count[i,] <- d.count
		
	if(i %% 100 == 0) cat(i," ")	
	}
	
	if(identical(LOOCV,T)){
		for(j in 1:ncol(assayData(cgh)$GL))
			{
			cat("\nLeave one out cross validation - sample",j,"\n")
			GL.j <- assayData(cgh)$GL[,-j]
			gain.j <- (GL.j >= 1)
			loss.j <-  (GL.j <= -1)
			amp.j <- (GL.j >= 2)
			del.j <-  (GL.j <= -2)
			
			for(i in 1:nrow(assayData(cgh)$GL))
				{	
				if(ncol(GL.j) > sum(gain.j[i,]) && sum(gain.j[i,]) > 0) LOOCV.p.gain[i,j] <- fisher.test(gain.j[i,], pheno[-j])$p.value
				if(ncol(GL.j) > sum(loss.j[i,]) && sum(loss.j[i,]) > 0) LOOCV.p.loss[i,j] <-  fisher.test(table(loss.j[i,], pheno[-j]))$p.value
				if(ncol(GL.j) > sum(amp.j[i,]) && sum(amp.j[i,]) > 0) LOOCV.p.amp[i,j] <-  fisher.test(table(amp.j[i,], pheno[-j]))$p.value
				if(ncol(GL.j) > sum(del.j[i,]) && sum(del.j[i,]) > 0) LOOCV.p.del[i,j] <-  fisher.test(table(del.j[i,], pheno[-j]))$p.value
			
				if(i %% 100 == 0) cat(i," ")
				}
			}
		}
			
	fData(cgh)$Fishers.p.gain <- pval.gain
	fData(cgh)$Fishers.adjp.gain <- p.adjust.2.01(pval.gain,p.val.adjustment)
	fData(cgh)$Fishers.p.loss <- pval.loss
	fData(cgh)$Fishers.adjp.loss <- p.adjust.2.01(pval.loss,p.val.adjustment)
	fData(cgh)$Fishers.p.amp <- pval.amp
	fData(cgh)$Fishers.adjp.amp <- p.adjust.2.01(pval.amp,p.val.adjustment)
	fData(cgh)$Fishers.p.del <- pval.del
	fData(cgh)$Fishers.adjp.del <- p.adjust.2.01(pval.del,p.val.adjustment)
	
	if(identical(LOOCV,T)){
		
		for(i in 1:ncol(assayData(cgh)$GL))
		{
			LOOCV.adjp.gain[,i] <- p.adjust.2.01(LOOCV.p.gain[,i],p.val.adjustment)
			LOOCV.adjp.loss[,i] <- p.adjust.2.01(LOOCV.p.loss[,i],p.val.adjustment)
			LOOCV.adjp.amp[,i] <- p.adjust.2.01(LOOCV.p.amp[,i],p.val.adjustment)
			LOOCV.adjp.del[,i] <- p.adjust.2.01(LOOCV.p.del[,i],p.val.adjustment)
		}
	
		fData(cgh)$LOOCV.sig.p.gain.count <- apply(LOOCV.p.gain,1, function(x) sum(x<pval.thresh, na.rm=T))
		fData(cgh)$LOOCV.sig.adjp.gain.count <- apply(LOOCV.adjp.gain,1, function(x) sum(x<pval.thresh, na.rm=T))
		fData(cgh)$LOOCV.sig.p.loss.count <- apply(LOOCV.p.loss,1, function(x) sum(x<pval.thresh, na.rm=T))
		fData(cgh)$LOOCV.sig.adjp.loss.count <- apply(LOOCV.adjp.loss,1, function(x) sum(x<pval.thresh, na.rm=T))
		fData(cgh)$LOOCV.sig.p.amp.count <- apply(LOOCV.p.amp,1, function(x) sum(x<pval.thresh, na.rm=T))
		fData(cgh)$LOOCV.sig.adjp.amp.count <- apply(LOOCV.adjp.amp,1, function(x) sum(x<pval.thresh, na.rm=T))
		fData(cgh)$LOOCV.sig.p.del.count <- apply(LOOCV.p.del,1, function(x) sum(x<pval.thresh, na.rm=T))
		fData(cgh)$LOOCV.sig.adjp.del.count <- apply(LOOCV.adjp.del,1, function(x) sum(x<pval.thresh, na.rm=T))
		
		fData(cgh)$LOOCV.max.p.gain <- apply(LOOCV.p.gain,1, function(x) max(x, na.rm=T))
		fData(cgh)$LOOCV.max.adjp.gain <- apply(LOOCV.adjp.gain,1, function(x) max(x, na.rm=T))
		fData(cgh)$LOOCV.max.p.loss <- apply(LOOCV.p.loss,1, function(x) max(x, na.rm=T))
		fData(cgh)$LOOCV.max.adjp.loss <- apply(LOOCV.adjp.loss,1, function(x) max(x, na.rm=T))
		fData(cgh)$LOOCV.max.p.amp <- apply(LOOCV.p.amp,1, function(x) max(x, na.rm=T))
		fData(cgh)$LOOCV.max.adjp.amp <- apply(LOOCV.adjp.amp,1, function(x) max(x, na.rm=T))
		fData(cgh)$LOOCV.max.p.del <- apply(LOOCV.p.del,1, function(x) max(x, na.rm=T))
		fData(cgh)$LOOCV.max.adjp.del <- apply(LOOCV.adjp.del,1, function(x) max(x, na.rm=T))
		
		fData(cgh)$LOOCV.max.p.gain[which(is.infinite(fData(cgh)$LOOCV.max.p.gain))]<- NA
		fData(cgh)$LOOCV.max.adjp.gain[which(is.infinite(fData(cgh)$LOOCV.max.adjp.gain))]<- NA
		fData(cgh)$LOOCV.max.p.loss[which(is.infinite(fData(cgh)$LOOCV.max.p.loss))]<- NA
		fData(cgh)$LOOCV.max.adjp.loss[which(is.infinite(fData(cgh)$LOOCV.max.adjp.loss))]<- NA
		fData(cgh)$LOOCV.max.p.amp[which(is.infinite(fData(cgh)$LOOCV.max.p.amp))]<- NA
		fData(cgh)$LOOCV.max.adjp.amp[which(is.infinite(fData(cgh)$LOOCV.max.adjp.amp))]<- NA
		fData(cgh)$LOOCV.max.p.del[which(is.infinite(fData(cgh)$LOOCV.max.p.del))]<- NA
		fData(cgh)$LOOCV.max.adjp.del[which(is.infinite(fData(cgh)$LOOCV.max.adjp.del))]<- NA
	}
		
	cgh$Fishers.pheno <- pheno
	
	cat("\nDone\n")

	gain.count.d <- data.frame(gain.count)
	names(gain.count.d) <- paste("gains",phenos,sep=".")
	loss.count.d <- data.frame(loss.count)
	names(loss.count.d) <- paste("losses",phenos,sep=".")
	amp.count.d <- data.frame(amp.count)
	names(amp.count.d) <- paste("amps",phenos,sep=".")
	del.count.d <- data.frame(del.count)
	names(del.count.d) <- paste("dels",phenos,sep=".")
	
	pheno.counts <- cbind(gain.count.d, loss.count.d, amp.count.d, del.count.d)
	
	glad <- assayData(cgh)$GLAD
	gains.all <- rowSums(gain)
	losses.all <- rowSums(loss)
	amps.all <- rowSums(amp)
	dels.all <- rowSums(del)
	
	cgh.name <- deparse(match.call()$cgh)
	
	if(!is.element("Fishers.results", list.files())) dir.create("Fishers.results")
	Fisher.out <- cbind(fData(cgh), gains.all, losses.all, amps.all, dels.all, pheno.counts)
	all.out <- paste("Fishers.results/",project, ".Fishers.results.all.probes.xls", sep="")
	cat("Results of Fishers Exact tests for gains, losses, amps and dels in",
		levels(pheno)[1:(length(levels(pheno))-1)],"and",levels(pheno)[length(levels(pheno))],"groups\n",file=all.out, append=F)
	cat("CGH object\t", cgh.name, "\n","No. of samples\t",
		ncol(cgh),"\n","Probes\t", nrow(cgh),"\n","p value adjustment\t",p.val.adjustment, "\n\n",file=all.out, sep="", append=T)
	write.table(Fisher.out, file=all.out, sep="\t",row.names=F, append=T)
	
	file.out <- paste("Fishers.results/",project, ".Fishers.results.xls", sep="")
	cat("Results of Fishers Exact tests for gains, losses, amps and dels in",
		levels(pheno)[1:(length(levels(pheno))-1)],"and",levels(pheno)[length(levels(pheno))],"groups\n",file=file.out, append=F)
	cat("CGH object\t", cgh.name, "\n","No. of samples\t",
		ncol(cgh),"\n","Probes\t", nrow(cgh),"\n","p value adjustment\t",p.val.adjustment, "\n\n",file=file.out, sep="", append=T)
	
	if(!is.null(p.val.adjustment)){
		cat("\nProbes differentially gained in",
		levels(pheno)[1:(length(levels(pheno))-1)],"and",levels(pheno)[length(levels(pheno))],"groups\n",file=file.out, append=T)
		cat("Gain threshold\t", cgh$thresholds$gainthresh, "\n",file=file.out, append=T)
		gains.out <- Fisher.out[which(fData(cgh)$Fishers.adjp.gain < pval.thresh),]
		write.table(gains.out, file=file.out, sep="\t",row.names=F, append=T)
		
		cat("\nProbes differentially lost in",
		levels(pheno)[1:(length(levels(pheno))-1)],"and",levels(pheno)[length(levels(pheno))],"groups\n",file=file.out, append=T)
		cat("Loss threshold\t", cgh$thresholds$lossthresh, "\n",file=file.out, append=T)
		losses.out <- Fisher.out[which(fData(cgh)$Fishers.adjp.loss < pval.thresh),]
		write.table(losses.out, file=file.out, sep="\t",row.names=F, append=T)
		
		cat("\nProbes differentially amplified in",
		levels(pheno)[1:(length(levels(pheno))-1)],"and",levels(pheno)[length(levels(pheno))],"groups\n",file=file.out, append=T)
		cat("Amp threshold\t", cgh$thresholds$ampthresh, "\n",file=file.out, append=T)
		amps.out <- Fisher.out[which(fData(cgh)$Fishers.adjp.amp < pval.thresh),]
		write.table(amps.out, file=file.out, sep="\t",row.names=F, append=T)
		
		cat("\nProbes differentially deleted in",
		levels(pheno)[1:(length(levels(pheno))-1)],"and",levels(pheno)[length(levels(pheno))],"groups\n",file=file.out, append=T)
		cat("Del threshold\t", cgh$thresholds$delthresh, "\n",file=file.out, append=T)
		dels.out <- Fisher.out[which(fData(cgh)$Fishers.adjp.del < pval.thresh),]
		write.table(dels.out, file=file.out, sep="\t",row.names=F, append=T)
	}
	
	if(identical(return.eset,T)) return(cgh)
}


listBreaksFisher <- function(cgh, contig=2, project, adjusted.p=T, pval.thresh=0.05, acgh.cnv.table="acgh.cnvs.txt", full.cnv.table="full.cnvs.txt", genes.table="all.genes.txt", cytoband.table="cytobands.txt", mirnas.table="mirnas.txt", probeID=fData(cgh)$bac.id)
{
	
	require(Biobase)
	if(adjusted.p == F) pvals <- cbind(fData(cgh)$Fishers.p.gain,fData(cgh)$Fishers.p.loss,fData(cgh)$Fishers.p.amp,fData(cgh)$Fishers.p.del)
	if(adjusted.p == T) pvals <- cbind(fData(cgh)$Fishers.adjp.gain,fData(cgh)$Fishers.adjp.loss,fData(cgh)$Fishers.adjp.amp,fData(cgh)$Fishers.adjp.del)
	
	pvals[which(pvals[,1] < pval.thresh),1] <- 2
	pvals[which(pvals[,2] < pval.thresh),2] <- 2
	pvals[which(pvals[,3] < pval.thresh),3] <- 2
	pvals[which(pvals[,4] < pval.thresh),4] <- 2
	
	gain.breaks <- listBreaks(cgh=cgh, GL.column=pvals[,1], contig=contig)
	gain.breaks <- gain.breaks[,which(gain.breaks[2,] == 2), drop=F]
	loss.breaks <- listBreaks(cgh=cgh, GL.column=pvals[,2], contig=contig)
	loss.breaks <- loss.breaks[,which(loss.breaks[2,] == 2) ,drop=F]
	amp.breaks <- listBreaks(cgh=cgh, GL.column=pvals[,3], contig=contig)
	amp.breaks <- amp.breaks[,which(amp.breaks[2,] == 2),drop=F]
	del.breaks <- listBreaks(cgh=cgh, GL.column=pvals[,4], contig=contig)
	del.breaks <- del.breaks[,which(del.breaks[2,] == 2),drop=F]

	breaks <- list(gain.breaks, loss.breaks,amp.breaks, del.breaks)
	
	names(breaks) <- c("gain.breaks", "loss.breaks", "amp.breaks", "del.breaks")
	
	if(!is.element("Fishers.results", list.files())) dir.create("Fishers.results")
	
	filename <- paste("Fishers.results/",project, ".Fishers.regions.xls", sep="")

	genes.table <- read.delim(genes.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	mirnas.table <- read.delim(mirnas.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	cytobands.table <- read.delim(cytoband.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	acgh.cnv.table <- read.delim(acgh.cnv.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	full.cnv.table <- read.delim(full.cnv.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	
	if(!is.null(cgh@thresholds))
	{
		gainthresh <- cgh$thresholds$gainthresh
		lossthresh <- cgh$thresholds$lossthresh
		ampthresh <- cgh$thresholds$ampthresh
		delthresh <- cgh$thresholds$delthresh
	}
#	if(!is.null(cgh@MADS))
#	{
#		gainthresh <- paste(exprs(cgh)ADS$gainMADS, "MADS")
#		lossthresh <- paste(exprs(cgh)ADS$lossMADS, "MADS")
#		ampthresh <- paste(exprs(cgh)ADS$ampMADS, "MADS")
#		delthresh <- paste(exprs(cgh)ADS$delMADS, "MADS")
#	}

	cgh.name <- deparse(match.call()$cgh)

	cat("Results of Fishers Exact tests for gains, losses, amps and dels in",
		paste(levels(cgh$Fishers.pheno), collapse=" and "),"groups\n",file=filename, append=F)
		
	cat("CGH object\t", cgh.name, "\n","No. of samples\t",
		ncol(cgh),"\n","Probes\t", nrow(cgh),"\n","p value adjustment\t",adjusted.p, "\n\n",file=filename, sep="", append=T)

	cat("\n\nRegions differentially gained in",
		paste(levels(cgh$Fishers.pheno), collapse=" and "),"groups\n",file=filename, append=T)
		cat("Gain threshold\t", cgh@thresholds$gainthresh, "\n",file=filename, append=T)

		cat(file=filename,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","BACs","length MB", "maxM","max.overlap","gain cases", paste("gains in",levels(cgh$Fishers.pheno)), "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
	
	GLAD <- matrix(0,nrow(cgh),4)	
	for(i in 1:nrow(cgh))
	{
		GLAD[i,1] <- length(assayData(cgh)$GL[i,which(assayData(cgh)$GL[i,] > 0)])
		GLAD[i,2] <- length(assayData(cgh)$GL[i,which(assayData(cgh)$GL[i,] < 0)])
		GLAD[i,3] <- length(assayData(cgh)$GL[i,which(assayData(cgh)$GL[i,] > 1)])
		GLAD[i,4] <- length(assayData(cgh)$GL[i,which(assayData(cgh)$GL[i,] < -1)])
	}
	
	GLAD <- as.data.frame(GLAD)
	names(GLAD) <- c("Gains", "Losses", "Amps", "Dels")

	if(ncol(breaks$gain.breaks) > 0) for (i in 1:ncol(breaks$gain.breaks))
		{
		end.index <- breaks$gain.breaks[1,i]
		start.index <- breaks$gain.breaks[1,i]-(breaks$gain.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]	
		iProbes <- breaks$gain.breaks[3,i]
		
		start <- fData(cgh)$start[start.index]	
		end <- fData(cgh)$end[end.index]
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.probe <- probeID[start.index]
		end.probe <- probeID[end.index]
		iM <- exprs(cgh)[start.index:end.index,]
		maxM <- max(iM[!is.na(iM)])

		istates <- GLAD[start.index:end.index,1]
		icases <- max(istates)

		GL <- assayData(cgh)$GL[start.index:end.index,,drop=F]
		iGL <- matrix(0, nrow(GL), ncol(GL))
		iGL[GL >= 1] <- 1
		gains <- apply(iGL,2,max)
		gain.cases <- sampleNames(cgh)[gains==1]
		pheno.gains <- unlist(lapply(split(gains, cgh$Fishers.pheno), sum))

		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end > start),]
		igenes <- igenes.max[which(igenes.max$start < end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end > start),]
		imirnas <- imirnas.max[which(imirnas.max$start < end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end > start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start < end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end > start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start < end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end > start),]
		ibands <- ibands.max[which(ibands.max$start < end),]
		ibands <- ibands$cytoband
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.probe, end.probe,iProbes, ilength,maxM, icases)
		
		cat(file=filename, row.result, append=T, sep="\t")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, gain.cases,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, pheno.gains,append=T, sep="\t")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, igenes,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, imirnas,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, i.acgh.cnvs,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, i.full.cnvs,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, ibands,append=T, sep="-")
		cat(file=filename, "\n",append=T, sep="\t")

		}
	
	cat("\n\nRegions differentially lost in",
		paste(levels(cgh$Fishers.pheno), collapse=" and "),"groups\n",file=filename, append=T)
		cat("Loss threshold\t", cgh@thresholds$lossthresh, "\n",file=filename, append=T)

		cat(file=filename,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","BACs","length MB", "minM","max.overlap","loss cases", paste("losses in",levels(cgh$Fishers.pheno)), "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")

	
	if(ncol(breaks$loss.breaks) > 0) for (i in 1:ncol(breaks$loss.breaks))
		{
		end.index <- breaks$loss.breaks[1,i]
		start.index <- breaks$loss.breaks[1,i]-(breaks$loss.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]	
		iProbes <- breaks$loss.breaks[3,i]
		
		start <- fData(cgh)$start[start.index]
		end <- fData(cgh)$end[end.index]
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.probe <- probeID[start.index]
		end.probe <- probeID[end.index]
		iM <- exprs(cgh)[start.index:end.index,]
		minM <- min(iM[!is.na(iM)])

		istates <- GLAD[start.index:end.index,2]
		icases <- max(istates)
		
		GL <- assayData(cgh)$GL[start.index:end.index,,drop=F]
		iGL <- matrix(0, nrow(GL), ncol(GL))
		iGL[GL <= -1] <- 1
		losses <- apply(iGL,2,max)
		loss.cases <- sampleNames(cgh)[losses==1]
		pheno.losses <- unlist(lapply(split(losses, cgh$Fishers.pheno), sum))
		
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end > start),]
		igenes <- igenes.max[which(igenes.max$start < end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end > start),]
		imirnas <- imirnas.max[which(imirnas.max$start < end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end > start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start < end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end > start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start < end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end > start),]
		ibands <- ibands.max[which(ibands.max$start < end),]
		ibands <- ibands$cytoband
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.probe, end.probe, iProbes, ilength,minM, icases)
		
		cat(file=filename, row.result, append=T, sep="\t")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, loss.cases,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, pheno.losses,append=T, sep="\t")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, igenes,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, imirnas,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, i.acgh.cnvs,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, i.full.cnvs,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, ibands,append=T, sep="-")
		cat(file=filename, "\n",append=T, sep="\t")

		}
	
	cat("\n\nRegions differentially amplified in",
		paste(levels(cgh$Fishers.pheno), collapse=" and "),"groups\n",file=filename, append=T)
		cat("Amplification threshold\t", cgh$thresholds$ampthresh, "\n",file=filename, append=T)

	cat(file=filename,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","BACs","length MB", "maxM","max.overlap","amplified cases", paste("amplifications in",levels(cgh$Fishers.pheno)), "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
	
	if(ncol(breaks$amp.breaks) > 0) for (i in 1:ncol(breaks$amp.breaks))
		{
		end.index <- breaks$amp.breaks[1,i]
		start.index <- breaks$amp.breaks[1,i]-(breaks$amp.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]	
		iProbes <- breaks$amp.breaks[3,i]
		
		start <- fData(cgh)$start[start.index]
		end <- fData(cgh)$end[end.index]
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.probe <- probeID[start.index]
		end.probe <- probeID[end.index]
		iM <- exprs(cgh)[start.index:end.index,]
		maxM <- max(iM[!is.na(iM)])

		istates <- GLAD[start.index:end.index,3]
		icases <- max(istates)
		
		GL <- assayData(cgh)$GL[start.index:end.index,,drop=F]
		iGL <- matrix(0, nrow(GL), ncol(GL))
		iGL[GL == 2] <- 1
		amps <- apply(iGL,2,max)
		amp.cases <- sampleNames(cgh)[amps==1]
		pheno.amps <- unlist(lapply(split(amps, cgh$Fishers.pheno), sum))
		
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end > start),]
		igenes <- igenes.max[which(igenes.max$start < end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end > start),]
		imirnas <- imirnas.max[which(imirnas.max$start < end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end > start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start < end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end > start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start < end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end > start),]
		ibands <- ibands.max[which(ibands.max$start < end),]
		ibands <- ibands$cytoband
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.probe, end.probe, iProbes, ilength,maxM, icases)
		
		cat(file=filename, row.result, append=T, sep="\t")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, amp.cases,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, pheno.amps,append=T, sep="\t")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, igenes,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, imirnas,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, i.acgh.cnvs,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, i.full.cnvs,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, ibands,append=T, sep="-")
		cat(file=filename, "\n",append=T, sep="\t")

		}
	
	cat("\n\nRegions differentially deleted in",
		paste(levels(cgh$Fishers.pheno), collapse=" and "),"groups\n",file=filename, append=T)
		cat("Deletion threshold\t", cgh$thresholds$delthresh, "\n",file=filename, append=T)

		cat(file=filename,c("chrom", "start", "end","start.MB", "end.MB","start.probe", "end.probe","BACs","length MB", "minM","max.overlap","deleted cases", paste("deletions in",levels(cgh$Fishers.pheno)), "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
	
	if(ncol(breaks$del.breaks) > 0) for (i in 1:ncol(breaks$del.breaks))
		{
		end.index <- breaks$del.breaks[1,i]
		start.index <- breaks$del.breaks[1,i]-(breaks$del.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]	
		iProbes <- breaks$del.breaks[3,i]
		
		start <- fData(cgh)$start[start.index]
		end <- fData(cgh)$end[end.index]
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.probe <- probeID[start.index]
		end.probe <- probeID[end.index]
		iM <- exprs(cgh)[start.index:end.index,]
		minM <- min(iM[!is.na(iM)])

		istates <- GLAD[start.index:end.index,4]
		icases <- max(istates)
	
		GL <- assayData(cgh)$GL[start.index:end.index,,drop=F]
		iGL <- matrix(0, nrow(GL), ncol(GL))
		iGL[GL == -2] <- 1
		dels <- apply(iGL,2,max)
		del.cases <- sampleNames(cgh)[dels==1]
		pheno.dels <- unlist(lapply(split(dels, cgh$Fishers.pheno), sum))		
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end > start),]
		igenes <- igenes.max[which(igenes.max$start < end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end > start),]
		imirnas <- imirnas.max[which(imirnas.max$start < end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end > start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start < end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end > start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start < end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end > start),]
		ibands <- ibands.max[which(ibands.max$start < end),]
		ibands <- ibands$cytoband
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.probe, end.probe, iProbes, ilength, minM, icases)
		
		cat(file=filename, row.result, append=T, sep="\t")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, del.cases,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, pheno.dels,append=T, sep="\t")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, igenes,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, imirnas,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, i.acgh.cnvs,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, i.full.cnvs,append=T, sep=", ")
		cat(file=filename, "\t",append=T, sep="")
		cat(file=filename, ibands,append=T, sep="-")
		cat(file=filename, "\n",append=T, sep="\t")

		}
}



p.adjust.2.01 <- function (p, method = p.adjust.methods, n = length(p)) 
{
    if (n == 1) 
        return(p)
    method <- match.arg(method)
    switch(method, holm = {
        i <- 1:n
        o <- order(p)
        ro <- order(o)
        pmin(1, cummax((n - i + 1) * p[o]))[ro]
    }, hochberg = {
        i <- n:1
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin((n - i + 1) * p[o]))[ro]
    }, hommel = {
        i <- 1:n
        s <- sort(p, index = TRUE)
        p <- s$x
        ro <- order(s$ix)
        q <- pa <- rep.int(min(n * p/(1:n)), n)
        for (j in (n - 1):2) {
            q1 <- min(j * p[(n - j + 2):n]/(2:j))
            q[1:(n - j + 1)] <- pmin(j * p[1:(n - j + 1)], q1)
            q[(n - j + 2):n] <- q[n - j + 1]
            pa <- pmax(pa, q)
        }
        pmax(pa, p)[ro]
    }, fdr = {
        i <- n:1
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(n/i * p[o]))[ro]
    }, bonferroni = pmin(n * p, 1), none = p)
}

latticePlotFishers <- function(cgh, device="PDF", project, plot.adjusted.p=T, LOOCV=F, pval.thresh=0.05, use_genomic_positions=T)
{
require(Biobase)
library(lattice)
library(grid)
library(marray)

	phenos <- unique(cgh$Fishers.pheno)
	GL.groups <-rep(list(NULL),length(phenos))
	if (use_genomic_positions==T) {
		MB <- calculateMB(cgh)
	} else { MB <- 1:length(featureNames(cgh)) }
	
	chlabs <- as.character(c(1:22,"X","Y"))
	chr.nms <- names(table(fData(cgh)$chrom))
	chr.nms[match("23",chr.nms,nomatch=0)] <- "X"
	chr.nms[match("24",chr.nms,nomatch=0)] <- "Y"
	ends <- sapply(split(MB,fData(cgh)$chrom),max)

	begs <- sapply(split(MB,fData(cgh)$chrom),min)
	mids <- (begs+ends)/2

	freq.matrix <- matrix(0,nrow=nrow(cgh), 5)
	GL.groups[[1]]<- assayData(cgh)$GL[,which(cgh$Fishers.pheno==phenos[1])]
	freq.matrix[,1] <- apply(GL.groups[[1]],1,function(x) sum(x >= 1))/ncol(GL.groups[[1]])
	freq.matrix[,2] <- 0-(apply(GL.groups[[1]],1,function(x) sum(x <= -1))/ncol(GL.groups[[1]]))
	freq.matrix[,3] <- apply(GL.groups[[1]],1,function(x) sum(x >= 2))/ncol(GL.groups[[1]])
	freq.matrix[,4] <- 0-(apply(GL.groups[[1]],1,function(x) sum(x <= -2))/ncol(GL.groups[[1]]))
	freq.matrix[,5] <- MB

	for (i in 2:length(phenos))
		{
		freq.matrix.i <- matrix(0,nrow=nrow(assayData(cgh)$GL), 5)
		GL.groups[[i]]<- assayData(cgh)$GL[,which(cgh$Fishers.pheno==phenos[i])]
		freq.matrix.i[,1] <- apply(GL.groups[[i]],1,function(x) sum(x >= 1))/ncol(GL.groups[[i]])
		freq.matrix.i[,2] <- 0-(apply(GL.groups[[i]],1,function(x) sum(x <= -1))/ncol(GL.groups[[i]]))
		freq.matrix.i[,3] <- apply(GL.groups[[i]],1,function(x) sum(x >= 2))/ncol(GL.groups[[i]])
		freq.matrix.i[,4] <- 0-(apply(GL.groups[[i]],1,function(x) sum(x <= -2))/ncol(GL.groups[[i]]))
		freq.matrix.i[,5] <- MB
		freq.matrix <- rbind(freq.matrix,freq.matrix.i)
		}
		
	freq.matrix.phenos <- rep(as.character(phenos[1]),nrow(cgh))
	for (i in 2:length(phenos))
		{
			freq.matrix.phenos <- c(freq.matrix.phenos,rep(as.character(phenos[i]),nrow(cgh)))
		}
	freq.matrix <- data.frame(freq.matrix.phenos,freq.matrix)
	names(freq.matrix) <-c("phenos","gains","losses","amps","dels","MB")
	freq.matrix[freq.matrix == 0] <- NA
	pvals <- matrix(1,nrow(cgh),4)
	
	if(identical(LOOCV,T)){
		
	
	if(identical(plot.adjusted.p,T)){
		pvals[,1] <- fData(cgh)$max.adjp.gain
		pvals[,2] <- fData(cgh)$max.adjp.loss
		pvals[,3] <- fData(cgh)$max.adjp.amp
		pvals[,4] <- fData(cgh)$max.adjp.del
		}
	if(identical(plot.adjusted.p,F)){
		pvals[,1] <- fData(cgh)$max.p.gain
		pvals[,2] <- fData(cgh)$max.p.loss
		pvals[,3] <- fData(cgh)$max.p.amp
		pvals[,4] <- fData(cgh)$max.p.del
		}	
		
		}else{
			
	if(identical(plot.adjusted.p,T)){
		pvals[,1] <- fData(cgh)$Fishers.adjp.gain
		pvals[,2] <- fData(cgh)$Fishers.adjp.loss
		pvals[,3] <- fData(cgh)$Fishers.adjp.amp
		pvals[,4] <- fData(cgh)$Fishers.adjp.del
		}
	if(identical(plot.adjusted.p,F)){
		pvals[,1] <- fData(cgh)$Fishers.p.gain
		pvals[,2] <- fData(cgh)$Fishers.p.loss
		pvals[,3] <- fData(cgh)$Fishers.p.amp
		pvals[,4] <- fData(cgh)$Fishers.p.del
		}
			
			
		}
	
	pvals[is.na(pvals)] <- 1
	pvals[is.infinite(pvals)] <- 1
	log.pvals <- -(logb(pvals,10))
	
	p.thresh <- -(logb(pval.thresh,10))/max(log.pvals)
	log.pvals <- log.pvals/max(log.pvals)
	
	log.pvals[log.pvals < p.thresh] <- NA
	log.pvals[,2] <- -(log.pvals[,2])
	log.pvals[,4] <- -(log.pvals[,4])
	
	if(identical(LOOCV,T)){
		Fisher.p.values <- rep("LOOCV.max.Fisher.p.values",nrow(cgh))
		}else{
		Fisher.p.values <- rep("log.Fisher.p.values",nrow(cgh))
		}
		
	pval.matrix <- data.frame(Fisher.p.values,log.pvals,MB)
	names(pval.matrix) <-c("phenos","gains","losses","amps","dels","MB")
		
	fisher.matrix <- rbind(pval.matrix,freq.matrix)
	
	if(identical(LOOCV,T)){
			gain.title <- paste(project,"- LOOCV Fishers tests of Gains and Losses")
			amp.title <- paste(project,"- LOOCV Fishers tests of Amplifications and Deletions")
			}else{
			gain.title <- paste(project,"- Fishers tests of Gains and Losses")
			amp.title <- paste(project,"- Fishers tests of Amplifications and Deletions")
			
			}
	
	r.pal <- maPalette(low="white", high="red", k= length(cgh$Fishers.pheno))
	g.pal <- maPalette(low="white", high="green", k= length(cgh$Fishers.pheno))
	
	if(!is.element("Fishers.plots", list.files())) dir.create("Fishers.plots")
	
	GL.plot <- xyplot(fisher.matrix$gains+fisher.matrix$losses~fisher.matrix$MB|fisher.matrix$phenos, type="h", layout=c(1,length(phenos)+1), col=c("dark green","red"), main=gain.title, xlab="Chromosome", ylab="Frequency",
	scales=list(x=list(at=mids,labels=chr.nms, cex=0.5)), lwd=3,  ylim=c(-1.1,1.1),
	panel = function(...)
	{panel.abline(h = 0, lty = 1);panel.abline(v=c(1,ends),col="grey");panel.xyplot(...)})
	
	if(device == "quartz") {quartz(title=paste(project,gain.title),width=6, height=10)}
	if(device == "PNG")
		{GL.file <- paste(project,"GL.Fishers","png", sep=".")
		trellis.device("png", color=T, file=paste("Fishers.plots/",GL.file,sep=""),width=1600, height=800, pointsize=12)}
	if(device == "JPEG")
		{GL.file <- paste(project,"GL.Fishers","jpeg", sep=".")
		trellis.device("jpeg", color=T, file=paste("Fishers.plots/",GL.file,sep=""),width=1600, height=800, pointsize=12)}
	if(device == "PDF")
		{GL.file <- paste(project,"GL.Fishers","pdf", sep=".")
		trellis.device("pdf", color=T, file=paste("Fishers.plots/",GL.file,sep=""),width=8, height=6, pointsize=8)}
	if(device == "PS")
		{GL.file <- paste(project,"GL.Fishers","ps", sep=".")
		trellis.device("postscript", color=T, file=paste("Fishers.plots/",GL.file,sep=""),width=10, height=5, pointsize=8)}
	
	plot(GL.plot)
	if(length(device) != 0 && device != "quartz") dev.off()
		
 	AD.plot <- xyplot(fisher.matrix$amps+fisher.matrix$dels~fisher.matrix$MB|fisher.matrix$phenos, type="h", layout=c(1,length(phenos)+1), col=c("dark green","red"), main=amp.title, xlab="Chromosome", ylab="Frequency",
 	scales=list(x=list(at=mids,labels=chr.nms, cex=0.5)), lwd=3, ylim=c(-1.1,1.1),
 	panel = function(...)
 	{panel.abline(h = 0, lty = 1);panel.abline(v=c(1,ends),col="grey");panel.xyplot(...)})
	
	if(device == "quartz") {quartz(title=paste(project,gain.title),width=6, height=10)}
	if(device == "PNG")
		{AD.file <- paste(project,"AD.Fishers","png", sep=".")
		trellis.device("png", color=T, file=paste("Fishers.plots/",AD.file,sep=""),width=1600, height=800, pointsize=12)}
	if(device == "JPEG")
		{AD.file <- paste(project,"AD.Fishers","jpeg", sep=".")
		trellis.device("jpeg", color=T, file=paste("Fishers.plots/",AD.file,sep=""),width=1600, height=800, pointsize=12)}
	if(device == "PDF")
		{AD.file <- paste(project,"AD.Fishers","pdf", sep=".")
		trellis.device("pdf", color=T, file=paste("Fishers.plots/",AD.file,sep=""),width=8, height=6, pointsize=8)}
	if(device == "PS")
		{AD.file <- paste(project,"AD.Fishers","ps", sep=".")
		trellis.device("postscript", color=T, file=paste("Fishers.plots/",AD.file,sep=""),width=10, height=5, pointsize=8)}
	
	plot(AD.plot)
 	if(length(device) != 0 && device != "quartz") dev.off()
	
}


pairedCGH <- function(cgh, control.sample, test.sample, file,
	acgh.cnv.table="acgh.cnvs.txt", full.cnv.table="full.cnvs.txt", genes.table="all.genes.txt", cytoband.table="cytobands.txt", mirnas.table="mirnas.txt")
{
	require(Biobase)
	
	if(!is.numeric(control.sample)) {
		control.sample <- match(control.sample, sampleNames(cgh))
		if(is.na(control.sample)) stop("Control sample is not found in the ExpressionSet")
		}
	if(!is.element(control.sample,1:ncol(cgh))) stop("Control sample is not found in the ExpressionSet")
	
	if(!is.numeric(test.sample)) {
		test.sample <- match(test.sample, sampleNames(cgh))
		if(is.na(test.sample)) stop("Test sample is not found in the ExpressionSet")
		}
	if(!is.element(test.sample,1:ncol(cgh))) stop("Test sample is not found in the ExpressionSet")
	
	diffs <- which(assayData(cgh)$GL[,control.sample] != assayData(cgh)$GL[,test.sample])
	c.states <- paste(assayData(cgh)$GL[,control.sample],assayData(cgh)$GL[,test.sample])
	rle.result <- rle(c.states)
	rbind(rle.result$lengths, rle.result$values)
	states.to <- cumsum(rle.result$lengths)
	states.from <- cumsum(rle.result$lengths)-rle.result$lengths-1
	states.from[1] <- 1
	result <- data.frame(states=rle.result$values, from=states.from, to=states.to)
	
	result$control <- assayData(cgh)$GL[states.to,control.sample]
	result$test <- assayData(cgh)$GL[states.to,test.sample]
	
	result <- result[which(result$control != result$test),] 
	result$diff <- abs(result$control-result$test)
	
	genes.table <- read.delim(genes.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	mirnas.table <- read.delim(mirnas.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	cytobands.table <- read.delim(cytoband.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	acgh.cnv.table <- read.delim(acgh.cnv.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	full.cnv.table <- read.delim(full.cnv.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	
	cat(file=file,c("Changes in aCGH states - ", sampleNames(cgh)[control.sample], " vs ", sampleNames(cgh)[test.sample], "\n"), append=F, sep="")
	cat(file=file,c("chrom", "start", "end","start.MB", "end.MB","start.bac", "end.bac","BACs","length MB", sampleNames(cgh)[control.sample], sampleNames(cgh)[test.sample],"diff.aCGH", "genes", 	"mirnas","acgh.cnvs","full.cnvs", "cytobands"),"\n", append=T, sep="\t")
		
	
	for(i in 1:nrow(result))
	{

		start.index <- result$from[i]	
		end.index <- result$to[i]
		iBACs <- end.index-start.index
		start <- fData(cgh)$start[start.index]	
		end <- fData(cgh)$end[end.index]		
		ichrom <- fData(cgh)$chrom[start.index]	
		start.MB <- round(start/1000000,2)
		end.MB <- round(end/1000000,2)
		ilength <- (end-start)/1000000
		start.bac <- fData(cgh)$bac.id[start.index]
		end.bac <- fData(cgh)$bac.id[end.index]
		
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end > start),]
		igenes <- igenes.max[which(igenes.max$start < end),]
		igenes <- igenes$symbol
		
		imirnas.chrom <- mirnas.table[which(mirnas.table$chrom == ichrom),]
		imirnas.max <- imirnas.chrom[which(imirnas.chrom$end > start),]
		imirnas <- imirnas.max[which(imirnas.max$start < end),]
		imirnas <- imirnas$mirna
		
		i.full.cnv.chrom <- full.cnv.table[which(full.cnv.table$chrom == ichrom),]
		i.full.cnv.max <- i.full.cnv.chrom[which(i.full.cnv.chrom$end > start),]
		i.full.cnvs <- i.full.cnv.max[which(i.full.cnv.max$start < end),]
		i.full.cnvs <- i.full.cnvs$cnv
		
		i.acgh.cnv.chrom <- acgh.cnv.table[which(acgh.cnv.table$chrom == ichrom),]
		i.acgh.cnv.max <- i.acgh.cnv.chrom[which(i.acgh.cnv.chrom$end > start),]
		i.acgh.cnvs <- i.acgh.cnv.max[which(i.acgh.cnv.max$start < end),]
		i.acgh.cnvs <- i.acgh.cnvs$cnv
		
		ibands.chrom <- cytobands.table[which(cytobands.table$chrom == ichrom),]
		ibands.max <- ibands.chrom[which(ibands.chrom$end > start),]
		ibands <- ibands.max[which(ibands.max$start < end),]
		ibands <- ibands$cytoband
		if(length(ibands)>1)
			{ibands <- paste(ibands[1],ibands[length(ibands)],sep="-")}
		
		row.result <- c(ichrom, start, end, start.MB, end.MB, start.bac, end.bac,iBACs, ilength,
			result$control[i], result$test[i], result$diff[i])
		
		cat(file=file, row.result, append=T, sep="\t")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, igenes,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, imirnas,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.acgh.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, i.full.cnvs,append=T, sep=", ")
		cat(file=file, "\t",append=T, sep="")
		cat(file=file, ibands,append=T, sep="-")
		cat(file=file, "\n",append=T, sep="\t")

	}

}

rescaleCGHtoMAD <- function(cgh, center="median", rescale=NULL)
#	Rescale CGH to make the MAD the same in each case
{
#	commands <- c(cgh$commands, deparse(match.call()))
	require(Biobase)
	
	if(is.element("smo", assayDataElementNames(cgh))) new.smo <- assayData(cgh)$smo

	
	if(is.null(pData(cgh)$MAD)) cgh <- calculateMAD(cgh)
	if(center=="mean") project.MAD <- mean(pData(cgh)$MAD)
	if(center=="median") project.MAD <- median(pData(cgh)$MAD)
	if(is.null(rescale))
		{
		for(i in 1:ncol(exprs(cgh)))
			{	
			exprs(cgh)[,i] <- exprs(cgh)[,i]*(project.MAD/pData(cgh)$MAD[i])
			if(is.element("smo", assayDataElementNames(cgh))) new.smo[,i] <- assayData(cgh)$smo[,i]*(project.MAD/pData(cgh)$MAD[i])
			}
		}else{
		for(i in 1:ncol(exprs(cgh)))
			{	
			exprs(cgh)[,i] <- exprs(cgh)[,i]*(rescale/pData(cgh)$MAD[i])
			if(is.element("smo", assayDataElementNames(cgh))) new.smo[,i] <- assayData(cgh)$smo[,i]*(rescale/pData(cgh)$MAD[i])
			}
		}
	
	if(is.element("smo", assayDataElementNames(cgh))) assayDataElement(cgh, "smo") <- new.smo
	
	if(is.element("GL", assayDataElementNames(cgh))) assayDataElement(cgh, "GL") <- NULL
	cgh <- calculateMAD(cgh)
	return(cgh)
}


medianCGH <- function(expression, cgh)
#	Calculate median CGH values for each probe in an expression object
{
	require(Biobase)
	options(warn=-1)
#	commands <- c(expression@commands,deparse(match.call()))
	cat("Matching expression and aCGH sampleNames\n")
	cat("Expression sampleNames\n")
	cat(sampleNames(expression),"\n")
	cat("aCGH sampleNames\n")
	cat(sampleNames(cgh),"\n")
	
	exp.match <- match(sampleNames(expression), sampleNames(cgh), nomatch=0)
	cgh <- cgh[,exp.match]
	
	cgh.match <- match(sampleNames(cgh), sampleNames(expression), nomatch=0)
	expression <- expression[,cgh.match]
	
	expression <- expression[order(fData(expression)$chrom, fData(expression)$start),]	
	cat(ncol(expression),"matching sampleNames\n")
	cat(sampleNames(expression),"\n")

	if(is.element("smo", assayDataElementNames(cgh))) cgh.smo <- matrix(0,nrow(expression),ncol(expression))

	cgh.exprs <- matrix(0,nrow(expression),ncol(expression))
	
	cat("Calculating median CGH values for expression probes\n")
	
	for (i in 1:nrow(expression))
	{
		iprobes.cgh.chrom <- fData(cgh)[which(fData(cgh)$chrom == fData(expression)$chrom[i]),,drop=F]
		
		iprobes.cgh.plus <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$end > fData(expression)$start[i]),,drop=F]
		iprobes.cgh.minus <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$start < fData(expression)$end[i]),,drop=F]
		
		iprobes.cgh.overlap <- fData(cgh)[which(is.element(featureNames(cgh), intersect(row.names(iprobes.cgh.plus),row.names(iprobes.cgh.minus)))),,drop=F]
		
		iprobes.cgh.more <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$start > fData(expression)$start[i]),,drop=F]
		iprobes.cgh.less <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$end < fData(expression)$end[i]),,drop=F]

		iprobes.cgh.more <- iprobes.cgh.more[order(iprobes.cgh.more$start),,drop=F]
		iprobes.cgh.less <- iprobes.cgh.less[rev(order(iprobes.cgh.less$start)),,drop=F]
		
		if(nrow(iprobes.cgh.overlap) <1)
			{
			iprobes.cgh <- rbind(iprobes.cgh.less[1,],iprobes.cgh.more[1,])
			iprobes.cgh <- iprobes.cgh[order(iprobes.cgh$start),]
			}else{iprobes.cgh <- iprobes.cgh.overlap[order(iprobes.cgh.overlap$start),]}
		
		iM <- as.data.frame(exprs(cgh)[which(is.element(featureNames(cgh),row.names(iprobes.cgh))),,drop=F])
		med.M <- apply(iM,2,median, na.rm=T)
		cgh.exprs[i,] <- as.vector(med.M)
		
		if(is.element("smo", assayDataElementNames(cgh)))
		{
		iSMO <- as.data.frame(assayData(cgh)$smo[which(is.element(featureNames(cgh),row.names(iprobes.cgh))),,drop=F])
		med.smo <- apply(iSMO,2,median, na.rm=T)
		cgh.smo[i,] <- as.vector(med.smo)
		}
		if(i %% 100 == 0) cat(i," ")
		cat(i)
		
	}
	
	expression.cgh <- new("BACE.exp.cgh", exprs=exprs(expression), featureData=featureData(expression), phenoData=phenoData(expression))
	
	assayData(expression.cgh) <- assayData(expression)
	row.names(cgh.exprs) <- featureNames(expression.cgh)
	colnames(cgh.exprs) <- sampleNames(expression.cgh)
	assayDataElement(expression.cgh, "cgh.exprs") <- cgh.exprs
	if(is.element("thresholds", slotNames(cgh))) expression.cgh@thresholds <- cgh@thresholds
	row.names(cgh.smo) <- featureNames(expression.cgh)
	colnames(cgh.smo) <- sampleNames(expression.cgh)
	assayDataElement(expression.cgh, "smo") <- cgh.smo
	
	if(is.element("smo", assayDataElementNames(cgh)))
	{

	iGL <- thresholdCGH(cgh=expression.cgh, 
						gainthresh=cgh@thresholds$gainthresh,
						lossthresh=cgh@thresholds$lossthresh,
						ampthresh=cgh@thresholds$ampthresh,
						delthresh=cgh@thresholds$delthresh,
						verbose=F)$GL
	row.names(iGL) <- featureNames(expression.cgh)
	colnames(iGL) <- sampleNames(expression.cgh)
	assayDataElement(expression.cgh, "GL") <- iGL
	}
	
#	expression.cgh$commands <- commands
	expression.cgh
}


listAmpGenes <- function(cgh, amp.count, genes.table="all.genes.txt", gene.id.match="symbol", pheno=NULL, project)
{
	require(Biobase)
	cat("Counting genes in regions amplified in",amp.count,"or more cases\n")
	glad <- matrix(0,nrow(cgh), 4)
	
	ampGL <- apply(assayData(cgh)$GL,1,function(x) sum(x==2))
	amps <- rep(NA, nrow(cgh))
	amps[which(ampGL >= amp.count)] <- 1

	amp.breaks <- listBreaks(cgh=cgh, GL.column=amps, contig=1)
	amp.breaks <- amp.breaks[,which(amp.breaks[2,] > 0),drop=F]

	genes.table <- read.delim(genes.table, sep="\t", header=T, row.names=NULL, stringsAsFactors=F)
	
	if(ncol(amp.breaks) > 0)
	{
		end.index <- amp.breaks[1,1]
		start.index <- amp.breaks[1,1]-(amp.breaks[3,1]-1)
		ichrom <- fData(cgh)$chrom[start.index]
		start <- fData(cgh)$start[start.index]	
		end <- fData(cgh)$end[end.index]
		
		genes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		genes.max <- genes.chrom[which(genes.chrom$end > start),,drop=F]
		genes <- genes.max[which(genes.max$start < end),,drop=F]
		
		if(ncol(amp.breaks) > 1) for (i in 1:ncol(amp.breaks))
		{		
		end.index <- amp.breaks[1,i]
		start.index <- amp.breaks[1,i]-(amp.breaks[3,i]-1)
		ichrom <- fData(cgh)$chrom[start.index]
		start <- fData(cgh)$start[start.index]	
		end <- fData(cgh)$end[end.index]
		
		igenes.chrom <- genes.table[which(genes.table$chrom == ichrom),]
		igenes.max <- igenes.chrom[which(igenes.chrom$end > start),,drop=F]
		igenes <- igenes.max[which(igenes.max$start < end),,drop=F]
		if(is.null(genes)) {genes <- igenes}else{
			genes <- rbind(genes,igenes)}
		# igenes is the subset of the genes.table for the amp.breaks
		}

		
		uni.genes <- genes[!duplicated(genes),]
		cat("Listing",nrow(uni.genes),"genes in amplified regions\n")
		
		maxM <- rep(NA,nrow(uni.genes))
		amp.cases <- rep(NA,nrow(uni.genes))
		amps <- rep(NA,nrow(uni.genes))
		if(!is.null(pheno))
			{
			phenos <- unique(pheno)
			amp.count <- matrix(0,nrow(uni.genes),length(phenos))
			}

	for (i in 1:nrow(uni.genes))

	{
		iprobes.cgh.chrom <- fData(cgh)[which(fData(cgh)$chrom == uni.genes$chrom[i]),,drop=F]
		
		iprobes.cgh.plus <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$end >= uni.genes$start[i]),,drop=F]
		iprobes.gene.Minus <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$start <= uni.genes$end[i]),,drop=F]
		
		iprobes.cgh.overlap <- fData(cgh)[which(is.element(featureNames(cgh), intersect(row.names(iprobes.cgh.plus),row.names(iprobes.gene.Minus)))),,drop=F]
		
		iprobes.gene.More <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$start > uni.genes$start[i]),,drop=F]
		iprobes.cgh.less <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$end < uni.genes$end[i]),,drop=F]

		iprobes.gene.More <- iprobes.gene.More[order(iprobes.gene.More$start),,drop=F]
		iprobes.cgh.less <- iprobes.cgh.less[rev(order(iprobes.cgh.less$start)),,drop=F]
		
		if(nrow(iprobes.cgh.overlap) <1)
			{
			iprobes.cgh <- rbind(iprobes.cgh.less[1,],iprobes.gene.More[1,])
			iprobes.cgh <- iprobes.cgh[order(iprobes.cgh$start),]
			}else{iprobes.cgh <- iprobes.cgh.overlap[order(iprobes.cgh.overlap$start),]}
			
		iM <- as.data.frame(exprs(cgh)[which(is.element(featureNames(cgh),row.names(iprobes.cgh))),,drop=F])
		maxM[i] <- max(iM, na.rm=T)
		
		
		iGL <- as.data.frame(assayData(cgh)$GL[which(is.element(featureNames(cgh),row.names(iprobes.cgh))),,drop=F])
		topGL <- apply(iGL, 2,max, na.rm=T)
		bottomGL <- apply(iGL, 2,min, na.rm=T)
		
		a.cases <- sampleNames(cgh)[topGL == 2]
		amps[i] <- length(a.cases)
		amp.cases[i] <- paste(a.cases, collapse=", ")
		if(!is.null(pheno))
		{
		phenos <- unique(as.character(pheno))
		a.count <- rep(0, length(phenos))
		for (p in 1:length(phenos))
		{
		a.count[p] <- sum(topGL[phenos[p]==pheno]==2)
		}
		amp.count[i,] <- a.count
		}		
		if(i %% 100 == 0) cat(i," ")
	}
		uni.genes$maxM <- maxM
		uni.genes$amp.cases <- amp.cases
		uni.genes$frequency.all <- amps		
		if(!is.null(pheno))
			{
			amp.count.d <- data.frame(amp.count)
			names(amp.count.d) <- paste("frequency",phenos,sep=".")
			uni.genes <- cbind(uni.genes,amp.count.d)
			}
		filename <- paste(project, "amp.genes", "xls", sep=".")
		write.table(uni.genes, file=filename, sep="\t", row.names=F)

	}

}

geneCGH <- function(cgh, genes.table="all.genes.txt", project, return.eset=T)
#	Function to list aCGH values and copy number for an external list of input genes
{

	require(Biobase)
	genes.table <- read.delim(genes.table, sep="\t", header=T, stringsAsFactors=F)
	row.names(genes.table) <- make.names(genes.table[,1], unique=T)
	uni.genes <- genes.table[!duplicated(row.names(genes.table)),]
	cat("Listing aCGH states for",nrow(uni.genes),"genes\n")
	
	if(is.element("smo", assayDataElementNames(cgh))) smo <- assayData(cgh)$smo
	
	cgh.exprs <- cgh.smooth <- matrix(NA,nrow(uni.genes),ncol(cgh))
	
	for (i in 1:nrow(uni.genes)){
		
		iprobes.cgh.chrom <- fData(cgh)[which(fData(cgh)$chrom == uni.genes$chrom[i]),,drop=F]
		
		iprobes.cgh.plus <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$end >= uni.genes$start[i]),,drop=F]
		iprobes.gene.Minus <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$start <= uni.genes$end[i]),,drop=F]
		
		iprobes.cgh.overlap <- fData(cgh)[which(is.element(featureNames(cgh), intersect(row.names(iprobes.cgh.plus),row.names(iprobes.gene.Minus)))),,drop=F]
		
		iprobes.gene.More <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$start > uni.genes$start[i]),,drop=F]
		iprobes.cgh.less <- iprobes.cgh.chrom[which(iprobes.cgh.chrom$end < uni.genes$end[i]),,drop=F]

		iprobes.gene.More <- iprobes.gene.More[order(iprobes.gene.More$start),,drop=F]
		iprobes.cgh.less <- iprobes.cgh.less[rev(order(iprobes.cgh.less$start)),,drop=F]
		
		if(nrow(iprobes.cgh.overlap) <1)
			{
			iprobes.cgh <- rbind(iprobes.cgh.less[1,],iprobes.gene.More[1,])
			iprobes.cgh <- iprobes.cgh[order(iprobes.cgh$start),]
			}else{iprobes.cgh <- iprobes.cgh.overlap[order(iprobes.cgh.overlap$start),]}
			
		iM <- as.data.frame(exprs(cgh)[which(is.element(featureNames(cgh),row.names(iprobes.cgh))),,drop=F])
		med.M <- apply(iM,2,median, na.rm=T)
		cgh.exprs[i,] <- as.vector(med.M)

		iSMO <- smo[which(is.element(featureNames(cgh),row.names(iprobes.cgh))),,drop=F]
		cgh.smooth[i,] <- as.vector(apply(iSMO,2,median, na.rm=T))

		if(i %% 100 == 0) cat(i," ")
		
	}
	
	cgh.fdata <- new("AnnotatedDataFrame", data=uni.genes)
	varMetadata(cgh.fdata)$labelDescription[grep("chrom",names(uni.genes))] <- "Chromosome"
	varMetadata(cgh.fdata)$labelDescription[grep("start",names(uni.genes))] <- "Start position (bp)"
	varMetadata(cgh.fdata)$labelDescription[grep("end",names(uni.genes))] <- "End position (bp)"
	varMetadata(cgh.fdata)$labelDescription[grep("symbol",names(uni.genes))] <- "Gene symbol"
	varMetadata(cgh.fdata)$labelDescription[grep("ensg",names(uni.genes))] <- "Ensembl gene ID"
	varMetadata(cgh.fdata)$labelDescription[grep("cytoband",names(uni.genes))] <- "Cytoband"
	
	genes.cgh <- new("BACE.cgh", exprs=cgh.exprs, phenoData=phenoData(cgh), featureData=cgh.fdata, annotation=cgh@annotation)
	
	if(is.element("smo", assayDataElementNames(cgh))) assayDataElement(genes.cgh, "smo") <- cgh.smooth
	
	GL.result <- thresholdCGH(genes.cgh, gainthresh=cgh@thresholds$gainthresh,
		lossthresh=cgh@thresholds$lossthresh,
		ampthresh=cgh@thresholds$ampthresh,
		delthresh=cgh@thresholds$delthresh, verbose=F)

	assayDataElement(genes.cgh, "GL") <- GL.result$GL
	
	genes.GLAD <- countGLAD(genes.cgh)
	
	genes.smo <- assayData(genes.cgh)$smo
	colnames(genes.smo) <- paste(sampleNames(genes.cgh), "smo", sep=".")
	genes.GL <- assayData(genes.cgh)$GL
	colnames(genes.GL) <- paste(sampleNames(genes.cgh), "GL", sep=".")
	
	genes.out <- cbind(genes.table, genes.GLAD, genes.smo, genes.GL)
	write.table(genes.out, file=paste(project, "genes.CGH.xls", sep="."), row.names=F, sep="\t", na="")
	
	if(return.eset==T){
		return(genes.cgh)
		}

}


expressionCGHCorrelation <- function(exp.cgh, correlation.method="pearson", adj.method="BH", pval.thresh=0.05)
#	exp.cgh is a list object with matched data.frames for expression and CGH
#	exp.cgh is a result of medianCGH
#	
{
	require(Biobase)
	exp.cgh <- exp.cgh[order(fData(exp.cgh)$chrom, fData(exp.cgh)$start),]
	
	options(warn=-1)
	require(multtest)
#	commands <- c(exp.cgh$commands,deparse(match.call()))
	
	cat("Correlation of gene expression with aCGH\n")
	
	pearson.p <- rep(NA,nrow(exp.cgh))
	pearson.cor <- rep(NA,nrow(exp.cgh))
	spearman.p <- rep(NA,nrow(exp.cgh))
	spearman.cor <- rep(NA,nrow(exp.cgh))

	for (i in 1:nrow(exp.cgh))
	
		{
			if(sum(is.finite(exprs(exp.cgh)[i,])) > 2)
			{
				ok <- is.finite(exprs(exp.cgh)[i,])
				if(correlation.method=="pearson"){
				pearson.score <- cor.test(exprs(exp.cgh)[i,ok],assayData(exp.cgh)$smo[i,ok], method="pearson", na.rm=T)
				pearson.p[i] <- pearson.score$p.value
				pearson.cor[i] <- pearson.score$estimate}
				if(correlation.method=="spearman"){
				spearman.score <- cor.test(exprs(exp.cgh)[i,ok],assayData(exp.cgh)$smo[i,ok], method="spearman", na.rm=T)
				spearman.p[i]<- spearman.score$p.value
				spearman.cor[i] <- spearman.score$estimate}
				if(i %% 100 == 0) cat(i," ")
			}
		
		}
		
	cat("\nDone\n")
	
	if(correlation.method=="pearson"){
	fData(exp.cgh)$pearson.p <- pearson.p
	fData(exp.cgh)$pearson.adjp <- p.adjust(pearson.p, method=adj.method)
	fData(exp.cgh)$pearson.cor <- pearson.cor
	cat(length(which(fData(exp.cgh)$pearson.adjp < pval.thresh)),"significantly correlated genes\n")
	}
	if(correlation.method=="spearman"){
	fData(exp.cgh)$spearman.p <- spearman.p
	fData(exp.cgh)$spearman.adjp <- p.adjust(spearman.p, method=adj.method)
	fData(exp.cgh)$spearman.cor <- spearman.cor
	cat(length(which(fData(exp.cgh)$spearman.adjp < pval.thresh)),"significantly correlated genes\n")
	}

#	exp.cgh$commands <- commands
	cat("\nDone\n")

	exp.cgh
}



wilcoxTestCGH <- function(exp.cgh, threshold.data=T, count=1, pval.thresh=0.05)
#	Function to run Wilcoxon rank sum tests 
#	to correlate gene expression with aCGH states
{
	require(Biobase)
	require(multtest)
#	commands <- c(exp.cgh$commands,deparse(match.call()))
	cat("Wilcoxon rank sum tests of expression based on copy number changes\n")
	exp.cgh <- exp.cgh[order(fData(exp.cgh)$chrom, fData(exp.cgh)$start),]
	
	wilcox.p <- matrix(NA,nrow(exp.cgh),4)
	gain.fold <- rep(NA,nrow(exp.cgh))
	loss.fold <- rep(NA,nrow(exp.cgh))
	amp.fold <- rep(NA,nrow(exp.cgh))
	del.fold <- rep(NA,nrow(exp.cgh))
	
	if(threshold.data == T)
	{
			
	for (i in 1:nrow(exp.cgh))
		{
		
		gains <- which(assayData(exp.cgh)$GL[i,] >= 1)
		if((sum(!is.na(exprs(exp.cgh)[i,gains])) >= count) && (sum(!is.na(exprs(exp.cgh)[i,-gains])) >= count)){
			wilcox.p[i,1] <- wilcox.test(exprs(exp.cgh)[i,gains],exprs(exp.cgh)[i,-gains])$p.value}
		
		losses <- which(assayData(exp.cgh)$GL[i,] <= -1)
		if((sum(!is.na(exprs(exp.cgh)[i,losses])) >= count) && (sum(!is.na(exprs(exp.cgh)[i,-losses])) >= count)){
			wilcox.p[i,2] <- wilcox.test(exprs(exp.cgh)[i,losses],exprs(exp.cgh)[i,-losses])$p.value}
		
		amps <- which(assayData(exp.cgh)$GL[i,]== 2)
		if((sum(!is.na(exprs(exp.cgh)[i,amps])) >= count) && (sum(!is.na(exprs(exp.cgh)[i,-amps])) >= count)){
			wilcox.p[i,3] <- wilcox.test(exprs(exp.cgh)[i,amps],exprs(exp.cgh)[i,-amps])$p.value}
		
		dels <- which(assayData(exp.cgh)$GL[i,]== -2)
		if((sum(!is.na(exprs(exp.cgh)[i,dels])) > count) && (sum(!is.na(exprs(exp.cgh)[i,-dels])) > count)){
			wilcox.p[i,4] <- wilcox.test(exprs(exp.cgh)[i,dels],exprs(exp.cgh)[i,-dels])$p.value}
		
		gain.fold[i] <- 2^(mean(exprs(exp.cgh)[i,gains],na.rm=T))/2^(mean(exprs(exp.cgh)[i,-gains],na.rm=T))
		loss.fold[i] <- 2^(mean(exprs(exp.cgh)[i,losses],na.rm=T))/2^(mean(exprs(exp.cgh)[i,-losses],na.rm=T))
		amp.fold[i] <- 2^(mean(exprs(exp.cgh)[i,amps],na.rm=T))/2^(mean(exprs(exp.cgh)[i,-amps],na.rm=T))
		del.fold[i] <- 2^(mean(exprs(exp.cgh)[i,dels],na.rm=T))/2^(mean(exprs(exp.cgh)[i,-dels],na.rm=T))
		
		if(i %% 100 == 0) cat(i," ")
		
		}
	}
	
	if(threshold.data == F)
	{
	for (i in 1:nrow(exprs(exp.cgh)))
		{

		if((sum(!is.na(exprs(exp.cgh)[i,assayData(exp.cgh)$GL[i,] == 1])) >= count) && (sum(!is.na(exprs(exp.cgh)[i,assayData(exp.cgh)$GL[i,] == 0])) >= count))
		{wilcox.p[i,1] <- wilcox.test(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 1)],exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 0)])$p.value
			}
		if((sum(!is.na(exprs(exp.cgh)[i,assayData(exp.cgh)$GL[i,] == -1])) >= count) && (sum(!is.na(exprs(exp.cgh)[i,assayData(exp.cgh)$GL[i,] == 0])) >= count))
		{wilcox.p[i,2] <- wilcox.test(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== -1)],exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 0)])$p.value
			}
		if((sum(!is.na(exprs(exp.cgh)[i,assayData(exp.cgh)$GL[i,] == 2])) >= count) && (sum(!is.na(exprs(exp.cgh)[i,assayData(exp.cgh)$GL[i,] == 0])) >= count))
		{wilcox.p[i,3] <- wilcox.test(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 2)],exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 0)])$p.value
			}
		if((sum(!is.na(exprs(exp.cgh)[i,assayData(exp.cgh)$GL[i,] == -2])) >= count) && (sum(!is.na(exprs(exp.cgh)[i,assayData(exp.cgh)$GL[i,] == 0])) >= count))
		{wilcox.p[i,4] <- wilcox.test(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== -2)],exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 0)])$p.value
			}
		
		gain.fold[i] <- 2^(mean(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 1)],na.rm=T))/2^(mean(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 0)],na.rm=T))
		loss.fold[i] <- 2^(mean(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== -1)],na.rm=T))/2^(mean(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 0)],na.rm=T))
		amp.fold[i] <- 2^(mean(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 2)],na.rm=T))/2^(mean(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 0)],na.rm=T))
		del.fold[i] <- 2^(mean(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== -2)],na.rm=T))/2^(mean(exprs(exp.cgh)[i,which(assayData(exp.cgh)$GL[i,]== 0)],na.rm=T))
		
		if(i %% 100 == 0) cat(i," ")
		}
	}
	
	cat("\nDone\n")
	
	wilcox.adjp <- wilcox.p
	
	gain.counts <- as.vector(apply(assayData(exp.cgh)$GL >= 1, 1, sum))
	loss.counts <- as.vector(apply(assayData(exp.cgh)$GL <= -1, 1, sum))
	amp.counts <- as.vector(apply(assayData(exp.cgh)$GL == 2, 1, sum))
	del.counts <- as.vector(apply(assayData(exp.cgh)$GL == -2, 1, sum))
	
	gain.counts[gain.counts < count] <- 0
	gain.counts[gain.counts >= count] <- 1
	loss.counts[loss.counts < count] <- 0
	loss.counts[loss.counts >= count] <- 1
	amp.counts[amp.counts < count] <- 0
	amp.counts[amp.counts >= count] <- 1
	del.counts[del.counts < count] <- 0
	del.counts[del.counts >= count] <- 1
	
	gain.breaks <- listBreaks(exp.cgh, GL.column=gain.counts, contig=1)
	gain.breaks <- gain.breaks[,which(gain.breaks[2,] > 0), drop=F]
	loss.breaks <- listBreaks(cgh=exp.cgh, GL.column=loss.counts, contig=1)
	loss.breaks <- loss.breaks[,which(loss.breaks[2,] > 0) ,drop=F]
	amp.breaks <- listBreaks(cgh=exp.cgh, GL.column=amp.counts, contig=1)
	amp.breaks <- amp.breaks[,which(amp.breaks[2,] > 0),drop=F]
	del.breaks <- listBreaks(cgh=exp.cgh, GL.column=del.counts, contig=1)
	del.breaks <- del.breaks[,which(del.breaks[2,] > 0),drop=F]
	
	if(ncol(gain.breaks) > 0) for (i in 1:ncol(gain.breaks))
		{		
		end.index <- gain.breaks[1,i]
		start.index <- gain.breaks[1,i]-(gain.breaks[3,i]-1)
		pvals.i <- wilcox.p[start.index:end.index,1]
		if(sum(!is.na(pvals.i)) >1)
		{wilcox.adjp[start.index:end.index,1] <- p.adjust(pvals.i, method="BH")}
		}
	if(ncol(loss.breaks) > 0) for (i in 1:ncol(loss.breaks))
		{		
		end.index <- loss.breaks[1,i]
		start.index <- loss.breaks[1,i]-(loss.breaks[3,i]-1)
		pvals.i <- wilcox.p[start.index:end.index,2]
		if(sum(!is.na(pvals.i)) >1)
		{wilcox.adjp[start.index:end.index,2] <- p.adjust(pvals.i, method="BH")}
		}
	if(ncol(amp.breaks) > 0) for (i in 1:ncol(amp.breaks))
		{		
		end.index <- amp.breaks[1,i]
		start.index <- amp.breaks[1,i]-(amp.breaks[3,i]-1)
		pvals.i <- wilcox.p[start.index:end.index,3]
		if(sum(!is.na(pvals.i)) >1)
		{wilcox.adjp[start.index:end.index,3] <- p.adjust(pvals.i, method="BH")}
		}
	if(ncol(del.breaks) > 0) for (i in 1:ncol(del.breaks))
		{		
		end.index <- del.breaks[1,i]
		start.index <- del.breaks[1,i]-(del.breaks[3,i]-1)
		pvals.i <- wilcox.p[start.index:end.index,4]
		if(sum(!is.na(pvals.i)) >1)
		{wilcox.adjp[start.index:end.index,4] <- p.adjust(pvals.i, method="BH")}
		}
	
	fData(exp.cgh)$gain.fold <- gain.fold
	fData(exp.cgh)$loss.fold <- loss.fold
	fData(exp.cgh)$amp.fold <- amp.fold
	fData(exp.cgh)$del.fold <- del.fold
	
	if(threshold.data == T)
	{
		fData(exp.cgh)$Wilcox.p.gain <- wilcox.p[,1]
		fData(exp.cgh)$Wilcox.adjp.gain <- wilcox.adjp[,1]
		fData(exp.cgh)$Wilcox.p.loss <- wilcox.p[,2]
		fData(exp.cgh)$Wilcox.adjp.loss <- wilcox.adjp[,2]
		fData(exp.cgh)$Wilcox.p.amp <- wilcox.p[,3]
		fData(exp.cgh)$Wilcox.adjp.amp <- wilcox.adjp[,3]
		fData(exp.cgh)$Wilcox.p.del <- wilcox.p[,4]
		fData(exp.cgh)$Wilcox.adjp.del <- wilcox.adjp[,4]
		cat(length(which(fData(exp.cgh)$Wilcox.p.amp < pval.thresh)),"genes overexpressed when amplified\n")
	}
	
	if(threshold.data == F)
	{
		fData(exp.cgh)$Wilcox.p.gain.vs.NC <- wilcox.p[,1]
		fData(exp.cgh)$Wilcox.adjp.gain.vs.NC <- wilcox.adjp[,1]
		fData(exp.cgh)$Wilcox.p.loss.vs.NC <- wilcox.p[,2]
		fData(exp.cgh)$Wilcox.adjp.loss.vs.NC <- wilcox.adjp[,2]
		fData(exp.cgh)$Wilcox.p.amp.vs.NC <- wilcox.p[,3]
		fData(exp.cgh)$Wilcox.adjp.amp.vs.NC <- wilcox.adjp[,3]
		fData(exp.cgh)$Wilcox.p.del.vs.NC <- wilcox.p[,4]
		fData(exp.cgh)$Wilcox.adjp.del.vs.NC <- wilcox.adjp[,4]
		cat(length(which(fData(exp.cgh)$Wilcox.p.amp.vs.NC < pval.thresh)),"genes overexpressed when amplified\n")
	}

#	exp.cgh$commands <- commands
	cat("\nDone\n")
	exp.cgh
}

expCGHCorrelationPlot <- function(exp.cgh, lookup.string, ann.column="symbol", cor.method="pearson")
{
	lookup <- paste("\\^", lookup.string, "\\$", sep="")
	if(cor.method == "pearson") exp.cgh <- exp.cgh[order(fData(exp.cgh)$pearson.p),]
	if(cor.method == "spearman") exp.cgh <- exp.cgh[order(fData(exp.cgh)$spearman.p),]
	index <- match(lookup.string, fData(exp.cgh)[,match(ann.column,names(fData(exp.cgh)))])
	if(length(index) < 1 ) stop("No probes matching input found")
	pdf(file=paste(lookup.string,"exp.cgh.plot.pdf", sep="."), width=6, height=5.5, pointsize=12)

	par(mar=c(4,4,4,6),xpd=T)

	plot(assayData(exp.cgh)$smo[index,],exprs(exp.cgh)[index,], ylab="Gene Expression", xlab="aCGH", pch=16, main=fData(exp.cgh)[index,match(ann.column,names(fData(exp.cgh)))], las=1, xlim=range(assayData(exp.cgh)$smo[index,])*1.2)
	text(assayData(exp.cgh)$smo[index,],exprs(exp.cgh)[index,], labels=sampleNames(exp.cgh), cex=0.5, adj=1.2)
	
	points(assayData(exp.cgh)$smo[index,which(assayData(exp.cgh)$GL[index,] == -2)],exprs(exp.cgh)[index,which(assayData(exp.cgh)$GL[index,] == -2)], pch=16, col="green")
	points(assayData(exp.cgh)$smo[index,which(assayData(exp.cgh)$GL[index,] == 2)],exprs(exp.cgh)[index,which(assayData(exp.cgh)$GL[index,] == 2)], pch=16, col="red")
	points(assayData(exp.cgh)$smo[index,which(assayData(exp.cgh)$GL[index,] == -1)],exprs(exp.cgh)[index,which(assayData(exp.cgh)$GL[index,] == -1)], pch=16, col="light green")
	points(assayData(exp.cgh)$smo[index,which(assayData(exp.cgh)$GL[index,] == 1)],exprs(exp.cgh)[index,which(assayData(exp.cgh)$GL[index,] == 1)], pch=16, col="pink")
	
	if(cor.method=="pearson"){
		legend("topleft", legend = 
		paste("Pearson cor = ",signif(fData(exp.cgh)$pearson.cor[index],3), "\n",
		"Pearson adj.p = ",signif(fData(exp.cgh)$pearson.adjp[index],3),
		sep=""), bty="n", yjust=1, xjust=1, cex=0.75)}
		
	if(cor.method=="spearman"){
		legend("topleft", legend = 
		paste("Spearman cor = ",signif(fData(exp.cgh)$spearman.cor[index],3), "\n",
		"Spearman adj.p = ",signif(fData(exp.cgh)$spearman.adjp[index],3),
		sep=""), bty="n", yjust=1, xjust=1, cex=0.75)}
		
	legend(max(assayData(exp.cgh)$smo[index,])*1.27,min(exprs(exp.cgh)[index,]),		legend=c("Amp", "Gain", "NC", "Loss", "Del"), fill=c("red", "pink", "black", "light green", "green")
		, bty="n", yjust=0, xjust=0)
		
	dev.off()
}


writeGSEAFiles <- function(eset, pheno, project, fvarLabel="symbol")
{
	require(Biobase)
	eset <- eset[,which(!is.na(pheno))]
	pheno <- pheno[which(!is.na(pheno))]
	gsea.data <- cbind.data.frame(fData(eset)[,fvarLabel], fData(eset)$description, exprs(eset))
	names(gsea.data) <- c("NAME", "Description", sampleNames(eset))
	write.table(gsea.data, file=paste(project,"data.txt", sep="."),row.names=F, sep="\t", quote=F)
	write.table(data.frame(fData(eset)[,fvarLabel]), file=paste(project,"CHIP", sep="."),row.names=F, col.names=F, sep="\t", quote=F)

	cls.file <- paste(project,"cls", sep=".")
	pheno <- gsub(" ", "", pheno)
	cat(ncol(eset), length(unique(pheno)),"1\n",sep=" ", file=cls.file, append=F)
	cat("#",unique(pheno), sep=" ", file=cls.file, append=T)
	cat("\n", file=cls.file, append=T)
	cat(as.character(pheno), sep=" ", file=cls.file, append=T)
}



matrixCorrelationPlot <- function(exp.cgh, cgh.table=assayData(exp.cgh)$smo, cor.method="pearson", project, device)
{
	require(Biobase)
	require(lattice)
	require(latticeExtra)
	require(grid)
	require(marray)
	
	exp.cgh <- exp.cgh[order(fData(exp.cgh)$chrom, fData(exp.cgh)$start),]
	expression.table <- exprs(exp.cgh)

	chroms <- chrom.labs <- unique(fData(exp.cgh)$chrom)
	chrom.labs[chroms == 23] <- "X"
	chrom.labs[chroms == 24] <- "Y"
	chrom.lengths <- as.vector(table(fData(exp.cgh)$chrom))
	chrom.props <- chrom.lengths/nrow(exp.cgh)

	split.levels <- seq(-1,1,by=0.1)
	m.pal <- maPalette(low = "blue", high = "red", mid="white", k =length(split.levels))
	
	if(device == "quartz") {quartz(title=paste(project,"correlation matrix"),width=6, height=6)}
	if(device == "PNG")
		{matrix.file <- paste(project,"matrix","png", sep=".")
		png(file=matrix.file,width=10000, height=10000, pointsize=256)}
	if(device == "JPEG")
		{matrix.file <- paste(project,"matrix","jpeg", sep=".")
		jpeg(file=matrix.file,width=10000, height=10000, pointsize=256, quality=100)}
		
	grid.newpage()
	
	title.vp <- viewport(x = unit(0, "npc"),
		y = unit(0.9,  "npc"),
		width = unit(1, "npc"),
		height = unit(0.1, "npc"),
		just = c("left", "bottom"))
		pushViewport(title.vp)
		grid.text("Correlation matrix", just="top")
		popViewport()

	grid.text(chrom.labs, x=0.1, y=0.15+(0.7*cumsum(chrom.props))-(0.35*chrom.props), gp=gpar(cex=min(1,12/length(chroms))))
	grid.text(chrom.labs, x=0.15+(0.7*cumsum(chrom.props))-(0.35*chrom.props), y=0.1, gp=gpar(cex=min(1,12/length(chroms))))

	grid.text("Gene Expression", x=0.5, y=0.05)
	grid.text("aCGH", x=0.05, y=0.5, rot=90)

	chrom.matrix <-
 	   viewport(width = 0.7, height = 0.7,
		layout = # necessary to fix aspect ratio
		grid.layout(length(chroms), length(chroms),
			widths = chrom.props,
			heights = rev(chrom.props),
			respect = TRUE),
			name = "chrom.matrix")

	pushViewport(chrom.matrix)

	for (i in 1:length(chroms)){

		y.mat.i <- t(as.matrix(cgh.table[which(fData(exp.cgh)$chrom == chroms[i]),]))
		cat(paste("Correlating chromosome",chroms[i],"aCGH values with gene expression values on chromosome\n"))
		for (j in 1:length(chroms)){

			x.mat.j <- t(as.matrix(expression.table[which(fData(exp.cgh)$chrom == chroms[j]),]))
			cor.tab <- cor(x.mat.j,y.mat.i, method = cor.method)
			cor.tab[!is.finite(cor.tab)] <- 0
			cor.grob <- levelplot(t(cor.tab), axes=F, col.regions=m.pal, at=split.levels)
			
			pushViewport(viewport(layout.pos.row = (length(chroms):1)[i], layout.pos.col = j,
				xscale = c(0.5,ncol(cor.tab)+0.5),yscale = c(0.5,nrow(cor.tab)+0.5)))
			grid.rect(gp=gpar(lwd=2))
			do.call("panel.levelplot", trellis.panelArgs(cor.grob, 1))
			
			
			popViewport()
			cat(chroms[j],"\t")
		if(j == length(chroms)) cat("\n")
		}
		
	}

	if(length(device) != 0 && device != "quartz") dev.off()
}


splitHeatmap <- function(exp.cgh, cgh.table="GL", chrom, start, end, project, pvals.to.plot="Wilcox.p.amp", pval.thresh=0.05, cex.label=NULL, heatmap.scale=2, symbol.width=NULL, plot.symbols=T, main=NULL, device=NULL, probeID="symbol")
{

	require(Biobase)
	require(lattice)
	require(grid)
	require(marray)
	
	exp.cgh <- centerGenes(exp.cgh)
	
	chr <- exp.cgh[which(fData(exp.cgh)$chrom == chrom),]
	chr.start <- chr[which(fData(chr)$start >= start),]
	amplicon <- chr.start[which(fData(chr.start)$end <= end),]
	amplicon <- amplicon[order(fData(amplicon)$chrom, fData(amplicon)$start),]

	topGL <- apply(assayData(amplicon)$GL, 2,max)
	sumGL <- apply(assayData(amplicon)$GL, 2,sum)
	
	ampliconM <- exprs(amplicon)[rev(order(fData(amplicon)$chrom, fData(amplicon)$start)),]
	
	amp.matrix <- ampliconM[,which(topGL==2)]
	sum.amp.matrix <- apply(amp.matrix, 2,sum)
	amp.matrix <- amp.matrix[,rev(order(sum.amp.matrix))]
	amp.matrix[is.na(amp.matrix)] <- 0
	
	nonamp.matrix <- ampliconM[,which(topGL!=2)]
	sum.nonamp.matrix <- apply(nonamp.matrix, 2,sum)
	nonamp.matrix <- nonamp.matrix[,rev(order(sum.nonamp.matrix))]
	nonamp.matrix[is.na(nonamp.matrix)] <- 0
	
	split.matrix <-cbind(amp.matrix,rep(NA,nrow(amp.matrix)),nonamp.matrix)
	
	if(cgh.table=="smo") matrix.smo <- as.data.frame(assayData(amplicon)$smo)
	if(cgh.table=="GL") matrix.smo <- as.data.frame(assayData(amplicon)$GL)
	
	matrix.smo <- matrix.smo[rev(order(fData(amplicon)$chrom, fData(amplicon)$start)),]
	
	amp.matrix.smo <- matrix.smo[,which(topGL==2)]
	amp.matrix.smo <- amp.matrix.smo[,rev(order(sum.amp.matrix))]
	
	nonamp.matrix.smo <- matrix.smo[,which(topGL!=2)]
	nonamp.matrix.smo <- nonamp.matrix.smo[,rev(order(sum.nonamp.matrix))]
			
	split.matrix.smo <-cbind(amp.matrix.smo,rep(NA,nrow(amp.matrix.smo)),nonamp.matrix.smo)
	
	pvals <- fData(amplicon)[,which(names(fData(amplicon))==pvals.to.plot)]
	log.pvals <- rev(-(logb(pvals,10)))
	log.thresh <- -(logb(pval.thresh,10))
	log.pvals[is.na(log.pvals)] <- 0
	bars <- log.pvals/max(log.pvals)
	bar.width=1/length(bars)
	line.thresh <- log.thresh/max(log.pvals)
	
	maxM <- (as.integer(max(abs(split.matrix[!is.na(split.matrix)]))))+1
	split.levels <- seq(-heatmap.scale-0.1,heatmap.scale+0.1,by=0.1)
	if(maxM > heatmap.scale) {
		split.matrix[split.matrix > heatmap.scale] <- heatmap.scale
		split.matrix[split.matrix < -heatmap.scale] <- -heatmap.scale
		}
	m.pal <- maPalette(low = "green", high = "red", mid="black", k =length(split.levels))
	gl.pal <- maPalette(low = "green", high = "red", mid="black", k =6)
	
	split.mat <- matrix(split.matrix, nrow(split.matrix), ncol(split.matrix))
	split.mat.smo <- as.matrix(split.matrix.smo, nrow(split.matrix.smo), ncol(split.matrix.smo))
	split.mat.grob <- levelplot(t(split.mat), axes=F, col.regions=m.pal, at=split.levels)
	split.mat.smo.grob <- levelplot(t(split.mat.smo), axes=F,col.regions=gl.pal, at=seq(-2.5,2.5,by=1))
	
	if(is.null(main)) main <- paste("chr", chrom, start, ":", end)
	
	if(length(device) !=0){
	if(!is.element("split.heatmaps", list.files())) dir.create("split.heatmaps")
			
	if(device == "quartz") {quartz(title=main,width=12, height=6)}
	if(device == "PNG")
		{png(file=paste("split.heatmaps/",project,".png",sep=""),width=1600, height=800, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("split.heatmaps/",project,".pdf",sep=""),width=12, height=6, pointsize=8)}
	if(device == "PS")
		{postscript(file=paste("split.heatmaps/",project,".ps",sep=""),width=10, height=5, pointsize=8)}}

	grid.newpage()
	
	if(is.null(cex.label)) cex.label <- -0.01*(nrow(amplicon))+1.2
	left.space <- right.space <- 2/ncol(split.mat)
	if(identical(plot.symbols,T) && nrow(split.mat) < 100)
		{
		if(is.null(symbol.width))
			{
			right.space <- -(nrow(split.mat)/500) +0.225
			}else
				{
				right.space <- symbol.width
				}
				
		}
	
	title <- viewport(x = 0.5, y = 0.95, w = 1, h = 0.1,name = "title")  
	pushViewport(title)
	grid.text(main, just="top")
	popViewport()

	matrix.vp <-
    viewport(width = 1, height = 0.8,
		layout = # necessary to fix aspect ratio
		grid.layout(3, 5,
			widths = c(0.5,left.space,0.5,right.space,0.5),
			heights = c(0.1,0.8,0.1),
			respect = TRUE),
			name = 'matrix.vp')

	pushViewport(matrix.vp)
	
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits))
    grid.text("aCGH")                  
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits))
    grid.text("Expression")                  
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 5,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits))
    grid.text("MWU")                  
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits))
	do.call("panel.levelplot", trellis.panelArgs(split.mat.smo.grob, 1))
	grid.rect(gp = gpar(col = 'black'))
	grid.lines(c(ncol(amp.matrix)/ncol(split.mat),ncol(amp.matrix)/ncol(split.mat)), c(0,1))
	grid.lines(c((ncol(amp.matrix)+1)/ncol(split.mat),(ncol(amp.matrix)+1)/ncol(split.mat)), c(0,1))
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits))
	do.call("panel.levelplot", trellis.panelArgs(split.mat.grob, 1))
	grid.rect(gp = gpar(col = 'black'))
	grid.lines(c(ncol(amp.matrix)/ncol(split.mat),ncol(amp.matrix)/ncol(split.mat)), c(0,1))
	grid.lines(c((ncol(amp.matrix)+1)/ncol(split.mat),(ncol(amp.matrix)+1)/ncol(split.mat)), c(0,1))
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits,, clip=T))
	for(i in 1:length(fData(amplicon)[,match(probeID, fvarLabels(exp.cgh))])){
		if(identical(plot.symbols,T) && nrow(split.mat) < 100) grid.text(rev(fData(amplicon)[,match(probeID, fvarLabels(exp.cgh))])[i],x=0.05, y=(i/length(fData(amplicon)[,match(probeID, fvarLabels(exp.cgh))])-0.5/length(fData(amplicon)[,match(probeID, fvarLabels(exp.cgh))])), just=c("left","centre"), gp=gpar(cex=cex.label))
		}
		
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 5,
                    xscale = split.mat.grob$x.limits,
                    yscale = split.mat.grob$y.limits))
	grid.rect(gp = gpar(col = 'black'))
	i <-1
	ybar <- 0
	while(i <= nrow(split.mat)){
		if(!is.na(rev(pvals)[i]) &&  rev(pvals)[i]< pval.thresh){bar.col="red"}else{(bar.col="blue")}
		grid.rect(0,ybar+(bar.width/2),bars[i],1/length(bars), just="left",
			gp = gpar(fill = bar.col, col="black"))
		grid.lines(x=line.thresh, gp=gpar(lty=2))
	i <- i +1
	ybar <- ybar + bar.width
	}
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits))
    grid.text("AMP", ncol(amp.matrix)/(2*ncol(split.matrix)))
    grid.text("NA", 1-(ncol(nonamp.matrix)/(2*ncol(split.matrix))))                     
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 3,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits))
	grid.text("AMP", ncol(amp.matrix)/(2*ncol(split.matrix)))
    grid.text("NA", 1-(ncol(nonamp.matrix)/(2*ncol(split.matrix))))                
	popViewport()
	
	pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 5,
                      xscale = split.mat.grob$x.limits,
                      yscale = split.mat.grob$y.limits))
    grid.text(paste("p.value threshold -",pval.thresh))                 
	popViewport()

	if(length(device) != 0 && device != "quartz") dev.off()

}




readGeneView <- function(design.file=NULL, TGS.threshold=2, f.data.file=NULL)
{
	require(Biobase)
	
	commands <- data.frame(commands.history="Analysis script: BACE", stringsAsFactors=F)
	commands <- rbind(commands,date())
	commands <- rbind(commands,deparse(match.call()))
	
	if(is.null(design.file)){design.file <- file.choose()}
	cat("Design file chosen was \"",basename(design.file),"\"\n", sep="")
	
	design <- read.delim(design.file, header=T, sep="\t", stringsAsFactors=F)
	design$file <- as.character(design$file)
	design$sampleNames <- as.character(design$sampleNames)

	commands <- rbind(commands,paste("Design file =", basename(design.file)))
	
	cat("Checking directory for results files...")
	directory.list <- list.files()
	missing <- design$file[which(!is.element(design$file,directory.list))]
	if(length(missing) != 0){
	cat(paste(missing, collapse="\n"),"\n")
	stop("These GeneView files in your design table are not in your directory\n")
	}else{cat("OK\n")}	
	
	geneview.list <- vector("list", nrow(design))
	
	for(i in 1:nrow(design)){
		geneview.list[[i]] <- read.table(file=design$file[i],sep="\t",header=T, check.names=F ,skip=1, stringsAsFactors=F)}
	
	TGS.list <- lapply(geneview.list, function(x)  x[,"gTotalGeneSignal"])
	TGS.matrix <- do.call(cbind, TGS.list)
	TGS.matrix[TGS.matrix <= 1] <- 1
	
	TGE.list <- lapply(geneview.list, function(x)  x[,"gTotalGeneError"])
	TGE.matrix <- do.call(cbind, TGE.list)
	
	flags.list <- lapply(geneview.list, function(x)  x[,"gIsGeneDetected"])
	flags <- do.call(cbind, flags.list)
	flags.matrix <- flags==0
	
	TGS.matrix.log <- log(TGS.matrix,2)
	TGS.matrix.log[TGS.matrix < TGS.threshold] <- NA
	flags.matrix[TGS.matrix < TGS.threshold] <- TRUE
	
	SystematicName.list[[7]][316] <- "banana"
	SystematicName.list <- lapply(geneview.list, function(x)  x[,"SystematicName"])
	if(sum(unlist(lapply(SystematicName.list, function(x) !identical(x, SystematicName.list[[1]])))) != 0)
		{stop("SystematicNames are not identical in results files")}else{
			featureNames.mirna <- make.names(SystematicName.list[[1]], unique=T)}
	
	sampleNames.mirna <- make.names(design$sampleNames, unique=T)
	
	colnames(TGS.matrix) <- colnames(TGS.matrix.log) <- colnames(TGE.matrix) <- colnames(flags.matrix) <- row.names(design) <- sampleNames.mirna
	row.names(TGS.matrix) <- row.names(TGS.matrix.log) <- row.names(TGE.matrix) <- row.names(flags.matrix) <- featureNames.mirna
	
	fData.mirna <- data.frame(probeID=featureNames.mirna)
	row.names(fData.mirna) <- featureNames.mirna
	eset <- new("ExpressionSet", exprs=TGS.matrix.log, featureData=new("AnnotatedDataFrame", data=fData.mirna))
	assayDataElement(eset, "error") <- TGE.matrix
	assayDataElement(eset, "flags") <- flags.matrix
	
	if(!is.null(f.data.file)){
		fdata <- read.table(file=fdata.file, header=T, sep="\t", stringsAsFactors=F)
		eset.fdata <- merge(fData(eset), fdata, by="probeID", all.x=T, all.y=F)
		eset.fdata <- eset.fdata[match(fData(eset)$probeID, eset.fdata$probeID, nomatch=0),]
		rownames(eset.fdata) <- featureNames(eset)
		eset.fdata.adf <- new("AnnotatedDataFrame", data=eset.fdata)
	
		cat(sum(is.element(fdata$probeID,fData(eset)$probeID)), "annotated probes\n")

  	  	varMetadata(eset.fdata.adf)[which(is.element(varLabels(eset.fdata.adf), fvarLabels(eset))),] <- fvarMetadata(eset)
    	varMetadata(eset.fdata.adf)$labelDescription[grep("chrom",names(eset.fdata))] <- "Chromosome"
		varMetadata(eset.fdata.adf)$labelDescription[grep("start",names(eset.fdata))] <- "Probe start position (bp)"
		varMetadata(eset.fdata.adf)$labelDescription[grep("end",names(eset.fdata))] <- "Probe end position (bp)"
		varMetadata(eset.fdata.adf)$labelDescription[grep("bac.id",names(eset.fdata))] <- "Probe BAC Clone ID"
		varMetadata(eset.fdata.adf)$labelDescription[grep("cytoband",names(eset.fdata))] <- "Cytoband"
		varMetadata(eset.fdata.adf)$labelDescription[grep("MB",names(eset.fdata))] <- "Mbp genome position"
		featureData(eset) <- eset.fdata.adf}
	
	
	if(validObject(eset)) cat ("valid ExpressionSet\n")
	eset
}


spitCIRCOS <- function(cgh, assayData.table="exprs"){
	if(!is.element("CIRCOS", list.files())) dir.create("CIRCOS")
	
	chroms <- fData(cgh)$chrom
	chroms <- paste(rep("hs", nrow(cgh)), chroms, sep="")
	chroms[chroms == "hs23"] <- "hsX"
	chroms[chroms == "hs24"] <- "hsY"
	for(i in 1:ncol(cgh)){
	circos.exprs <- paste(chroms, fData(cgh)$start, fData(cgh)$end, assayData(cgh)[[assayData.table]][,i])
	circos.gains <- circos.exprs[assayData(cgh)$GL[,i] >= 1]
	circos.losses <- circos.exprs[assayData(cgh)$GL[,i] <= -1]
	circos.amps <- circos.exprs[assayData(cgh)$GL[,i] > 1]
	circos.dels <- circos.exprs[assayData(cgh)$GL[,i] < -1]
	
	write.table(circos.exprs, paste("CIRCOS", paste(sampleNames(cgh)[i],"circos.aCGH.exprs.txt", sep="."),sep="/"), row.names=F, sep="\t", na="", quote=F, col.names=F)
	write.table(circos.amps, paste("CIRCOS", paste(sampleNames(cgh)[i],"circos.aCGH.amps.txt", sep="."),sep="/"), row.names=F, sep="\t", na="", quote=F, col.names=F)
	write.table(circos.dels, paste("CIRCOS", paste(sampleNames(cgh)[i],"circos.aCGH.dels.txt", sep="."),sep="/"), row.names=F, sep="\t", na="", quote=F, col.names=F)
	write.table(circos.gains, paste("CIRCOS", paste(sampleNames(cgh)[i],"circos.aCGH.gains.txt", sep="."),sep="/"), row.names=F, sep="\t", na="", quote=F, col.names=F)
	write.table(circos.losses, paste("CIRCOS", paste(sampleNames(cgh)[i],"circos.aCGH.losses.txt", sep="."),sep="/"), row.names=F, sep="\t", na="", quote=F, col.names=F)
	}
}


correlationMatrixPlot <- function(eset, assayDataTable= exprs(eset), dist.method="correlation", clust.method="ward", cor.method="pearson", colours=c("red", "green"), project, cor.dim="samples", fdata.ID=NULL, ptcex=0.5, main="Correlation plot"){
	
	require(Biobase)
	require(genefilter)
	require(latticeExtra)
	require(marray)
	if(is.null(fData.ID)) {feature.ID <- featureNames(eset)}else{
		feature.ID <- fData(eset)[,match(fdata.ID, names(fData(eset)))]}
		
	if(sum(is.na(assayDataTable)) != 0) stop("NAs in M table\nImpute values before clustering\n")
	if(cor.dim == "samples"|cor.dim == "columns"|cor.dim == 2) {
		M <- assayDataTable
		colnames(M) <- sampleNames(eset)}else{
			M <- t(assayDataTable)
			colnames(M) <- feature.ID}
	if(dist.method == "correlation"){
			atr <- hclust(dist(1-cor(M, method=cor.method)), method=clust.method)
			dd <- dist2(1-cor(M, method=cor.method))
			}else {
			atr <- hclust(dist(t(M), method=dist.method), method=clust.method)
			dd <- dist2(t(M))
			}
	
	diag(dd) <- 0
	dd.row <- as.dendrogram(atr)
	row.ord <- order.dendrogram(dd.row)
	
	legend <- list(top=list(fun=dendrogramGrob, args=list(x=dd.row, side="top")))
	
	mat.plot <- levelplot(dd[row.ord, row.ord], scales=list(x=list(rot=90), cex=ptcex), xlab="", ylab="", legend=legend, col.regions=maPalette(low = "red", high = "green", mid="white", k=100), main=main)
	
	pdf(paste(project, "matrix.correlation.pdf", sep="."))
	plot(mat.plot)
	dev.off()
	
	if(cor.dim == "samples"|cor.dim == "columns"|cor.dim == 2) {write.table(data.frame(sampleNames=sampleNames(eset), cor(M, method=cor.method)), paste(project, "correlations.xls", sep="."), sep="\t", row.names=F)}else{
		write.table(data.frame(probeID=feature.ID, cor(M, method=cor.method)), paste(project, "correlations.xls", sep="."), sep="\t", row.names=F)
		}
	
	
}

ExpressioncorrelationMulti <- function(eset, method="pearson", project){
	
	options(warn=-1)
	if(!is.element(paste(project, "expression.correlations", sep="."), list.files())) dir.create(paste(project, "expression.correlations", sep="."))

	correlation.log <- data.frame(fData(eset)$symbol, positively.correlated=rep(NA, nrow(fData(eset))), negatively.correlated=rep(NA, nrow(fData(eset))))
	
	for(i in 1:nrow(eset)){

	exp.cor <- exp.p <- rep(NA, length(fData(eset)))
	
	cat("Correlations for gene", i, fData(eset)$symbol[i], "\n")
	
	for(j in 1:nrow(eset)){
	cor <- cor.test(exprs(eset)[i,], exprs(eset)[j,], method=method, na.rm=T)
	exp.cor[j] <- cor$estimate
	exp.p[j] <- cor$p.value
	if(j %% 1000 == 0) cat(j," ")
	}
	
	exp.adjp <- p.adjust(exp.p, method="BH")
	cor.results <- cbind.data.frame(fData(eset), exp.cor, exp.p, exp.adjp)
	cor.results <- cor.results[which(cor.results$exp.p < 0.05),]
	cor.positive <- cor.results[cor.results$exp.cor >0,]
	cor.positive <- cor.positive[order(cor.positive$exp.p),]
	
	cor.negative <- cor.results[cor.results$exp.cor <0,]
	cor.negative <- cor.negative[order(cor.negative$exp.p),]
	
	positive.file <- paste(paste(project, "expression.correlations", sep="."),paste(fData(eset)$symbol[i], ".positively.correlated.xls",sep=""), sep="/")
	negative.file <- paste(paste(project, "expression.correlations", sep="."),paste(fData(eset)$symbol[i], ".negatively.correlated.xls",sep=""), sep="/")
	
	cat(c("Genes whose expression positively correlates with", fData(eset)$symbol[i], "Expression\n"), file=positive.file, append=F)
	write.table(cor.positive, sep="\t", na="", row.names=F, file=positive.file, append=T)
	
	cat(c("Genes whose expression negatively correlates with", fData(eset)$symbol[i], "Expression\n"), file=negative.file, append=F)
	write.table(cor.negative, sep="\t", na="", row.names=F, file=negative.file, append=T)
	
	correlation.log[i,2] <- nrow(cor.positive)
	correlation.log[i,3] <- nrow(cor.negative)

	}
	write.table(correlation.log, file=paste(paste(project, "expression.correlations", sep="."), "correlation.log.xls", sep="/"), row.names=F, na="", sep="\t")
}


waterfallExpression <- function(eset, pheno, pheno.colours=NULL, project, plot.names=NULL){
	
	pheno <- factor(pheno)
	pheno.end <- unlist(strsplit(deparse(match.call()$pheno),"\\$"))
	pheno.name <- pheno.end[length(pheno.end)]
	
	if(is.null(pheno.colours)) pheno.colours <- rep(brewer.pal(12,"Paired"),3)[1:(length(pheno))]
	
	if(!is.element(paste(project, "plots", sep="."), list.files())) dir.create(paste(project, "plots", sep="."))
	setwd(paste(project, "plots", sep="."))
	
	if(!is.element(paste(pheno.name, "boxplots", sep="."), list.files())) dir.create(paste(pheno.name, "boxplots", sep="."))
	if(!is.element(paste(pheno.name, "waterfalls", sep="."), list.files())) dir.create(paste(pheno.name, "waterfalls", sep="."))

for(i in 1:nrow(eset)){
	
	pdf(paste(paste(pheno.name, "boxplots", sep="."), "/",fData(eset)$symbol[i], ".boxplot.pdf", sep=""), width=12, height=12)
	boxplot(exprs(eset)[i,]~pheno, xlab=pheno.name, ylab="Z score", main= fData(eset)$symbol[i])
	legend("topleft", paste(levels(pheno), table(pheno)), bty="n")
	dev.off()
	
	pdf(paste(paste(pheno.name, "waterfalls", sep="."), "/",fData(eset)$symbol[i], ".waterfall.pdf", sep=""), width=12, height=6)
	barplot(exprs(eset)[i,order(exprs(eset)[i,])], col=pheno.colours[as.numeric(pheno[order(exprs(eset)[i,])])], main=fData(eset)$symbol[i], las=2, ylab="Expression", names.arg=plot.names)
	legend("topleft", levels(pheno), bty="n", fill=pheno.colours[1:length(levels(pheno))])
	dev.off()
	
	cat(i, "\t")
	}
	setwd("..")
}







#	Set class and definitions and methods the way you're supposed to!!!'


# Create a more robust class, with initialization and validation methods
# to ensure assayData contains specific matricies
#	Classes methods and definitions

setClass("BACE.cgh", contains="ExpressionSet", representation(thresholds="list", MADS="list", commands="character"))

setMethod("initialize", "BACE.cgh",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   exprs = new("matrix"),
                   commands = character(),
                   Fishers.groups = character(),
                   thresholds=list(),
                   ... ) {
            callNextMethod(.Object,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation,
                           exprs=exprs,
                           commands=commands,
                           thresholds=thresholds,
                           ...)
          })

setValidity("BACE.cgh", function(object) {
  assayDataValidMembers(assayData(object), c("exprs"))
})

#	promptClass() to write a shell of class documentation


setClass("BACE.exp", contains="ExpressionSet", representation(commands="character"))

setMethod("initialize", "BACE.exp",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   exprs = new("matrix"),
                   commands = character(),
                   ... ) {
            callNextMethod(.Object,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation,
                           exprs=exprs,
                           commands=commands,
                           ...)
          })

setValidity("BACE.exp", function(object) {
  assayDataValidMembers(assayData(object), c("exprs"))
})


setClass("BACE.exp.cgh", contains="ExpressionSet", representation(thresholds="list", MADS="list", commands="character"))

setMethod("initialize", "BACE.exp.cgh",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   exprs = new("matrix"),
                   commands = character(),
                   thresholds=list(),
                   ... ) {
            callNextMethod(.Object,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation,
                           exprs=exprs,
                           commands=commands,
                           thresholds=thresholds,
                           ...)
          })

setValidity("BACE.exp.cgh", function(object) {
  assayDataValidMembers(assayData(object), c("exprs"))
})

setClass("BACE.meth", contains="ExpressionSet", representation(methylated="matrix",unmethylated="matrix",pvals="matrix",commands="character"))

setMethod("initialize", "BACE.meth",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   commands = character(),
                   unmethylated = new("matrix"),
                   methylated = new("matrix"),
                   exprs = new("matrix"),
                   pvals = new("matrix"),
                   ... ) {
            callNextMethod(.Object,
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation,
                           unmethylated=unmethylated, methylated=methylated, exprs=exprs, pvals=pvals, commands=commands,
                           ...)
          })

setValidity("BACE.meth", function(object) {
  assayDataValidMembers(assayData(object), c("exprs"))
})

        
#	setClass(x, representation = slots, contains = inheritance)

#	setMethod("show" etc , class,
#	function(x){
#		body of function
#		)}
		
#	as
#	is
#	SetAs
#	class unions
#	setClassUnion() - how to use this??


#	extends()
#	setOldClass()

#	functions for initialize methods for classes
#	setMethod("initialize", "BACE.cgh", function(.Object, x,y,...){
#	body method
#	})
#setClass
#initialize
#new

#	validity

#make a function for valid object


#	validBACE.cgh <- function(cgh){
#	order(chrom, start) == 1:nrow(cgh)
#	}
	
#	setClass("BACE.aCGH", contains="ExpressionSet", representation(thresholds="list", commands="character"), validity=validBACE.cgh)

#methods are usually called BACE.plot, BACE.show etc




readBACE.cgh <- function(exprsFile, smoDataFile, featureDataFile, phenoDataFile, experimentDataFile, notesFile, path, annotation,
	exprsArgs = list(sep = sep, header = header, row.names = row.names, quote = quote),
	smoArgs = list(sep = sep, header = header, row.names = row.names, quote = quote),
	featureDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, stringsAsFactors = stringsAsFactors),
	phenoDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, stringsAsFactors = stringsAsFactors),
	experimentDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, stringsAsFactors = stringsAsFactors),
	sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE, row.names = 1L, warn=-1) 
{
    require(Biobase)
    if (missing(exprsFile)) 
        stop("exprs can not be missing!")
    exprsArgs$file = exprsFile
    ex = as.matrix(do.call(read.table, exprsArgs))
	
	if (!missing(phenoDataFile)) {
        phenoDataArgs$file = phenoDataFile
        pd = do.call(read.AnnotatedDataFrame, phenoDataArgs)
        if (!identical(sampleNames(pd), colnames(ex))) 
            stop("Column names of expression matrix must be identical to\n", 
                "the sample names of the phenodata table.\n", 
                "You could use 'options(error=recover)' to compare the", 
                "values of 'sampleNames(pd)' and 'colnames(ex)'.\n")
    }
    else {
        pd <- data.frame(sampleNames=colnames(ex), stringsAsFactors=F)
        row.names(pd) <- colnames(ex)
		pd <- new("AnnotatedDataFrame", data=pd)
     }
    
    obj = new("BACE.cgh", exprs = ex, phenoData = pd)
    
    if (!missing(smoDataFile)) {
        smoArgs$file = smoDataFile
        smo = as.matrix(do.call(read.table, smoArgs))
        if (!identical(sampleNames(pd), colnames(smo))) 
            stop("Column names of smo matrix must be identical to\n", 
                "the sample names of the phenodata table.\n", 
                "You could use 'options(error=recover)' to compare the", 
                "values of 'sampleNames(pd)' and 'colnames(ex)'.\n")
        if (!identical(row.names(ex), row.names(smo))) 
            stop("Row names of smo matrix must be identical to\n", 
                "the row names of the expression matrix.\n", 
                "You could use 'options(error=recover)' to compare the", 
                "values of 'sampleNames(pd)' and 'colnames(ex)'.\n")
        assayDataElement(obj, "smo") <- smo
    }
    
   if (!missing(featureDataFile)) {
        featureDataArgs$file = featureDataFile
        fd = do.call(read.AnnotatedDataFrame, featureDataArgs)
        if (!identical(featureNames(obj), row.names(fd@data))) 
            stop("Row names of featureData must be identical to\n", 
                "the row names of the expression matrix.\n", 
                "You could use 'options(error=recover)' to compare the", 
                "values of 'sampleNames(pd)' and 'colnames(ex)'.\n")
        
        varMetadata(fd)$labelDescription[grep("chrom", varLabels(fd))] <- "Chromosome"
		varMetadata(fd)$labelDescription[grep("start", varLabels(fd))] <- "Probe start position (bp)"
		varMetadata(fd)$labelDescription[grep("end", varLabels(fd))] <- "Probe end position (bp)"
		varMetadata(fd)$labelDescription[grep("bac.id", varLabels(fd))] <- "BAC clone idnetifier"
		varMetadata(fd)$labelDescription[grep("cytoband", varLabels(fd))] <- "Cytoband"
        
        featureData(obj) <- fd
        
    }
    
    if (!missing(experimentDataFile)) 
        experimentDataArgs$file = experimentDataFile
    if (!is.null(experimentDataArgs$file)) 
        experimentData(obj) <- do.call(read.MIAME, experimentDataArgs)
    if (!missing(annotation)) 
        annotation(obj) <- annotation
    if (!missing(notesFile)) 
        notes(obj) <- readLines(notesFile)
    validObject(obj)
    obj
}












cghHeatmap.Al <- function(cgh, cluster="GL", heatmap.matrix="GL", project, cex.labels=1, heatmap.scale=NULL, plot.symbols=T, main=NULL, device=NULL, dist.method="correlation", clust.method="ward", cor.method="pearson", phenotypes=NULL, pheno.colours=NULL, plot.sample.names=T, plot.matrix=T)
{
	require(Biobase)
	require(lattice)
	require(latticeExtra)
	require(grid)
	require(marray)
	
	cgh <- cgh[order(fData(cgh)$chrom, fData(cgh)$start),]
	
	heatmap.matrix.name <- heatmap.matrix
	heatmap.matrix <- assayDataElement(cgh, heatmap.matrix)
	
	if(!is.null(cluster) && !identical(cluster,F)){
		
		cluster <- assayDataElement(cgh, cluster)
		
		if(dist.method == "correlation"){
			atr <- hclust(dist(1-cor(cluster, method=cor.method)), method=clust.method)
			}else{atr <- hclust(dist(t(cluster), method=dist.method), method=clust.method)}
		atr.dg <- as.dendrogram(atr)
		col.order <- atr$order
		atr.grob <- dendrogramGrob(atr.dg, side="top", size=3)
		
	}else{col.order <- 1:ncol(cgh)}
	
	if(!is.null(phenotypes)){
		
		#	phenotypes
		pheno.labels <- as.matrix(phenotypes)
		
		pheno.matrix <- matrix(1, nrow(pheno.labels), ncol(pheno.labels))
	
		phenos <- unique(pheno.labels[!is.na(pheno.labels[,1]),1])
		phenos <- phenos[order(phenos)]
		if(ncol(pheno.labels) > 1)
		{
			for(i in 2:ncol(pheno.labels))
			{
			phenos.i <- unique(pheno.labels[!is.na(pheno.labels[,i]),i])
			phenos.i <- phenos.i[order(phenos.i)]
			phenos <- c(phenos,phenos.i[which(!is.element(phenos.i,phenos))])
			}
		}
	
		pheno.matrix <- matrix(as.numeric(match(pheno.labels, phenos)),nrow(pheno.labels), ncol(pheno.labels))
	
		pheno.matrix <- pheno.matrix[col.order,,drop=F]
		pheno.labels <- pheno.labels[col.order,,drop=F]
	
		leg.colours <- unique(as.vector(pheno.matrix[!is.na(pheno.matrix)]))
		leg.labels <- unique(as.vector(pheno.labels[!is.na(pheno.labels)]))
	
	
		if(is.null(pheno.colours)) {leg.pal <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Set2"), brewer.pal(8,"Dark2"))[1:(length(leg.colours))]}else{
		leg.pal <- pheno.colours}
		leg.pal <- leg.pal[1:length(leg.colours)]	
	
		if(ncol(pheno.matrix) ==1 ){pheno.matrix <- cbind(pheno.matrix, pheno.matrix)}
		pheno.mat.grob <- levelplot(pheno.matrix, col.regions=leg.pal, at=0:length(leg.colours)+0.5,colorkey=F)
	
		# reverse legend labels - ie alphabetical downwards
		leg.mat <- as.matrix(leg.colours[rev(order(leg.colours))],length(leg.colours),1)
		leg.grob <- levelplot(t(leg.mat), col.regions=leg.pal, at=0:length(leg.colours)+0.5, colorkey=F)
		leg.labels <- leg.labels[rev(order(leg.colours))]
		
		}
	
	heat.matrix <- heatmap.matrix[rev(1:nrow(cgh)),col.order]
	maxM <- (as.integer(max(abs(heat.matrix[!is.na(heat.matrix)]))))+1
	
	if(heatmap.matrix.name == "GL") heatmap.scale <- 2
	if(is.null(heatmap.scale)) heatmap.scale <- min(abs(range(heat.matrix)))
	
	if(maxM > heatmap.scale) {
		heat.matrix[heat.matrix > heatmap.scale] <- heatmap.scale
		heat.matrix[heat.matrix < -heatmap.scale] <- -heatmap.scale
		}

	if(heatmap.matrix.name == "GL"){
		split.levels=seq(1:6)-3.5
		m.pal <- maPalette(low = "blue", mid="white", high = "red", k =6)
		}else{
			split.levels <- seq(min(heat.matrix, na.rm=T),max(heat.matrix, na.rm=T),by=0.1)
			m.pal <- maPalette(low = "blue", mid="white", high = "red", k =length(split.levels))
		}
	 
	heat.mat <- as.matrix(heat.matrix, nrow(heat.matrix), ncol(heat.matrix))
	heat.mat.grob <- levelplot(t(heat.mat), axes=F, col.regions=m.pal, at=split.levels)
	s.labels <- sampleNames(cgh)[col.order]
	
	colorkey.mat <- matrix(split.levels, nrow=2, ncol=length(split.levels), byrow=T)
	colorkey.grob <- levelplot(colorkey.mat[,-1], axes=F, col.regions=m.pal, at=split.levels)
	
	if(length(device) !=0){
	if(!is.element("heatmaps", list.files())) dir.create("heatmaps")
			
	if(device == "quartz") {quartz(title=main,width=4, height=6)}
	if(device == "PNG")
		{png(file=paste("heatmaps/",project,".heatmap.png",sep=""),width=600, height=900, pointsize=12)}
	if(device == "JPEG")
		{jpeg(file=paste("heatmaps/",project,".heatmap.jpeg",sep=""),width=600, height=900, pointsize=12)}
	if(device == "PDF")
		{pdf(file=paste("heatmaps/",project,".heatmap.pdf",sep=""),width=6, height=9, pointsize=12)}
	if(device == "PS")
		{postscript(file=paste("heatmaps/",project,".heatmap.ps",sep=""),width=6, height=9, pointsize=12)}}
	
	grid.newpage()
	
	title.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.9,  "npc"),
	width = unit(1, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "bottom"))
	pushViewport(title.vp)
	grid.text(main, just="top", gp=gpar(cex=1))
	popViewport()
	
	colorkey.vp <- viewport(x = unit(0.825, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.025,  "npc"),
	height = unit(0.1,  "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,(ncol(colorkey.mat)/2)-0.5),
	yscale=c(0.5,ncol(colorkey.mat)-0.5))
	pushViewport(colorkey.vp)
	do.call("panel.levelplot", trellis.panelArgs(colorkey.grob, 1))
	popViewport()
	
	colorkey.leg.vp <- viewport(x = unit(0.85, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.1,  "npc"),
	height = unit(0.1,  "npc"),
	just = c("left", "bottom"))
	pushViewport(colorkey.leg.vp)
	
		if(heatmap.matrix.name == "GL"){
			glad.leg <- c("Del", "Loss", "NC", "Gain", "Amp")
			for(i in 1:5) grid.text(glad.leg[i],y=(i/5-0.5/5),
		x=0.05, just=c("left","centre"), gp=gpar(cex=min(0.5,15/5)))
			}else{
				for(i in 1:length(split.levels)){
			grid.text(split.levels[i],y=(i/length(split.levels)-0.5/length(split.levels)),
		x=0.05, just=c("left","centre"), gp=gpar(cex=min(0.5,15/length(split.levels))))
			}
		}
	popViewport()
	
	if(!is.null(phenotypes)) {
		phenobar.height <- 0.05
		
		
		if(ncol(pheno.labels) ==1 ) {phenobar.legend.height <- phenobar.height/2*ncol(pheno.labels)}else{
			phenobar.legend.height <- phenobar.height/ncol(pheno.labels)}
			
		phenobar.legend <-
   			viewport(x=0 + phenobar.legend.height/2, y = unit(0.8-phenobar.height+phenobar.legend.height/2,  "npc"),
   			width = 1, height = unit(phenobar.legend.height*length(leg.labels), "npc"),
			layout = # necessary to fix aspect ratio
			grid.layout(1, 4,
			widths = unit(c(0.2,0.6,phenobar.legend.height, 0.2-phenobar.legend.height),"npc"),
			heights = rep(phenobar.legend.height*length(leg.labels),3),
			respect = TRUE),
			just = c("left", "bottom"),
			name = 'phenobar.legend')
		
		
		
		pushViewport(phenobar.legend)
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3,
                      xscale = leg.grob$x.limits,
                      yscale = leg.grob$y.limits))
		do.call("panel.levelplot", trellis.panelArgs(leg.grob, 1))
		popViewport()
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4,
                      xscale = leg.grob$x.limits,
                      yscale = leg.grob$y.limits))
		
		leg.fontsize <- min(8,(300/length(leg.labels)))
		for(i in 1:length(leg.labels)){
		grid.text(leg.labels[i],x=0.1, y=(i/length(leg.labels)-0.5/length(leg.labels)), 		just="left",gp=gpar(fontsize=leg.fontsize))}
		popViewport()
		
		popViewport()

		
		phenobar.vp <-
   		viewport(x=0, y = unit(0.8-phenobar.height,  "npc"),
   		width = 1, height = phenobar.height,
		layout = # necessary to fix aspect ratio
		grid.layout(1, 3,
			widths = unit(c(0.2,0.6,0.2),"npc"),
			heights = rep(phenobar.height,3),
			respect = TRUE),
			just = c("left", "bottom"),
			name = 'phenobar.vp',
			xscale = pheno.mat.grob$x.limits,
            yscale = pheno.mat.grob$y.limits)
		
			
		pushViewport(phenobar.vp)
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
		
		phenotypes.fontsize <- min(8,(300/ncol(pheno.labels)))
		for(i in 1:ncol(pheno.labels)){
		grid.text(names(phenotypes)[i],x=0.9, y=(i/ncol(pheno.labels)-0.5/ncol(pheno.labels)), 		just=c("right","centre"),gp=gpar(fontsize=phenotypes.fontsize))}
		popViewport()
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2,
                      xscale = pheno.mat.grob$x.limits,
                      yscale = pheno.mat.grob$y.limits))
   		grid.rect(gp = gpar(col = 'black', lwd=1))
   		do.call("panel.levelplot", trellis.panelArgs(pheno.mat.grob, 1))
		
	#	plot horoizontal lines
	#	line.breaks <- seq(0,1,by=1/ncol(pheno.labels))
	#	for(i in 1:length(line.breaks)) grid.lines(y=line.breaks[i], gp = gpar(col = 'black', lwd=1))

		popViewport()
		popViewport()
		
				
		}else{phenobar.height <- 0}
		
		
	symbol.vp <- viewport(x = unit(0, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.2, "npc"),
	height = unit(0.7-phenobar.height, "npc"),
	just = c("left", "bottom"))
	pushViewport(symbol.vp)
	
			chrom.lengths <- rev(as.vector(table(fData(cgh)$chrom))/length(fData(cgh)$chrom))
			i <- 1
			chrom.start <- 0
			while(i <= length(chrom.lengths)){
			grid.rect(0.7,chrom.start,width=0.2,height=chrom.lengths[i], just=c("centre", "bottom"), gp=gpar(fill=i%%2))
			chroms <- unique(fData(cgh)$chrom)
			chroms[chroms == 23] <- "X"
			chroms[chroms == 24] <- "Y"
			grid.text(rev(chroms)[i], 0.5,chrom.start+(chrom.lengths[i]/2), gp=gpar(cex=min(1,12/length(chroms))))
			chrom.start <- chrom.start+chrom.lengths[i]
			i <- i+1
			}
	popViewport()
	
	if(plot.sample.names==T){
	sample.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(0.1, "npc"),
	just = c("left", "top"))
	pushViewport(sample.vp)
	for(i in 1:length(s.labels)){
		grid.text(s.labels[i],x=(i/length(s.labels)-0.5/length(s.labels)),
		y=0.95, just=c("right","centre"), gp=gpar(cex=cex.labels*min(1,30/length(s.labels))), rot=90)}
	popViewport()
	}
	
	matrix.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.1,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(0.7-phenobar.height, "npc"),
	just = c("left", "bottom"),
	xscale=c(0.5,ncol(heat.mat)+0.5),
	yscale=c(0.5,nrow(heat.mat)+0.5))
	pushViewport(matrix.vp)
	grid.rect()
	if(plot.matrix == T) do.call("panel.levelplot", trellis.panelArgs(heat.mat.grob, 1))
	popViewport()
	
	if(!is.null(cluster) && !identical(cluster,F)){
	
	atr.vp <- viewport(x = unit(0.2, "npc"),
	y = unit(0.8,  "npc"),
	width = unit(0.6, "npc"),
	height = unit(3, "lines"),
	just = c("left", "bottom"))
	pushViewport(atr.vp)
	grid.draw(atr.grob)
	popViewport()

	}
	if(length(device) != 0 && device != "quartz") dev.off()
}



