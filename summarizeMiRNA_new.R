#summarize_mature_miRNA <- function(data, miRNA_all_loc, genome, chr_length=NULL, offset=0, verbose=F) {
#	counts <- numeric()
#	leftover <- AlignedRead()
#	leftover@alignData@varMetadata <- data@alignData@varMetadata
#	
#	miRNA_all_loc <- miRNA_all_loc[order(miRNA_all_loc[,1], miRNA_all_loc[,2], miRNA_all_loc[,3]),]
#	read_range <- RangedData(IRanges(start=position(data), width=width(data)), space=chromosome(data))
#	dum <- findOverlaps(read_range, ann_range, type="within")
#	counts <-  table(as.matrix(dum)[,2])
#	names <- miRNA_all_loc[as.numeric(names(counts)),"name"]
#	counts_new <- aggregate(as.numeric(counts), list(as.character(names)), sum)
#	
#	summary_count <- array(0, dim=length(levels(miRNA_all_loc$name)))
#	names(summary_count) <- levels(miRNA_all_loc$name)
#	ind <- match(counts_new[,1], names(summary_count))
#	summary_count[ind] <- counts_new[,2]
#
##	if(length(unique(as.matrix(dum)[,1])) >0) {
	#	leftover <- data[-unique(as.matrix(dum)[,1])]
#	} else {
#		leftover <- data
#	}
#	gc()
#	return(list(summary_count, leftover))
#}#
#}#

summarize_mature_miRNA<- function(data_new, miRNA_all_loc, genome, chr_length=NULL, offset=0, verbose=F) {
	counts <- numeric()
	#leftover <- AlignedRead()
	#leftover@alignData@varMetadata <- data_new@alignData@varMetadata
	#browser()
	miRNA_all_loc <- miRNA_all_loc[order(miRNA_all_loc[,1], miRNA_all_loc[,2], miRNA_all_loc[,3]),]
	#read_range <- RangedData(IRanges(start=start(data_new), width=width(data_new)), space=as.character(rname(data_new)))
	#Now using GappedAligment instead of shortread
	ann_range <- RangedData(IRanges(start=miRNA_all_loc$start - offset, end=miRNA_all_loc$end + offset), space=miRNA_all_loc$chr, name=miRNA_all_loc$name)
	dum <- findOverlaps(data_new, ann_range, type="within")
	counts <-  table(subjectHits(dum))
	names <- miRNA_all_loc[as.numeric(names(counts)),"name"]
	counts_new <- aggregate(as.numeric(counts), list(as.character(names)), sum)
	
	summary_count <- array(0, dim=length(levels(miRNA_all_loc$name)))
	names(summary_count) <- levels(miRNA_all_loc$name)
	ind <- match(counts_new[,1], names(summary_count))
	summary_count[ind] <- counts_new[,2]

	if(length(unique(queryHits(dum))) >0) {
		leftover <- data_new[-unique(queryHits(dum)),]
	} else {
		leftover <- data_new
	}
	gc()
	return(list(summary_count, leftover))
}
