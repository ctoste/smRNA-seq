getAnnotation <- function(TC, ann, extend=10) {
	result <- character()
	for(i in levels(TC[,1])) {
		#print(i)
		#browser()
		TC_chr <- TC[TC[,1] == i,]
		ann_chr <- ann[ann[,1] == i,]	
		TC_range <- IRanges(start=as.numeric(as.character(TC_chr[,2])), end=as.numeric(as.character(TC_chr[,3])))
		ann_range <- IRanges(start=as.numeric(as.character(ann_chr[,2]))-extend, end=as.numeric(as.character(ann_chr[,3]))+extend)
		dum <- findWithin(TC_range, ann_range)
		if (length(dum) == 0) {
			next
		}
		#print(paste("Found", nrow(dum)))
		result <- rbind(result, cbind(TC_chr[as.numeric(as.matrix(dum)[,1]),], ann_chr[as.numeric(as.matrix(dum)[,2]),]))
	}
	return(result)
}

#foo <- getAnnotation(CNS_TC_all_2, miRNA_all, 10)
