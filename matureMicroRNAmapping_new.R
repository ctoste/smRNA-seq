matureMicroRNAmapping_new <- function(fastaFile, genome, type=c("matureMiRNA", "matureStarMiRNA")) {
###Improved speed by changing datastructure from matrix (rbind) to list

	#library(time)
	library(ShortRead)
	library(parallel)
	browser()
	type <- match.arg(type)
	if (!is.element(type, c("matureMiRNA", "matureStarMiRNA")))
		stop ("'type' need to be 'matureMiRNA', or 'matureStarMiRNA'")
	pckName <- available.packages(contriburl = contrib.url("http://www.bioconductor.org/packages/release/data/annotation"))[,1]
	pckName.BSgenome <- pckName[grep("BSgenome",pckName)]
	libName <- as.character(pckName.BSgenome[grep(paste(genome, "$", sep=""), pckName.BSgenome)])

	if (is.na(libName)) {
		stop("Genome not supported")
	} 
	if (length(libName)>1) {
		stop("library name is not unique")
	} 
	
	require(libName, character.only = TRUE)
	seqnames <- seqnames(get(strsplit(libName, "\\.")[[1]][2]))
	dum <- readRNAStringSet(file=fastaFile, "fasta") 
	ind <- grep(tolower(substr(strsplit(libName, "\\.")[[1]][2], 1,3)), names(dum))
	dict0 <- DNAStringSet(dum[ind])

	count <- 0
	result_all <- mclapply(seqnames, function(seqname) {
	#for (seqname in seqnames) {
		#browser()
		subject <- get(strsplit(libName, "\\.")[[1]][2])[[seqname]]
		#result[[seqname]] <- list()
		result <- list()
		print(seqname)
		
		#prev <- progressBar()
		for(j in seq_len(length(dict0))) {
			#print(j)
			patternID <- strsplit(names(dict0)[j], " ")[[1]][1]
			pattern <- dict0[[j]]
			#index <- match(seqname[j], seqnames(genome))
			dum <- matchPattern(pattern, subject)
			if(length(dum) !=0) {
				#browser()
					count <- count + 1
					result[[count]] <- data.frame(chr = rep(seqname, length(dum)), start = start(dum), 
						end = end(dum), strand = rep("+", length(dum)), name = rep(patternID, length(dum)), type=rep(type, length(dum)))

			}
			rcpattern <- reverseComplement(pattern)
			dum1<- matchPattern(rcpattern, subject)
			if(length(dum1) !=0) {
					count <- count + 1
					result[[count]] <- data.frame(chr = rep(seqname, length(dum1)), start = start(dum1), 
						end = end(dum1), strand = rep("-", length(dum1)), name = rep(patternID, length(dum1)), type=rep(type, length(dum1)))
			}
			#prev <- progressBar(j/length(dict0), prev)
		}
		return(result)
	}, mc.cores=12)
	#browser()
	a <- do.call("rbind", lapply(result_all, function(X) do.call("rbind", X)))
	return(a)
}
#dum <- matureMicroRNAmapping_new("mature.fa", genome="hg19")

