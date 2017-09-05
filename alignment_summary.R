Alignment_summary <- function(data_dir, format, ext, out) {
	fnames <- list.files(data_dir, ext)
	result <- array(dim=c(4,length(fnames)))
	if (format=="novoalign") {
		system(paste("tail *", ext, " > temp.txt", sep=""))
		dum <- read.delim("temp.txt", header=F)
		ind <- grep("==", dum[,1])
		for (i in 1:length(ind)) {
			fname <-  strsplit(as.character(dum[ind[i],1]), split=" ")[[1]][2]
			total_read <- as.integer(strsplit(as.character(dum[ind[i]+3,1]), ":")[[1]][2])
			aligned_read <- as.integer(strsplit(as.character(dum[ind[i]+4,1]), ":")[[1]][2])
			filtered <- as.integer(strsplit(as.character(dum[ind[i]+7,1]), ":")[[1]][2]) + as.integer(strsplit(as.character(dum[ind[i]+8,1]), ":")[[1]][2])
			aligned_per <- aligned_read/(total_read-filtered)
			result[,i] <- c(total_read, total_read - filtered, aligned_read, aligned_per)
		}
		
	} else if (format=="bowtie") {
		
	}	
	colnames(result) <- sapply(fnames, function(X) strsplit(X, split="\\.")[[1]][1])
	rownames(result) <- c("Total Read", "Reads Passed Filtered", "Aligned Read", "Aligned %")
	write.table(result, out, sep="\t", quote=F)
	return(result)
}
#Alignment_summary(".", "novoalign", ".txt", "aligned_TTTTT_new.txt")
