preprocess_miRNA <- function (data_dir = getwd(), miRNA_file, RNA_repeat_file, genome = "hg18", offset = 5, format = "bam") {
    library(rtracklayer)
	#browser()
    datadir <- system.file("extdata", package = "smRNAseq3")
    annotation_new <- function(reads, ann) {
        reads_pos <- summarize_mature_miRNA(reads[as.character(strand(reads)) == 
            "+"], ann[ann[, 4] == "+", ], genome, offset = offset)
        reads_neg <- summarize_mature_miRNA(reads[as.character(strand(reads)) == 
            "-"], ann[ann[, 4] == "-", ], genome, offset = offset)
        reads_all <- rbind(cbind(names(reads_pos[[1]]), reads_pos[[1]]), 
            cbind(names(reads_neg[[1]]), reads_neg[[1]]))
        reads_mature_count <- tapply(as.numeric(reads_all[, 2]), 
            reads_all[, 1], sum, na.rm = T)
        return(list(summary_count = reads_mature_count, leftover = c(reads_pos[[2]], 
            reads_neg[[2]])))
    }
    print("Obtaining annotation information from UCSC...")
    miRNA_all <- read.delim(file = miRNA_file)
	RNAgene <- read.delim(file=RNA_repeat_file, skip=1)
    if (genome == "hg18") {
        session <- browserSession()
        genome(session) <- genome
        foo <- getTable(ucscTableQuery(session, trackNames(session)[grep("wgRna", trackNames(session))]))
        miRNA_gene <- foo[foo$type == "miRNA", c(2, 3, 4, 7, 5, 10)]
        miRNA_gene$type <- "miRNA_gene"
        colnames(miRNA_gene) <- colnames(miRNA_all)
        piRNA_aligned <- readAligned("/net/isi-solexa/ifs/bic/xwu/annotation/hg18/piRNA_bowtie_aligned.bam", 
            type = "BAM")
        piRNA_aligned <- piRNA_aligned[!is.na(chromosome(piRNA_aligned))]
        piRNA <- data.frame(chr = chromosome(piRNA_aligned), 
            start = position(piRNA_aligned), end = position(piRNA_aligned) + 
                width(piRNA_aligned), strand = strand(piRNA_aligned), 
            name = as.character(piRNA_aligned@id), type = "piRNA")
        #RNAgene <- read.delim(file = gzfile(paste(datadir, "hg18_RNA_gene.txt.gz", sep = "/")))
		#RNAgenenew <- data.frame(RNAgene[RNAgene$repClass != "miRNA", c(6:8, 10:12)])
		RNAgenenew <- data.frame(RNAgene[RNAgene$type != "miRNA", c(1:3, 6,4,8)])
        colnames(RNAgenenew) <- colnames(miRNA_all)
        RNAgenenew_all <- rbind(miRNA_gene, piRNA, RNAgenenew)
        refseq <- getTable(ucscTableQuery(session, trackNames(session)[grep("refGene", 
            trackNames(session))]))
        refseq_NM <- refseq[grep("NM", refseq$name), ]
        tx_range <- RangedData(IRanges(start = refseq_NM$txStart, 
            end = refseq_NM$txEnd), space = refseq_NM$chrom, 
            refseq_NM[, -c(1, 3, 5, 6)])
    }
    else if (genome == "hg19") {
        session <- browserSession()
        genome(session) <- genome
        sno_miRNA <- getTable(ucscTableQuery(session, trackNames(session)[grep("wgRna", 
            trackNames(session))]))
        dum <- sno_miRNA[sno_miRNA$type != "miRNA", c(2, 3, 4, 7, 5, 10)]
        dum$type = "snoRNA"
        miRNA_gene <- read.delim("/isi-solexa/bic/xwu/annotation/mirBase_v18/hg19_hsa.gff3", 
            comment.char = "#", header = F)
        foo <- strsplit(as.character(miRNA_gene[, 9]), split = c(";|="))
        miRNA_gene <- cbind(paste("chr", miRNA_gene[, 1], sep = ""), 
            miRNA_gene[, c(4, 5, 7)], sapply(foo, function(X) X[4]), 
            miRNA_gene[, 3])
        colnames(miRNA_gene) <- colnames(miRNA_all)
        colnames(dum) <- colnames(miRNA_all)
        sno_miRNA_new <- rbind(dum, miRNA_gene)
       # RNAgene <- read.delim(file = "/net/isi-solexa/ifs/bic/xwu/annotation/hg19/hg19_RNA_repeat_053111.txt", skip = 1)
        RNAgenenew <- data.frame(RNAgene[, c(6:8, 10:12)])
        colnames(RNAgenenew) <- colnames(miRNA_all)
        RNAgenenew_all <- rbind(sno_miRNA_new, RNAgenenew)
        refseq <- getTable(ucscTableQuery(session, trackNames(session)[grep("refGene", 
            trackNames(session))]))
        refseq_NM <- refseq[grep("NM", refseq$name), ]
        tx_range <- RangedData(IRanges(start = refseq_NM$txStart, 
            end = refseq_NM$txEnd), space = refseq_NM$chrom, 
            refseq_NM[, -c(1, 3, 5, 6)])
    }
    else if (genome == "mm9") {
        session <- browserSession()
        genome(session) <- genome
        miRNA_gene <- getTable(ucscTableQuery(session, trackNames(session)[grep("miRNA", 
            trackNames(session))]))
        miRNA_gene$type <- "miRNA_gene"
        miRNA_gene_new <- miRNA_gene[, c(2, 3, 4, 7, 5, 8)]
        colnames(miRNA_gene_new) <- colnames(miRNA_all)
        #RNAgene <- read.delim(file = "/net/isi-solexa/ifs/bic/xwu/annotation/mm9/mm9_RNA_gene.txt", skip = 1)
        RNAgenenew <- data.frame(RNAgene[, c(6:8, 10:12)])
        colnames(RNAgenenew) <- colnames(miRNA_all)
        RNAgenenew_all <- rbind(miRNA_gene_new, RNAgenenew)
        refseq <- getTable(ucscTableQuery(session, trackNames(session)[grep("refGene", 
            trackNames(session))]))
        refseq_NM <- refseq[grep("NM", refseq$name), ]
        tx_range <- RangedData(IRanges(start = refseq_NM$txStart, 
            end = refseq_NM$txEnd), space = refseq_NM$chrom, 
            refseq_NM[, -c(1, 3, 5, 6)])
    }
	else if (genome == "mm10") {
        session <- browserSession()
        genome(session) <- genome
        #miRNA_gene <- getTable(ucscTableQuery(session, trackNames(session)[grep("miRNA", 
        #    trackNames(session))]))
        #miRNA_gene$type <- "miRNA_gene"
        #miRNA_gene_new <- miRNA_gene[, c(2, 3, 4, 7, 5, 8)]
        #colnames(miRNA_gene_new) <- colnames(miRNA_all)
        #RNAgene <- read.delim(file = "/net/isi-solexa/ifs/bic/xwu/annotation/mm9/mm9_RNA_gene.txt", skip = 1)
        RNAgenenew <- data.frame(RNAgene[, c(6:8, 10:12)])
        colnames(RNAgenenew) <- colnames(miRNA_all)
        RNAgenenew_all <- rbind(miRNA_all, RNAgenenew)
        refseq <- getTable(ucscTableQuery(session, trackNames(session)[grep("refGene", 
            trackNames(session))]))
        refseq_NM <- refseq[grep("NM", refseq$name), ]
        tx_range <- RangedData(IRanges(start = refseq_NM$txStart, 
            end = refseq_NM$txEnd), space = refseq_NM$chrom, 
            refseq_NM[, -c(1, 3, 5, 6)])
    }
	else if (genome == "rn4") {
        session <- browserSession()
        genome(session) <- genome
        miRNA_gene <- getTable(ucscTableQuery(session, trackNames(session)[grep("miRNA", 
            trackNames(session))]))
        miRNA_gene$type <- "miRNA_gene"
        miRNA_gene_new <- miRNA_gene[, c(2, 3, 4, 7, 5, 10)]
        colnames(miRNA_gene_new) <- colnames(miRNA_all)
        #RNAgene <- read.delim(file = "/net/isi-solexa/ifs/bic/xwu/annotation/rn4/rn4_RNA_repeat.txt", skip = 1)
        RNAgenenew <- data.frame(RNAgene[, c(6:8, 10:12)])
        colnames(RNAgenenew) <- colnames(miRNA_all)
        RNAgenenew_all <- rbind(miRNA_gene_new, RNAgenenew)
        refseq <- getTable(ucscTableQuery(session, trackNames(session)[grep("refGene", 
            trackNames(session))]))
        refseq_NM <- refseq[grep("NM", refseq$name), ]
        tx_range <- RangedData(IRanges(start = refseq_NM$txStart, 
            end = refseq_NM$txEnd), space = refseq_NM$chrom, 
            refseq_NM[, -c(1, 3, 5, 6)])
    }
	
	else if (genome == "galGal4") {
        session <- browserSession()
        genome(session) <- genome
        miRNA_gene <- read.delim("/net/isi-solexa/ifs/bic/xwu/annotation/mirBase_v19/gga.gff3", 
            comment.char = "#", header = F)
        foo <- strsplit(as.character(miRNA_gene[, 9]), split = c(";|="))
        miRNA_gene <- cbind(paste("chr", miRNA_gene[, 1], sep = ""), 
            miRNA_gene[, c(4, 5, 7)], sapply(foo, function(X) X[4]), 
            miRNA_gene[, 3])
        colnames(miRNA_gene) <- colnames(miRNA_all)
		
        #RNAgene <- read.delim(file = "/net/isi-solexa/ifs/bic/xwu/annotation/ggal4/ggal4_RNA_gene.txt", skip = 1)
        RNAgenenew <- data.frame(RNAgene[, c(6:8, 10:12)])
        colnames(RNAgenenew) <- colnames(miRNA_all)
        RNAgenenew_all <- rbind(miRNA_gene, RNAgenenew)
        refseq <- read.delim("/isi-solexa/bic/xwu/annotation/ggal4/ggal4_refseq_gene_060414.txt")
        refseq_NM <- refseq[grep("NM", refseq$name), ]
        tx_range <- RangedData(IRanges(start = refseq_NM$txStart, 
            end = refseq_NM$txEnd), space = refseq_NM$chrom, 
            refseq_NM[, -c(1, 3, 5, 6)])
    }
    print("Done")
	browser()
	fnames <- list.files(data_dir, paste(format, "$", sep = ""))
    name <- sapply(fnames, function(X) strsplit(X, "\\.")[[1]][1])
    count_summary <- array(dim = c(length(unique(RNAgenenew_all$type)) + 4, length(fnames)))
    miRNA_table <- integer()
    ncRNA_table <- integer()
    setwd(data_dir)
    for (i in 1:length(fnames)) {
        print(paste("Processing", fnames[i]))
		data <- readBamGappedAlignments(fnames[i])
        ##################
		#Sometimes, CIGAR had problem, so above won't run
		##################
		#data <- readAligned(fnames[i], , type="BAM")
		#data <- as(data, "GRanges")
		#browser()
		print("Finished Reading Sequences.")
        print("counting mature miRNAs...")
        data_pos <- summarize_mature_miRNA(data[as.character(strand(data)) == "+"], miRNA_all[miRNA_all$strand == "+", ], 
						genome, offset = offset)
        data_neg <- summarize_mature_miRNA(data[as.character(strand(data)) == "-"], miRNA_all[miRNA_all$strand == "-", ], 
						genome, offset = offset)
		data_all <- rbind(cbind(names(data_pos[[1]]), data_pos[[1]]), cbind(names(data_neg[[1]]), data_neg[[1]]))
        data_mature_count <- tapply(as.numeric(data_all[, 2]), data_all[, 1], sum, na.rm = T)
        print("Done")
        
		print("counting other ncRNAs...")
        data_RNA_gene <- annotation_new(c(data_pos[[2]], data_neg[[2]]), RNAgenenew_all)
        ind <- match(names(data_RNA_gene[[1]]), RNAgenenew_all$name)
        type <- RNAgenenew_all$type[ind]
        RNA_gene_summary <- tapply(data_RNA_gene[[1]], type, sum)
        print("Done")
        #browser()
		print("counting refseq RNAs")
        data_range <- RangedData(IRanges(start = start(data_RNA_gene[[2]]), 
            width = width(data_RNA_gene[[2]])), space = as.character(seqnames(data_RNA_gene[[2]])))
        data_in_refseq <- findOverlaps(data_range, tx_range)
        refseq_summary <- length(unique(queryHits(data_in_refseq)))
        intergenic <- length(data_RNA_gene[[2]]) - refseq_summary
        print("Done")
		#browser()
        count_summary[, i] <- c(miRNA = length(data) - length(data_pos[[2]]) - 
            length(data_neg[[2]]), RNA_gene_summary, Refseq = refseq_summary, 
            intergenic, Aligned_Total = length(data))
        miRNA_table <- cbind(miRNA_table, data_mature_count)
        ncRNA_table <- cbind(ncRNA_table, data_RNA_gene[[1]])
    }
    colnames(count_summary) <- name
    rownames(count_summary) <- c("Mature miRNA", names(RNA_gene_summary), 
        "Refseq", "Intergenic", "Total Aligned")
    colnames(miRNA_table) <- name
    colnames(ncRNA_table) <- name
    return(list(count_summary, miRNA_table, ncRNA_table))
}
