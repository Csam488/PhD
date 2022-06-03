#! /usr/bin/env Rscript
#installs optparse if not already installed
suppressWarnings(suppressPackageStartupMessages(library(optparse)))
suppressPackageStartupMessages(library(StructuralVariantAnnotation))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
  make_option(c("-m", "--mode"), type="character", default="BEDPE",
              help="specify output mode: BED, BEDPE, or BOTH [default %default]", 
              metavar="character"),
  make_option(c("-f", "--filter"), type="character", default="PASS",
              help="filters to keep separated by , [default %default]", 
              metavar="character"),
  make_option(c("-ro", "--minimum_overlap"), default="0.2",
              help="minimum reciperical overlap [default %default]"),
  make_option(c("-ms", "--minimum_size"), type="integer", default="50",
              help="minimum SV size [default %default]", 
              metavar="integer"),
  make_option(c("-pd", "--pading"), type="integer", default="1000",
              help="pading used to detect TRLX overlaps [default %default]", 
              metavar="integer")
 )
 
#Get command-line arguments & options
#args <- parse_args(OptionParser(usage="%prog [options] covMatrix.bed OUTDIR", option_list=option_list),positional_arguments=TRUE)
args <- parse_args(OptionParser(usage="%prog [options] INDIR OUTDIR", option_list=option_list),positional_arguments=TRUE)
opts <- args$options

if(length(args$args) != 2) 
  {cat("Incorrect number of required positional arguments\n\n")
   stop()}

INDIR <- args$args[1]
OUTDIR <- args$args[2]
M <- opts$mode
PASS_FLAGS <- c(strsplit(opts$filter, ","))
MIN_RIP <- as.numeric(opts$minimum_overlap)
MIN_SIZE <- as.numeric(opts$minimum_size)
PD <- as.numeric(opts$pading)

vcf = readVcf(INDIR)

if(M %in% c('BED','BOTH')){
	##FILTER## <-----TO DO
	begr = sort(breakendRanges(vcf),by=~sourceId)
	begr$score = begr$QUAL
	export(begr, con=OUTDIR)

}
if(M %in% c('BEDPE','BOTH')){
	bpgr = sort(breakpointRanges(vcf),by=~sourceId)
	bedpe = breakpointgr2bedpe(bpgr)
	bedpe$name = gsub('.{1}$', '', bedpe$name)
	bedpe$filter=bpgr[paste(bedpe$name,'o', sep = "")]$FILTER
	bedpe$svlen=bpgr[paste(bedpe$name,'o', sep = "")]$svLen
	bedpe$break1=gsub('\\w+','X',gsub('\\d+:\\d+','',gsub('[XYMT]','0',bpgr[paste(bedpe$name,'h', sep = "")]$ALT)))
	bedpe$break2=gsub('\\w+','X',gsub('\\d+:\\d+','',gsub('[XYMT]','0',bpgr[paste(bedpe$name,'o', sep = "")]$ALT)))
	bedpe$breaksig=gsub('X\\]\\]:X\\]\\]','HEADS',gsub('\\[\\[X:\\[\\[X','TAILS',gsub('X\\[\\[:\\]\\]X','BEFORE',gsub('\\]\\]X:X\\[\\[','AFTER',paste(bedpe$break1,bedpe$break2,sep=":")))))
	bedpe$break1=NULL
	bedpe$break2=NULL
	#DEAL WITH MT,X,Y
	TRLX <- subset(bedpe, chrom1 != chrom2)
	CNV <- subset(subset(bedpe, chrom1 == chrom2), breaksig %in% c("BEFORE","AFTER"))
	OTH <- subset(subset(bedpe, chrom1 == chrom2), breaksig %in% c("HEADS","TAILS"))


	############FILTER CNV################
	CNV_CAND <- subset(CNV, filter %in% PASS_FLAGS & abs(svlen) >= MIN_SIZE)

	############ALIGN OTH################
	#Find heads tials overlaps
	HEADS <- subset(OTH, breaksig == 'HEADS')
	TAILS <- subset(OTH, breaksig == 'TAILS')
	HEADS_RANGES <- GRanges(seqnames = Rle(HEADS$chrom1), ranges = IRanges(HEADS$start1, end = HEADS$end2), strand = Rle(strand('+')), score=HEADS$score, filter=HEADS$filter, signature=HEADS$breaksig, name=HEADS$name)
	TAILS_RANGES <- GRanges(seqnames = Rle(TAILS$chrom1), ranges = IRanges(TAILS$start1, end = TAILS$end2), strand = Rle(strand('+')), score=TAILS$score, filter=TAILS$filter, signature=TAILS$breaksig, name=TAILS$name)
	#OTH_gr <- subsetByOverlaps(OTH_gr,OTH_gr[OTH_gr$filter %in% c("PASS","NO_ASSEMBLY")])
	OVERLAPS <- data.frame(findOverlaps(HEADS_RANGES,TAILS_RANGES))
	OVERLAPS$name_H= mcols(HEADS_RANGES[OVERLAPS$queryHits])$name
	OVERLAPS$name_T= mcols(TAILS_RANGES[OVERLAPS$subjectHits])$name
	OVERLAPS$filter_H = mcols(HEADS_RANGES[OVERLAPS$queryHits])$filter
	OVERLAPS$filter_T = mcols(TAILS_RANGES[OVERLAPS$subjectHits])$filter
	OVERLAPS$chr = as.data.frame(seqnames(HEADS_RANGES[OVERLAPS$queryHits]))$value
	OVERLAPS$size_H = width(ranges(HEADS_RANGES[OVERLAPS$queryHits]))
	OVERLAPS$size_T = width(ranges(TAILS_RANGES[OVERLAPS$subjectHits]))
	OVERLAPS$st_H = start(ranges(HEADS_RANGES[OVERLAPS$queryHits]))
	OVERLAPS$sp_H = end(ranges(HEADS_RANGES[OVERLAPS$queryHits]))
	OVERLAPS$st_T = start(ranges(TAILS_RANGES[OVERLAPS$subjectHits]))
	OVERLAPS$sp_T = end(ranges(TAILS_RANGES[OVERLAPS$subjectHits]))
	OVERLAPS$overlap = width(pintersect(ranges(HEADS_RANGES[OVERLAPS$queryHits]),ranges(TAILS_RANGES[OVERLAPS$subjectHits])))
	OVERLAPS$overlap_H <- width(pintersect(ranges(HEADS_RANGES[OVERLAPS$queryHits]),ranges(TAILS_RANGES[OVERLAPS$subjectHits]))) / width(ranges(HEADS_RANGES[OVERLAPS$queryHits]))
	OVERLAPS$overlap_T <- width(pintersect(ranges(HEADS_RANGES[OVERLAPS$queryHits]),ranges(TAILS_RANGES[OVERLAPS$subjectHits]))) / width(ranges(TAILS_RANGES[OVERLAPS$subjectHits]))
	OVERLAPS$Score = mcols(HEADS_RANGES[OVERLAPS$queryHits])$score + mcols(TAILS_RANGES[OVERLAPS$subjectHits])$score

	############FILTER OTH################
	#Filter based on minimum reciperical overlaps
	OVERLAPS <- subset(OVERLAPS, overlap_H >= MIN_RIP & overlap_T >= MIN_RIP)
	#Filter based on filter field
	OVERLAPS <- subset(OVERLAPS, filter_H %in% PASS_FLAGS | filter_T %in% PASS_FLAGS)
	#Filter based on size field
	OVERLAPS <- subset(OVERLAPS, size_H >= MIN_SIZE | size_T >= MIN_SIZE)
	#Remove any non-overlaps still present
	OVERLAPS <- subset(OVERLAPS, !(sp_T > sp_H & st_T > sp_H))
	OVERLAPS <- subset(OVERLAPS, !(sp_T < st_H & st_T < st_H))

	#Identify duplicate entries and filter based on highest overlap
	HEADS_D <- subset(setDT(OVERLAPS)[, .N, by=queryHits], N > 1)
	#print(subset(OVERLAPS, queryHits == 91))
	SV_CAND <- subset(OVERLAPS, queryHits %in% subset(setDT(OVERLAPS)[, .N, by=queryHits], N == 1)$queryHits)

	for (i in 1:length(HEADS_D$queryHits)){
		#Give only the most well aligned pairs
		tmp <- subset(OVERLAPS, queryHits == HEADS_D$queryHits[i])
		SV_CAND <- rbind(SV_CAND,data.frame(tmp[tmp$overlap_H*tmp$overlap_T == max(tmp$overlap_H*tmp$overlap_T),]))
	}

	TAILS_D <- subset(setDT(SV_CAND)[, .N, by=subjectHits], N > 1)
	#print(TAILS_D)
	SV_CAND <- subset(OVERLAPS, subjectHits %in% subset(setDT(SV_CAND)[, .N, by=subjectHits], N == 1)$subjectHits)
	for (i in 1:length(TAILS_D$subjectHits)){
		#Give only the most well aligned pairs
		tmp <- subset(OVERLAPS, subjectHits == TAILS_D$subjectHits[i])
		#print(tmp)
		SV_CAND <- rbind(SV_CAND,tmp[tmp$overlap_H*tmp$overlap_T == max(tmp$overlap_H*tmp$overlap_T),])
	}

	############DEFINE OTH################
	SV_CAND$type = "INV"
	for (i in 1:(nrow(SV_CAND))){
		if (abs(SV_CAND[i]$st_H-SV_CAND[i]$st_T) >= MIN_SIZE){
			if (SV_CAND[i]$st_H > SV_CAND[i]$st_T){
				if (abs(SV_CAND[i]$sp_H-SV_CAND[i]$sp_T) >= MIN_SIZE){
					if (SV_CAND[i]$sp_H > SV_CAND[i]$sp_T){
						SV_CAND[i]$type = "dupINVdup" #dupINVdup
					}else{
						SV_CAND[i]$type = "dupINVdel" #dupINVdel
					}
				}else{
					SV_CAND[i]$type = "dupINV" #dupINV
				}
			}else{
				if (abs(SV_CAND[i]$sp_H-SV_CAND[i]$sp_T) >= MIN_SIZE){
					if (SV_CAND[i]$sp_H > SV_CAND[i]$sp_T){
						SV_CAND[i]$type = "delINVdup" #delINVdup
					}else{
						SV_CAND[i]$type = "delINVdel" #delINVdel
					}
				}else{
					SV_CAND[i]$type = "delINV" #delINV
				}
			}
		}else{
			if (abs(SV_CAND[i]$sp_H-SV_CAND[i]$sp_T) >= MIN_SIZE){
				if (SV_CAND[i]$sp_H > SV_CAND[i]$sp_T){
					SV_CAND[i]$type = "INVdup" #INVdup
				}else{
					SV_CAND[i]$type = "INVdel" #INVdel
				}
			}
		}
	}

	############ALIGN TRLX################
	############FILTER TRLX################

	#TRLX <- subset(TRLX, filter %in% PASS_FLAGS)
	CHRS <- unique(c(as.character(TRLX$chrom1),as.character(TRLX$chrom2)))
	for (i in 1:(length(CHRS)-1)){
		TRLX_W <- subset(TRLX, chrom1 == CHRS[i] | chrom2 == CHRS[i])
		TRLX <- subset(TRLX, chrom1 != CHRS[i] & chrom2 != CHRS[i])
		if (nrow(TRLX_W) > 0 ){
			CHRS_W <- unique(c(as.character(TRLX_W$chrom1),as.character(TRLX_W$chrom2)))
			for (w in 1:(length(CHRS_W)-1)){
				if (CHRS_W[w] != CHRS[i]){
					TRLX_W_W <- subset(TRLX_W, chrom1 == CHRS_W[w] | chrom2 == CHRS_W[w])
					TRLX_W <- subset(TRLX_W, chrom1 != CHRS_W[w] & chrom2 != CHRS_W[w])
					if (nrow(TRLX_W_W) > 1 ){
						TRLX_W_W_F <- subset (TRLX_W_W, chrom1 == CHRS[i])
						TRLX_W_W_R <- subset(TRLX_W_W, chrom1 != CHRS[i])
						TRLX_W_W_R$ch <- TRLX_W_W_R$chrom1
						TRLX_W_W_R$s <- TRLX_W_W_R$start1
						TRLX_W_W_R$e <- TRLX_W_W_R$end1
						TRLX_W_W_R$chrom1 <- TRLX_W_W_R$chrom2
						TRLX_W_W_R$start1 <- TRLX_W_W_R$start2
						TRLX_W_W_R$end1 <- TRLX_W_W_R$end2
						TRLX_W_W_R$chrom2 <- TRLX_W_W_R$ch
						TRLX_W_W_R$start2 <- TRLX_W_W_R$s
						TRLX_W_W_R$end2 <- TRLX_W_W_R$e
						TRLX_W_W_R$ch <- NULL
						TRLX_W_W_R$s <- NULL
						TRLX_W_W_R$e <- NULL
						TRLX_W_W <- rbind(TRLX_W_W_F, TRLX_W_W_R)
						F_W_RANGES <- GRanges(seqnames = Rle(TRLX_W_W$chrom1), ranges = IRanges(TRLX_W_W$start1-PD, end = TRLX_W_W$end1+PD), strand = Rle(strand('+')), score=TRLX_W_W$score, filter=TRLX_W_W$filter, signature=TRLX_W_W$breaksig, name=TRLX_W_W$name)
						R_W_RANGES <- GRanges(seqnames = Rle(TRLX_W_W$chrom2), ranges = IRanges(TRLX_W_W$start2-PD, end = TRLX_W_W$end2+PD), strand = Rle(strand('+')), score=TRLX_W_W$score, filter=TRLX_W_W$filter, signature=TRLX_W_W$breaksig, name=TRLX_W_W$name)
						F_OVERLAPS <- subset(subset(data.frame(findOverlaps(F_W_RANGES,F_W_RANGES)), queryHits < subjectHits), mcols(F_W_RANGES[queryHits])$filter %in% PASS_FLAGS | mcols(F_W_RANGES[subjectHits])$filter %in% PASS_FLAGS)
						R_OVERLAPS <- subset(subset(data.frame(findOverlaps(R_W_RANGES,R_W_RANGES)), queryHits < subjectHits), mcols(R_W_RANGES[queryHits])$filter %in% PASS_FLAGS | mcols(R_W_RANGES[subjectHits])$filter %in% PASS_FLAGS)
						F_OVERLAPS$sig <- paste(F_OVERLAPS$queryHits, F_OVERLAPS$subjectHits, sep="_")
						R_OVERLAPS$sig <- paste(R_OVERLAPS$queryHits, R_OVERLAPS$subjectHits, sep="_")
						TRLX_OVERLAPS <- intersect(F_OVERLAPS$sig,R_OVERLAPS$sig)
						if (length(TRLX_OVERLAPS) > 0){
							for (trl in TRLX_OVERLAPS){
								idx <- unlist(strsplit(trl,"_"))
								tmp <- TRLX_W_W[as.numeric(idx[1]),]
								tmp$pstart1 <- TRLX_W_W[as.numeric(idx[2]),]$start1
								tmp$pend1 <- TRLX_W_W[as.numeric(idx[2]),]$end1
								tmp$pstart2 <- TRLX_W_W[as.numeric(idx[2]),]$start2
								tmp$pend2 <- TRLX_W_W[as.numeric(idx[2]),]$end2
								tmp$pname <- TRLX_W_W[as.numeric(idx[2]),]$name
								tmp$pfilter <- TRLX_W_W[as.numeric(idx[2]),]$filter
								tmp$pbreaksig <- TRLX_W_W[as.numeric(idx[2]),]$breaksig
								if (exists('TRLX_CAND')){
									TRLX_CAND <- rbind(TRLX_CAND, tmp)
								}else{
									TRLX_CAND <- tmp
								}
							}
						}
						
						if (exists('TRLX_CAND')){
							TRLX_W_W <- subset(subset(TRLX_W_W, filter %in% PASS_FLAGS), !(name %in% TRLX_CAND$name | name %in% TRLX_CAND$pname))
						}else{
							TRLX_W_W <- subset(TRLX_W_W, filter %in% PASS_FLAGS)
						}
						if (length(TRLX_W_W) > 0 & nrow(TRLX_W_W) > 0){
							#DEAL WITH ORPAHNS
							if (exists('TRLX_ORPHAN')){
								TRLX_ORPHAN <- rbind(TRLX_ORPHAN, TRLX_W_W)
							}else{
								TRLX_ORPHAN <- TRLX_W_W
							}
						}
					}else if (nrow(TRLX_W_W) == 1){
						TRLX_W_W <- subset(TRLX_W_W, filter %in% PASS_FLAGS)
						if(length(TRLX_W_W) > 0 & nrow(TRLX_W_W) > 0){
							#DEAL WITH ORPHANS
							if (exists('TRLX_ORPHAN')){
								TRLX_ORPHAN <- rbind(TRLX_ORPHAN, TRLX_W_W)
							}else{
								TRLX_ORPHAN <- TRLX_W_W
							}
						}
					}
				}
			}
		}
	}
	DUP_CAND <-subset(CNV_CAND, breaksig == "BEFORE")
	DEL_CAND <-subset(CNV_CAND, breaksig == "AFTER")
	DUP_CAND <- data.frame(DUP_CAND$chrom1, DUP_CAND$start1, DUP_CAND$end2, "-", "-", "-", DUP_CAND$name, DUP_CAND$score, DUP_CAND$filter, "DUP",stringsAsFactors=FALSE)
	colnames(DUP_CAND) <- c("chrom1","start1","end1","chrom2","start2","stop2","name","score","filter","type")
	DEL_CAND <- data.frame(DEL_CAND$chrom1, DEL_CAND$start1, DEL_CAND$end2, "-", "-", "-", DEL_CAND$name, DEL_CAND$score, DEL_CAND$filter, "DEL",stringsAsFactors=FALSE)
	colnames(DEL_CAND) <- c("chrom1","start1","end1","chrom2","start2","stop2","name","score","filter","type")
	TRLX_CAND$fullname <- paste(TRLX_CAND$name, TRLX_CAND$pname, sep="|")
	TRLX_CAND$fullfilter <- paste(TRLX_CAND$filter, TRLX_CAND$pfilter, sep="|")
	TRLX_CAND <- data.frame(TRLX_CAND$chrom1, TRLX_CAND$start1, TRLX_CAND$pend1, TRLX_CAND$chrom2, TRLX_CAND$start2, TRLX_CAND$pend2, TRLX_CAND$fullname, TRLX_CAND$score, TRLX_CAND$fullfilter, "TRLX")
	colnames(TRLX_CAND) <- c("chrom1","start1","end1","chrom2","start2","stop2","name","score","filter","type")
	TRLX_ORPHAN <- data.frame(TRLX_ORPHAN$chrom1, TRLX_ORPHAN$start1, TRLX_ORPHAN$end1, TRLX_ORPHAN$chrom2, TRLX_ORPHAN$start2, TRLX_ORPHAN$end2, TRLX_ORPHAN$name, TRLX_ORPHAN$score, TRLX_ORPHAN$filter, "ORPHAN_TRLX")
	colnames(TRLX_ORPHAN) <- c("chrom1","start1","end1","chrom2","start2","stop2","name","score","filter","type")
	SV_CAND$full_name <- paste(SV_CAND$name_H, SV_CAND$name_T, sep="|")
	SV_CAND$full_filter <- paste(SV_CAND$filter_H, SV_CAND$filter_T, sep="|")
	SV_CAND <- data.frame(SV_CAND$chr, SV_CAND$st_H, SV_CAND$sp_T, SV_CAND$chr, SV_CAND$st_T, SV_CAND$sp_H, SV_CAND$full_name, SV_CAND$Score, SV_CAND$full_filter, SV_CAND$type)
	colnames(SV_CAND) <- c("chrom1","start1","end1","chrom2","start2","stop2","name","score","filter","type")
	bedpe <- rbind(DUP_CAND, DEL_CAND, TRLX_CAND, TRLX_ORPHAN, SV_CAND)
	###OUTPUT###
	write.table(bedpe, file=OUTDIR, sep="\t", quote=FALSE, col.names=FALSE)
}
