### R code from vignette source 'exomeCopy.Rnw'

###################################################
### code chunk number 1: exomeCopy.Rnw:18-19
###################################################
options(width = 70)


###################################################
### code chunk number 2: exomeCopy.Rnw:55-77 (eval = FALSE)
###################################################
## library(exomeCopy)
## target.file <- "targets.bed"
## bam.files <- c("/path/to/file1.bam", "/path/to/file2.bam", "/path/to/file3.bam")
## sample.names <- c("sample1","sample2","sample3")
## reference.file <- "/path/to/reference_genome.fa"
## target.df <- read.delim(target.file, header = FALSE)
## target <- GRanges(seqname = target.df[,1], IRanges(start = target.df[,2] + 1, end = target.df[,3]))
## counts <- target
## for (i in 1:length(bam.files)) {
##   mcols(counts)[[sample.names[i]]] <- countBamInGRanges(bam.files[i], target)
## }
## counts$GC <- getGCcontent(target, reference.file)
## counts$GC.sq <- counts$GC^2
## counts$bg <- generateBackground(sample.names, counts, median)
## counts$log.bg <- log(counts$bg + .1)
## counts$width <- width(counts)
## fit.list <- lapply(sample.names, function(sample.name) {
##   lapply(seqlevels(target), function(seq.name) { 
##     exomeCopy(counts[seqnames(counts) == seq.name], sample.name, X.names = c("log.bg", "GC", "GC.sq","width"), S = 0:4, d = 2)
##   })
## })
## compiled.segments <- compileCopyCountSegments(fit.list)


###################################################
### code chunk number 3: exomeCopy.Rnw:90-93
###################################################
library(exomeCopy)
gr <- GRanges(seqname="seq1",IRanges(start=1,end=345))
subdivideGRanges(gr)


###################################################
### code chunk number 4: exomeCopy.Rnw:99-106
###################################################
  plot(0,0,xlim=c(0,500),ylim=c(0,25),type="n",yaxt="n",ylab="",xlab="width of input GRanges object",main="Effect of subdivideGRanges")
  abline(v=1:5*100,col="grey")
  for (i in 1:24) {
    gr <- GRanges(seqname="chr1",IRanges(start=1,width=(i*20)))
    sbd.gr <- subdivideGRanges(gr)
    arrows(start(sbd.gr),rep(i,length(sbd.gr)),end(sbd.gr),rep(i,length(sbd.gr)),length=.04,angle=90,code=3)
  }


###################################################
### code chunk number 5: exomeCopy.Rnw:112-119
###################################################
target.file <- system.file("extdata","targets.bed",package="exomeCopy")
target.df <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end")) 
target <- GRanges(seqname=target.df$seqname,IRanges(start=target.df$start+1,end=target.df$end))
target
target <- reduce(target, min.gapwidth=0)
target.sub <- subdivideGRanges(target)
target.sub


###################################################
### code chunk number 6: exomeCopy.Rnw:128-137
###################################################
bam.file <- system.file("extdata","mapping.bam",package="exomeCopy")
scanBamHeader(bam.file)[[1]]$targets
seqlevels(target.sub)
toy.counts <- target.sub
sample.df <- data.frame(samples="sample1",bam.files=bam.file,stringsAsFactors=FALSE)
for (i in 1:nrow(sample.df)) {
  mcols(toy.counts)[[sample.df$samples[i]]] <- countBamInGRanges(sample.df$bam.files[i],target.sub)
}
toy.counts


###################################################
### code chunk number 7: exomeCopy.Rnw:143-146
###################################################
reference.file <- system.file("extdata","reference.fa",package="exomeCopy")
toy.counts$GC <- getGCcontent(target.sub, reference.file)
toy.counts


###################################################
### code chunk number 8: exomeCopy.Rnw:155-158
###################################################
data(exomecounts)
length(exomecounts)
head(exomecounts)


###################################################
### code chunk number 9: exomeCopy.Rnw:164-165
###################################################
plot(start(exomecounts),exomecounts$HG00551,xlim=c(0.8e6,1.8e6),xlab="genomic position",ylab="counts",main="HG00551 read counts in exonic ranges")


###################################################
### code chunk number 10: exomeCopy.Rnw:173-176
###################################################
chr1a <- exomecounts
seqlevels(chr1a) <- "chr1a"
example.counts <- c(exomecounts, chr1a)


###################################################
### code chunk number 11: exomeCopy.Rnw:183-187
###################################################
exome.samples <- grep("HG.+",colnames(mcols(example.counts)),value=TRUE)
example.counts$bg <- generateBackground(exome.samples, example.counts, median)
example.counts$log.bg <- log(example.counts$bg + .1)
example.counts$bg.var <- generateBackground(exome.samples, example.counts, var)


###################################################
### code chunk number 12: exomeCopy.Rnw:192-194
###################################################
summary(example.counts$bg)
example.counts <- example.counts[example.counts$bg > 0,]


###################################################
### code chunk number 13: exomeCopy.Rnw:199-201
###################################################
example.counts$GC.sq <- example.counts$GC^2
example.counts$width <- width(example.counts)


###################################################
### code chunk number 14: exomeCopy.Rnw:245-257
###################################################
simulateCNV <- function(x,indices,multiply,prob) {
  x[indices] <- x[indices] + multiply * rbinom(length(indices),prob=prob,size=x[indices])
  return(x)
}
set.seed(2)
cnv.probs <- rep(c(.99,.5,.5,.95),each=2)
cnv.mult <- rep(c(-1,1),each=4)
cnv.starts <- rep(c(1,301,601,901),each=2)
for (i in 1:8) {
  mcols(example.counts)[[exome.samples[i]]] <- simulateCNV(mcols(example.counts)[[exome.samples[i]]],cnv.starts[i]:(cnv.starts[i] + 99),multiply=cnv.mult[i],prob=cnv.probs[i])
  mcols(example.counts)[[exome.samples[i]]] <- simulateCNV(mcols(example.counts)[[exome.samples[i]]],1000 + cnv.starts[i]:(cnv.starts[i] + 99),multiply=cnv.mult[i],prob=cnv.probs[i])
}


###################################################
### code chunk number 15: exomeCopy.Rnw:266-268
###################################################
  fit <- exomeCopy(example.counts[seqnames(example.counts) == "chr1"],sample.name=exome.samples[3],X.names=c("log.bg","GC","GC.sq","width"),S=0:6,d=2)
  show(fit)


###################################################
### code chunk number 16: exomeCopy.Rnw:273-274
###################################################
  copyCountSegments(fit)


###################################################
### code chunk number 17: exomeCopy.Rnw:280-282
###################################################
  cnv.cols <- c("red","orange","black","deepskyblue","blue","blue2","blue4")
  plot(fit,col=cnv.cols)


###################################################
### code chunk number 18: exomeCopy.Rnw:290-293
###################################################
runExomeCopy <- function(sample.name,seqs) {
  lapply(seqs,function(seq.name) exomeCopy(example.counts[seqnames(example.counts) == seq.name],sample.name,X.names=c("log.bg","GC","GC.sq","width"),S=0:4,d=2))
}


###################################################
### code chunk number 19: exomeCopy.Rnw:298-308
###################################################
seqs <- c("chr1","chr1a")
names(seqs) <- seqs
samples <- exome.samples[1:8]
names(samples) <- samples
t0 <- as.numeric(proc.time()[3])
fit.list <- lapply(samples,runExomeCopy,seqs)
t1 <- as.numeric(proc.time()[3])
time.elapsed <- as.numeric(t1 - t0)
paste(round(time.elapsed),"seconds for",length(samples),"samples,",round(sum(width(example.counts))/1e3),"kb of target")
paste("~",round(time.elapsed / 60 / (8 * sum(width(example.counts))) * 32e6,1)," minutes for 1 sample, 32 Mb of target",sep="")


###################################################
### code chunk number 20: exomeCopy.Rnw:315-317
###################################################
compiled.segments <- compileCopyCountSegments(fit.list)
CNV.segments <- compiled.segments[compiled.segments$copy.count != 2]


###################################################
### code chunk number 21: exomeCopy.Rnw:325-328
###################################################
CNV.segments[1:6]
table(CNV.segments$nranges)
CNV.segments <- CNV.segments[CNV.segments$nranges > 5]


###################################################
### code chunk number 22: exomeCopy.Rnw:334-336
###################################################
CNV.overlaps.matrix <- as.matrix(findOverlaps(CNV.segments,drop.self=TRUE))
head(CNV.overlaps.matrix)


###################################################
### code chunk number 23: exomeCopy.Rnw:342-346
###################################################
par(mfrow=c(2,1),mar=c(4,3,2,1))
cnv.cols <- c("red","orange","black","deepskyblue","blue")
plotCompiledCNV(CNV.segments=CNV.segments, seq.name="chr1", col=cnv.cols)
plotCompiledCNV(CNV.segments=CNV.segments, seq.name="chr1a", col=cnv.cols)


###################################################
### code chunk number 24: exomeCopy.Rnw:354-356
###################################################
fit.list[[1]][[1]]@init.par$beta.hat
fit.list[[1]][[1]]@final.par$beta


###################################################
### code chunk number 25: exomeCopy.Rnw:361-362
###################################################
sessionInfo()


