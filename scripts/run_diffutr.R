library(diffUTR)
library(getopt)
library(GenomicRanges)

options<-matrix(c(
  "help", "h", "0", "logical","help",
  "bam_txt", "b", "1", "character","bam text",
  "out_path", "o", "1", "character","out path",
  "anno_path","a","1","character","annotation path",
  "cov_cutoff","c","2","integer","coverage cutoff",
  "reads_len","r","2","integer","reads length"
), ncol=5, byrow=T)
args=getopt(options)
if(is.null(args$bam_txt) | is.null(args$out_path) | is.null(args$anno_path)){
  cat(getopt(options, usage=T), "\n")
  q() 
}

anno_path = args$anno_path
cov_cutoff = args$cov_cutoff
reads_len = args$reads_len
bam_txt = args$bam_txt

diffutr_path = paste0(args$out_path,"/")
out_path = paste0(diffutr_path,"out/")
data_path = paste0(diffutr_path,"process/")

bin_bed = read.table(anno_path,header = T,sep = "\t")
bins = GRanges(seqnames = Rle(bin_bed$seqnames),
				ranges =  IRanges(bin_bed$start, end=bin_bed$end),
				strand = Rle(bin_bed$strand),
				gene_name = bin_bed$gene_name,
				type = bin_bed$type,
				gene = bin_bed$gene_name,
				bin_id = bin_bed$bin_id,
				geneAmbiguous = bin_bed$ambiguous)
if(is.null(names(bins))) names(bins) <- 1:length(bins$gene_name)

bam_table = read.table(bam_txt,header = F, sep = '\t')
class = unique(bam_table[,2])
con1 = class[1]
con2 = class[2]
con1_bam <- bam_table[,1][bam_table[,2]==con1]
con2_bam <- bam_table[,1][bam_table[,2]==con2]
com_bam <- c(con1_bam, con2_bam)
rse <- countFeatures(bamfiles=com_bam,
                    bins=bins,
					strandSpecific=1,
					nthreads=4,
					isPairedEnd=T,
					readLength=reads_len)

stage <- c(rep("control",length(con1_bam)),rep("treatment",length(con2_bam)))
rse@colData$condition <- stage
rse <- diffSpliceWrapper(se=rse, design = ~condition)
saveRDS(rse,paste0(data_path,con2,"_vs_",con1,".res.rds"))

utr3 <- geneLevelStats(se=rse, includeTypes="3UTR", returnSE=F, minDensityRatio = cov_cutoff, minWidth = 20)
write.table(utr3,paste0(out_path,con2,"_vs_",con1,".txt"),sep = "\t",quote = F)

