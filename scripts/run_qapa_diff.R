library(getopt)

options<-matrix(c(
  "help", "h", "0", "logical","help",
  "input_pau", "i", "1", "character","input pau result",
  "prefix_txt", "p", "1", "character","prefix table text",
  "cond_compare", "c", "1", "character","condition compare",
  "out_path", "o", "1", "character","out path"
), ncol=5, byrow=T)

args=getopt(options)

if(is.null(args$input_pau) | is.null(args$prefix_txt) | is.null(args$cond_compare) | is.null(args$out_path)){
  cat(getopt(options, usage=T), "\n")
  q() 
}

input_pau = args$input_pau
prefix_txt = args$prefix_txt
out_path = paste0(args$out_path,"/")
cond_compare <- args$cond_compare

a <- read.table(input_pau,header = T, sep = "\t",check.names = F)
a_topcol <- a[,1:12]
a_pau <- a[,grepl("PAU",colnames(a))]
a_tpm <- a[,grepl("TPM",colnames(a))]

prefix_table = read.table(prefix_txt,header = F, sep = '\t',check.names = F)
cond_split = unlist(strsplit(cond_compare, "_vs_"))
cond_list = c(cond_split[2],cond_split[1])
cond1 = cond_list[1]
cond2 = cond_list[2]
cond1_prefix <- prefix_table[,1][prefix_table[,2] == cond_list[1]]
cond2_prefix <- prefix_table[,1][prefix_table[,2] == cond_list[2]]
gene_df <- function(cond1, cond2, cond1_prefix, cond2_prefix, a_tpm, a_pau){
  cond1_pau <- paste0(cond1_prefix,".PAU")
  cond2_pau <- paste0(cond2_prefix,".PAU")
  cond1_tpm <- paste0(cond1_prefix,".TPM")
  cond2_tpm <- paste0(cond2_prefix,".TPM")
  a_cond1 <- data.frame(a_pau[,cond1_pau])
  a_cond2 <- data.frame(a_pau[,cond2_pau])
  a_cond1_tpm <- data.frame(a_tpm[,cond1_tpm])
  a_cond2_tpm <- data.frame(a_tpm[,cond2_tpm])
  a_merge_tpm <- cbind(a_cond1_tpm, a_cond2_tpm)
  
  a_merge <- cbind(a_cond1, a_cond2)
  a_merge_flt1 <- a_merge
  a_merge_flt2 <- a_merge_flt1
  n1 = ncol(a_cond1)
  n2 = ncol(a_cond2)
  cond1_mean_name = paste0(cond1,"_meanPAU")
  cond2_mean_name = paste0(cond2,"_meanPAU")
  a_merge[,cond1_mean_name] <- rowMeans(a_merge[,1:n1])
  a_merge[,cond2_mean_name] <- rowMeans(a_merge[,(n1+1):(n1+n2)])
  a_merge$DPAU <- a_merge[,cond2_mean_name] - a_merge[,cond1_mean_name]
  a_topcol_flt <- a_topcol[as.character(rownames(a_merge)),]
  pau_res <- cbind(a_topcol_flt,a_merge)
  return(pau_res)
}
pau_test <- gene_df(cond1, cond2, cond1_prefix, cond2_prefix, a_tpm, a_pau)
out_file = paste0(cond2,"_vs_",cond1, "_pau_result.txt")
write.table(pau_test, paste0(out_path, out_file), quote = F, sep = "\t", row.names = F)
