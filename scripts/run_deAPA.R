library(deAPA)
library(getopt)

options<-matrix(c(
  "help", "h", "0", "logical","help",
  "input_file", "i", "1", "character","input file by predictAPA",
  "output_file", "o", "1", "character","out path file",
  "cov_cutoff", "c", "1", "numeric","coverage cutoff"
), ncol=5, byrow=T)

args=getopt(options)

if(is.null(args$input_file) | is.null(args$output_file) | is.null(args$cov_cutoff)){
  cat(getopt(options, usage=T), "\n")
  q() 
}

input_file <- args$input_file
output_file <- args$output_file
cov_cut_off <- args$cov_cutoff

deAPA(input_file=input_file, 
		output_file=output_file, 
		group1=1, 
		group2=2, 
		least_qualified_num_in_group1=1, 
		least_qualified_num_in_group2=1, 
		coverage_cutoff=cov_cut_off)
