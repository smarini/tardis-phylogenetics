library(optparse)
library(DECIPHER)

option_list = list(
  make_option(c("-i", "--alignment"), type="character", default="aln.fa", 
              help="aligned genome fasta file name", metavar="character"),
  make_option(c("-d", "--distance"), type="character", default='jc.distance.precalc.rds', 
              help="output genome distance mastrix file name, should be an csv or a rds", metavar="character"),
  make_option(c("-c", "--n.cores"), type="numeric", default=1, 
              help="number of cores for parallel comp.", metavar="numeric")
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

### TAKE GENOMES
aligned.genomes = readDNAStringSet(opt$alignment, format='fasta')
aligned.genomes = replaceAmbiguities(aligned.genomes) # replace ambiguities

jc <- DistanceMatrix(aligned.genomes,
                     includeTerminalGaps = FALSE,
                     penalizeGapLetterMatches = TRUE,
                     penalizeGapGapMatches = FALSE,
                     correction = "Jukes-Cantor",
                     processors = opt$n.cores,
                     verbose = TRUE)

message('distances calculated')

if(grepl("\\.rds$", opt$distance) | grepl("\\.RDS$", opt$distance)){
  saveRDS(jc, file = opt$distance)
}else{
  write.csv(jc, file = opt$distance, quote = FALSE)
}
  
