#!/usr/bin/env Rscript

# libraries
library(optparse)
library(Biostrings)

option_list = list(
  make_option("--base.dir", type="character", default="data/example_group", 
              help="data folder for the project; will be used to generate group sub data folders"),
  make_option("--genome.file", type="character", default="data/example_group.dataset", 
              help="file with genomes, format fasta"),
  make_option("--metadata", type="character", default='data/example_group/metadata.csv', 
              help="metadata file, should be a csv with the date column called Collection.date in the %d/%m/%Y format; and column Group"),
  make_option("--group.parameters", type="character", default='data/example_group/parametes.group.csv', 
              help="parameter file; one line per group; unspecified options will go to default")
  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

genomes = readDNAStringSet(opt$genome.file)
metadata = read.table(opt$metadata, sep = ',', stringsAsFactors = FALSE, header = TRUE)
metadata = metadata[sapply(names(genomes), function(x,df){which(metadata$Accession.ID == x)}, df=df),] # make sure id order is respected

parameters = read.csv(opt$group.parameters, stringsAsFactors=FALSE)

# generate project folder
for(p in 1:dim(parameters)[1]){
  data.dir = paste(opt$base.dir, 'data', parameters$group[p], sep = '/')
  system(paste('mkdir -p', data.dir)) # here we store stuff
  writeXStringSet(genomes[names(genomes) %in% metadata$Accession.ID], filepath = paste(data.dir, 'aln.fa', sep = '/') )
  write.csv(metadata[metadata$Group %in% parameters$group[p],], file = paste(data.dir,  'metadata.csv', sep = '/') )
  
  options.tardis = parameters[p,-c(which(colnames(parameters) == 'group'), which(is.na(parameters[p,]))), drop = FALSE]
  parameter.file = c(
    paste0("params.data_set = \"", parameters$group[p], "\""),
    paste(colnames(options.tardis), '=', options.tardis)
  )
  writeLines(parameter.file, con=paste0(data.dir, '/', parameters$group[p], '.config'))
}

# genomes of groups not listed in parameters will not be subsampled
# and will be fully included in the solution
not.to.subsample = unique(metadata$Group)[!unique(metadata$Group) %in% parameters$group]

for(p in not.to.subsample){
  data.dir = paste(opt$base.dir, 'data', p, sep = '/')
  system(paste('mkdir -p', data.dir)) # here we store stuff
  writeXStringSet(genomes[names(genomes) %in% metadata$Accession.ID], filepath = paste(data.dir, 'aln.fa', sep = '/') )
  write.csv(metadata[metadata$Group %in% parameters$group[p],], file = paste(data.dir,  'metadata.csv', sep = '/') )
}
