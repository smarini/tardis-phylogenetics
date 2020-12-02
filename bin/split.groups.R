#!/usr/bin/env Rscript

# libraries
library(optparse)
library(Biostrings)

option_list = list(
  make_option("--base.dir", type="character", default="data/example_group", 
              help="data folder for the project; will be used to generate group sub data folders"),
  make_option("--genome.file", type="character", default="data/example_group/aln.group.fa", 
              help="file with genomes, format fasta"),
  make_option("--metadata", type="character", default='data/example_group/metadata.group.csv', 
              help="metadata file, should be a csv with the date column called Collection.date in the %d/%m/%Y format; and column Group"),
  make_option("--group.parameters", type="character", default='data/example_group/parameters.group.csv', 
              help="parameter file; one line per group; unspecified options will go to default")
  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

genomes = readDNAStringSet(opt$genome.file)
metadata = read.table(opt$metadata, sep = ',', stringsAsFactors = FALSE, header = TRUE)
metadata = metadata[sapply(names(genomes), function(x,df){which(metadata$Accession.ID == x)}, df=df),] # make sure id order is respected

parameters = read.csv(opt$group.parameters, stringsAsFactors=FALSE)

# generate project folder
for(p in unique(metadata$Group)){
  data.dir = paste(opt$base.dir, 'data', p, sep = '/')
  system(paste('mkdir -p', data.dir)) # here we store stuff
  metadata.tmp = metadata[metadata$Group %in% p,]
  writeXStringSet(genomes[names(genomes) %in% metadata.tmp$Accession.ID], filepath = paste(data.dir, 'aln.fa', sep = '/') )
  write.csv(metadata.tmp, file = paste(data.dir, 'metadata.csv', sep = '/'), quote = FALSE, row.names = FALSE )
  
  if(p %in% parameters$group){
    options.tardis = parameters[which(parameters$group == p),-c(which(colnames(parameters) == 'group'), which(is.na(parameters[which(parameters$group == p),]))), drop = FALSE]
    parameter.file = c(
      paste0("params.data_set = \"", parameters$group[p], "\""),
      paste(colnames(options.tardis), '=', options.tardis)
    )
    parameter.file.name = paste0(data.dir, '/', p, '.config')
    writeLines(parameter.file, con=parameter.file.name)
    write(parameter.file.name, "")
  }
}
