#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(ggplot2)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--shiny"), type="logical", default="FALSE", 
              help="R shiny or NextFlow"),
  make_option(c("-s", "--data.set"), type="character", default="Norway", 
              help="dataset name"),
  make_option(c("-d", "--distance"), type="character", default='jc.distance.precalc.csv', 
              help="output genome distance matrix file name, should be an csv"),
  make_option(c("--generations"), type="numeric", default=10, 
              help="number of generations"),
  make_option(c("--n.subsamples"), type="numeric", default=3, 
              help="number of top ranking subsampling solutions to consider as output"),
  make_option(c("--tot.batches"), type="numeric", default=1, 
              help="number of overall batches per generation"),
  make_option(c("--metadata"), type="character", default='metadata.csv', 
              help="metadata file, should be a csv with the date column called Collection.date in the %d/%m/%Y format"),
  make_option(c("--dist.opt"), type="character", default="highest", 
              help="Genetic distance optimized towards its highest/average point")
              )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(!opt$shiny){
  source('bin/read.distance.matrix.R')
  config <- read.table(paste(opt$data.set, 'config', sep='.'), header=FALSE, sep="=", row.names=1, strip.white=TRUE, na.strings="NA", stringsAsFactors=FALSE)
  config <- structure(as.character(config[,1]), names = rownames(config))
  config.default = setNames(c(5,1,3),
                            c("params.tot.batches",
                              "params.ngenerations",
                              "params.n.subsamples"))
  config = c(config, config.default[!config.default %in% config])
  batches = config['params.tot.batches']
  gen = as.numeric(config['params.ngenerations'])
  n.subsamples = config['params.n.subsamples']
}else{
  source('../bin/read.distance.matrix.R')
  batches = opt$tot.batches
  gen = opt$generations
  n.subsamples= opt$n.subsamples
}

metadata = read.table(opt$metadata, sep = ',', stringsAsFactors = FALSE, header = TRUE)

if(grepl("\\.rds$", opt$distance) | grepl("\\.RDS$", opt$distance)){
  dist.m = readRDS(opt$distance)
}else{
  dist.m = scan(opt$distance, sep = ',', what=character())
  dist.m = matrix(dist.m, nrow=dim(metadata)[1]+1)
  genome.names <- dist.m[1,2:dim(dist.m)[1]]
  dist.m = apply((dist.m[2:dim(dist.m)[1],2:dim(dist.m)[1]]), 1, as.numeric)
  colnames(dist.m) <- rownames(dist.m) <- genome.names
  }

mean.d = mean(dist.m)
median.d = median(dist.m)
sd.d = sd(dist.m)

fitnesses = NULL
indeces = NULL

for (i in 1:batches){
  if (!opt$shiny){
    in.file = paste0('output/', opt$data.set, '/GA.', opt$data.set, '.', gen-1, '.', i)
    }else{
      in.file = paste0('output/GA.', gen-1, '.', i)
      }
  
  tmp_fit = as.matrix(read.csv(paste(in.file, 'indeces.fitness.csv', sep = '.'), sep = ','))[2,]
  fitnesses = c(fitnesses, tmp_fit)
    
  tmp_ind = as.matrix(read.csv(paste(in.file, 'indeces.subsamples.csv', sep = '.'), sep = ',', stringsAsFactors = FALSE, header = FALSE))
  indeces = rbind(indeces, tmp_ind)
  }
  
names(fitnesses) = as.character(1:length(fitnesses))
selected.individuals = sort(fitnesses, decreasing = TRUE)[1:n.subsamples]

for (i in 1:n.subsamples){
  selected.genomes = colnames(dist.m)[indeces[as.numeric(names(selected.individuals[i])),]]
  if (!opt$shiny){
    out.file = paste0('output/', opt$data.set, '/subsample.GA.', i, '.csv')
  }else{
    out.file = paste0('output/subsample.GA.', i, '.csv')
  }
  writeLines( selected.genomes, out.file )
}

df.to.plot = data.frame(generation=numeric(), value=numeric(), group=character(), stringsAsFactors = FALSE)

for(g in 1:gen-1){
  gen.output = NULL
  for (i in 1:batches){
    if (!opt$shiny){
      in.file = paste0('output/', opt$data.set, '/GA.', opt$data.set, '.', g, '.', i)
    }else{
      in.file = paste0('output/GA.', g, '.', i)
    }
    
    tmp = as.matrix(read.csv(paste(in.file, 'indeces.fitness.csv', sep = '.'), sep = ','))
    gen.output = cbind(gen.output, tmp)
    }
  df.to.plot = df.to.plot %>% add_row(generation=g, value=mean(gen.output[1,]), group='mean fitness')
  df.to.plot = df.to.plot %>% add_row(generation=g, value=max(gen.output[1,]), group='best fitness')
  df.to.plot = df.to.plot %>% add_row(generation=g, value=mean(gen.output[2,]), group='mean genetic diversity')
  df.to.plot = df.to.plot %>% add_row(generation=g, value=mean(gen.output[3,]), group='mean temporal spread')
  }
write.csv(df.to.plot, 'output/per.gen.stats.csv', quote = FALSE, row.names = FALSE)

p.fit <- ggplot(df.to.plot[df.to.plot$group %in% c('mean fitness', 'best fitness'),], aes(x=generation, y=value, colour=group)) + 
  geom_line() +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values=c('brown1', 'darkblue')) +
  ylab('fitness') +
  theme(legend.position = 'bottom', legend.title = element_blank())

p.gen.div <- ggplot(df.to.plot[df.to.plot$group %in% c('mean genetic diversity', 'best genetic diversity'),], aes(x=generation, y=value, colour=group)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept=mean.d) +
  geom_hline(yintercept=mean.d+sd.d, linetype="dashed") +
  geom_hline(yintercept=min(mean.d-sd.d, 0), linetype="dashed") +
  geom_hline(yintercept=median.d, color='cornflowerblue') +
  geom_text(aes(x=0, y = 0.9*median.d), label="Median", color='cornflowerblue', size=3, nudge_x = 0.5) +
  geom_text(aes(x=0, y = 0.9*mean.d), label="Mean", size=3, nudge_x = 0.5, color = 'black') +
  scale_color_manual(values=c('darkblue')) +
  theme_minimal() +
  ylab('genomic diversity') +
  theme(legend.position = 'bottom', legend.title = element_blank())

p.tem.spr <- ggplot(df.to.plot[df.to.plot$group %in% c('mean temporal spread', 'best temporal spread'),], aes(x=generation, y=value, colour=group)) +
  geom_line() +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c('darkblue')) +
  ylab('temporal distribution') +
  theme(legend.position = 'bottom', legend.title = element_blank())

ggsave(filename='output/per.gen.stats.png', plot=do.call("grid.arrange", c(list(p.fit, p.gen.div, p.tem.spr), ncol=3)), width=30, height=8, units = "cm", dpi = 600)
