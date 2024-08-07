#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(ggplot2)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-s", "--data_set"), type="character", default="example.dataset", 
              help="dataset name"),
  make_option(c("-d", "--distances"), type="character", default='jc.distance.precalc.csv', 
              help="genome distance matrix file name, should be an csv"),
  make_option(c("--ngenerations"), type="numeric", default=10, 
              help="number of generations"),
  make_option(c("--n.batches"), type="numeric", default=1, 
              help="number of overall batches per generation"),
  make_option(c("--metadata"), type="character", default='metadata.csv', 
              help="metadata file, should be a csv with the date column called Collection.date in the %d/%m/%Y format"),
  # make_option(c("--distopt"), type="character", default="max", 
  #             help="Genetic distance optimized towards its highest/average point"),
  make_option(c("--outdir"), type="character", default=getwd(), 
              help="output directory")
              )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

metadata = read.table(opt$metadata, sep = ',', stringsAsFactors = FALSE, header = TRUE)

if(grepl("\\.rds$", opt$distances) | grepl("\\.RDS$", opt$distances)){
  dist.m = readRDS(opt$distances)
}else{
  dist.m = scan(opt$distances, sep = ',', what=character())
  dist.m = matrix(dist.m, nrow=dim(metadata)[1]+1)
  genome.names <- dist.m[1,2:dim(dist.m)[1]]
  dist.m = apply((dist.m[2:dim(dist.m)[1],2:dim(dist.m)[1]]), 1, as.numeric)
  colnames(dist.m) <- rownames(dist.m) <- genome.names
  }

mean.d = mean(dist.m)
median.d = median(dist.m)
sd.d = sd(dist.m)
sd.up = mean.d+sd.d
if(mean.d-sd.d>0){
  sd.down = mean.d-sd.d
}else{
  sd.down = 0
}

fitnesses = NULL
indeces = NULL

df.to.plot = data.frame(generation=numeric(), value=numeric(), group=character(), stringsAsFactors = FALSE)

for(g in 1:opt$ngenerations-1){
  gen.output = NULL
  for (i in 1:opt$n.batches){
    in.file = paste0(opt$outdir, '/GA.', opt$data_set, '.', g, '.', i)

    tmp = as.matrix(read.csv(paste(in.file, 'indeces.fitness.csv', sep = '.'), sep = ','))
    gen.output = cbind(gen.output, tmp)
    }
  df.to.plot = df.to.plot %>% add_row(generation=g, value=mean(gen.output[1,]), group='mean fitness')
  df.to.plot = df.to.plot %>% add_row(generation=g, value=max(gen.output[1,]), group='best fitness')
  df.to.plot = df.to.plot %>% add_row(generation=g, value=mean(gen.output[2,]), group='mean genetic diversity')
  df.to.plot = df.to.plot %>% add_row(generation=g, value=mean(gen.output[3,]), group='mean temporal spread')
  }
out.file = paste0(opt$outdir, '/', opt$data_set, '.per.gen.stats')
write.csv(df.to.plot, paste(out.file, 'csv', sep = '.'), quote = FALSE, row.names = FALSE)

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
  geom_hline(yintercept=sd.up, linetype="dashed") +
  geom_hline(yintercept=sd.down, linetype="dashed") +
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

ggsave(filename=paste(out.file, 'png', sep = '.'), plot=do.call("grid.arrange", c(list(p.fit, p.gen.div, p.tem.spr), ncol=3)), width=30, height=8, units = "cm", dpi = 600)
