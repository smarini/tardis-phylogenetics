#!/usr/bin/env Rscript

# libraries
library(doRNG)
library(optparse)

read.distance.matrix = function(dist.file){
  
  if(grepl("\\.rds$", dist.file) | grepl("\\.RDS$", dist.file)){
    dist.m = readRDS(dist.file)
  }else{
    dist.m = scan(opt$distance, sep = ',', what=character())
    dist.m = matrix(dist.m, nrow=dim(metadata)[1]+1)
    genome.names <- dist.m[1,2:dim(dist.m)[1]]
    dist.m = apply((dist.m[2:dim(dist.m)[1],2:dim(dist.m)[1]]), 1, as.numeric)
    colnames(dist.m) <- rownames(dist.m) <- genome.names
  }
  
  return(dist.m)
}

option_list = list(
  make_option(c("--data.set"), type="character", default="example.dataset", 
              help="dataset name"),
  make_option(c("--distance"), type="character", default='data/example/jc.distance.precalc.csv', 
              help="genome distance matrix file, should be a csv"),
  make_option(c("--metadata"), type="character", default='data/example/metadata.csv', 
              help="metadata file, should be a csv with the date column called Collection.date in the %d/%m/%Y format"),
  make_option(c("--n.cores"), type="numeric", default=1, 
              help="number of cores for parallel comp."),
  make_option(c("--frac.new"), type="numeric", default=0.14, 
              help="Fraction of random individuals"),
  make_option(c("--frac.evo"), type="numeric", default=0.85, 
              help="Fraction of evolved individuals"),
  make_option(c("--frac.eli"), type="numeric", default=0.01, 
              help="fraction of élite individuals"),
  make_option(c("--n.samples"), type="numeric", default=10, 
              help="genomes per individual"),
  make_option(c("--gen.size"), type="numeric", default=100, 
              help="population per batch"),
  make_option(c("--w.div"), type="numeric", default=0.5, 
              help="genomic diversity weight"),
  make_option(c("--w.tem"), type="numeric", default=0.5, 
              help="temporal distribution weight"),
  make_option(c("--generation"), type="numeric", default=1, 
              help="generation"),
  make_option(c("--batch"), type="numeric", default=1, 
              help="bacth"),
  make_option(c("--tot.batches"), type="numeric", default=1, 
              help="number of overall batches per generation"),
  make_option(c("--basedir"), type="character", default=getwd(), 
              help="base directory"),
  make_option(c("--dist.opt"), type="character", default="max", 
              help="Genetic distance optimized towards its max, mean, or median point"),
  make_option(c("--seeds"), type="character", default="data/seeds.txt", 
              help="file with a seeds per line/generation"),
  make_option(c("--out.dir"), type="character", default=getwd(), 
              help="output directory")
              )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

seeds = as.numeric(readLines(opt$seeds))
metadata = read.table(opt$metadata, sep = ',', stringsAsFactors = FALSE, header = TRUE)
dist.m = read.distance.matrix(opt$distance)

### CHECK METADATA
metadata = metadata[sapply(colnames(dist.m), function(x,df){which(metadata$Accession.ID == x)}, df=df),] # make sure id order is respected
metadata$Collection.date = as.Date(metadata$Collection.date, format = "%d/%m/%Y")
if(any(is.na(metadata$Collection.date))) { stop("NAs found in dates, check date format") }

cl <- parallel::makeCluster(opt$n.cores)
doParallel::registerDoParallel(cl)

set.seed(seeds[(opt$generation*opt$tot.batches)+opt$batch])

pairs.per.subsample = ((opt$n.samples^2)-opt$n.samples)/2

if(opt$dist.opt == 'max') { max.d = sum(sort(dist.m, decreasing = TRUE)[1:pairs.per.subsample]) # estimate max distance based on largest paired distances
}else{
  if(opt$dist.opt == 'mean' | opt$dist.opt == 'median') { 
    ### pairwise genetic diversity to percentile
    distances.perc = ecdf(dist.m)
    mean.perc = distances.perc(mean(dist.m))
  }
}

### Best possible time spread
first.date = min(metadata$Collection.date)
last.date = max(metadata$Collection.date)
ideal.time = seq(from = first.date,
                 to = last.date,
                 by =  (last.date-first.date)/(opt$n.samples-1))
  
### worst possible time spread
max.time.spread = as.numeric(sum(abs(ideal.time-rep(first.date, opt$n.samples))))
  
message('generating new individuals')
  
# read previous gen
if (opt$frac.evo > 0 | opt$frac.eli > 0 | opt$generation > 0){
  old.gen.fit = NULL
  old.gen.indeces = NULL
  for(i in 1:opt$tot.batches){
    tmp = as.matrix(read.csv(paste0(opt$out.dir, '/GA.', opt$data.set, '.', opt$generation-1, '.', i, '.indeces.fitness.csv'), sep = ','))[1,]
    old.gen.indeces = rbind(old.gen.indeces, as.matrix(read.table(paste0(opt$out.dir, '/GA.', opt$data.set,  '.',  opt$generation-1, '.', i, '.indeces.subsamples.csv'), sep = ',')) )
  }
  old.gen.fit = c(old.gen.fit, tmp)
  names(old.gen.fit) = as.numeric(1:length(old.gen.fit))
}

# pre allocate the random subsamples for this generation
subsamples = matrix(nrow = opt$gen.size, ncol = opt$n.samples, -1) # allocate a number of runs equal to r.tries * opt$n.samples
  
# randomly generated ones
if(opt$frac.new > 0){
  sampled.indeces = foreach(i=1: round( dim(subsamples)[1]*opt$frac.new ))  %dorng% {
    sampled.genomes = sample(metadata$Accession.ID)[1:opt$n.samples]
    sampled.indeces = which(rownames(dist.m) %in% sampled.genomes)
  }  
  subsamples[1: ( dim(subsamples)[1]*opt$frac.new),] = do.call(rbind, lapply(sampled.indeces, matrix, ncol=opt$n.samples, byrow=TRUE))
  message('randomly generated individuals ready')
}

# evolved ones
if(opt$frac.evo > 0){
  i = ( dim(subsamples)[1]*opt$frac.new )
  sampled.indeces = foreach( j=(i+1):( round(dim(subsamples)[1]*opt$frac.evo)+i )) %dorng% {
  # pick up two individuals, tournament with 5 individuals
    selected_a = sort(sample(old.gen.fit)[1:5])[5]
    selected_b = sort(sample(old.gen.fit)[1:5])[5]
    # mating
    keep.us = intersect(old.gen.indeces[as.numeric(names(selected_a)),], old.gen.indeces[as.numeric(names(selected_b)),])
    mixed = c(
              setdiff(old.gen.indeces[as.numeric(names(selected_a)),], old.gen.indeces[as.numeric(names(selected_b)),]),
              setdiff(old.gen.indeces[as.numeric(names(selected_b)),], old.gen.indeces[as.numeric(names(selected_a)),])
              )
    
    if(length(mixed) == 0){ sampled.indeces = keep.us
    }else{
      sampled.indeces = c( keep.us, sample(mixed)[1:(opt$n.samples-length(keep.us))] )       
    }
    if(runif(1) < 0.08) {sampled.indeces[sample(length(sampled.indeces))[1]]=sample(setdiff(1:dim(dist.m)[1], sampled.indeces))[1] } # mutation
    return(sampled.indeces)
  } # end of foreach
  
  sampled.indeces = do.call(rbind, lapply(sampled.indeces, matrix, ncol=opt$n.samples, byrow=TRUE))
  subsamples[( i+1 ):( (dim(subsamples)[1]*opt$frac.evo)+i ),] = sampled.indeces
  message('evolved individuals ready')
}

# elite ones
if(opt$frac.eli > 0){
  j= round(dim(subsamples)[1]*opt$frac.evo)+i
  all.elite = as.numeric(names(sort(old.gen.fit, decreasing = TRUE)[1:(opt$gen.size*opt$tot.batches*opt$frac.eli)]))
  message(paste("best fitness last gen was", old.gen.fit[all.elite[1]]))
  start.elite.pos = ((opt$batch-1)*length(all.elite)/opt$tot.batches)+1
  end.elite.pos = opt$batch*length(all.elite)/opt$tot.batches
  subsamples[(j+1):( (dim(subsamples)[1]*opt$frac.eli)+j ),] = old.gen.indeces[all.elite[start.elite.pos:end.elite.pos],]
  message('elite individuals ready')
}

# calculate fitness
output <- foreach(i=1:dim(subsamples)[1] ) %dorng% {
  r.sample = subsamples[i,] # for simplicity
  timeframe = metadata$Collection.date[r.sample]
    
  # how good this sample is in terms of diversity?
  d = combn(r.sample, 2)
  gen.diversity = 0
  for(k in 1:dim(d)[2]){ gen.diversity = gen.diversity + dist.m[d[,k][1],d[,k][2]] }
    
  # how good this sample is in terms of time spread?
  time.spread = as.numeric(sum(abs(sort(timeframe)-ideal.time)))
    
  # fitness score with rescaling
  tem.diversity = 1 - (time.spread/max.time.spread) +10^(-64)
  avg.gen.diversity = gen.diversity/dim(d)[2]
  if(opt$dist.opt == 'max') { gen.diversity = (gen.diversity/max.d) +10^(-64)  }
  if(opt$dist.opt == 'median') {
    gen.diversity = 1 - abs(distances.perc(avg.gen.diversity) - 0.5)
  }
  if(opt$dist.opt == 'mean') {
    gen.diversity = 1 - abs(distances.perc(avg.gen.diversity) - mean.perc) #
  }
  fitness = c( i, ((gen.diversity * opt$w.div + tem.diversity * opt$w.tem)/2), (avg.gen.diversity), (tem.diversity) )
}

  
message('all fitnesses calculated')
  
fitnesses = unlist(output)
# write output
out.file = paste0(opt$out.dir, '/GA.', opt$data.set, '.',  opt$generation, '.', opt$batch)

write.table(matrix(fitnesses, nrow = 4),
            file = paste(out.file, 'indeces.fitness.csv', sep = '.'),
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ',') 
    
write.table(subsamples,
            file = paste(out.file, 'indeces.subsamples.csv', sep = '.'),
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ',')

