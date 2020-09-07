
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