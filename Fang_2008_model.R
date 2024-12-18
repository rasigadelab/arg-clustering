
library(data.table)
library(foreach)
library(doParallel)
library(ggplot2)
library(stringr)

# Persistence drives gene clustering in bacterial genomes
# Gang Fang 1 , Eduardo PC Rocha 1,2 and Antoine Danchin* 1
# BMC Genomics 2008
# deletion -> "random batch of contiguous genes deletion in each chromosome"
# if persistent gene deleted, cell death (chromosome does not pass to next generation)
# insertion -> "randomly picked up a position in each surviving progeny and inserted there a batch of contiguous non-peristent genes"
# steady state -> "the inadequate amount in second generation was restored by picking up genes randomly from the surviving bacteria"

#### Activate cluster ####
# Make cluster
detectCores() # number of cores of the computer
cl <- makeCluster(30) # number of cores allowed to parallel computation
registerDoParallel(cl)
# Get info
getDoParWorkers()
getDoParName()
getDoParVersion()
# Stop cluster
stopImplicitCluster()
stopCluster(cl)

#### SIMULATION ####
### PARAMETERS ###
# Fang et al 2008 parameters #
nb_bacteria <- 5000 # one chromosome each
nb_genes <- 4000 # 
nb_persistent_genes <- 400 # genes which removal leads to cell death
length_batch <- 3 # number genes in the batches of contiguous genes for insertion/deletion events
nb_generations <- 70000 # number of simulation steps

# Custom parameters - test #
nb_bacteria <- 1E5 # one chromosome each
nb_genes <- 100 #
nb_persistent_genes <- 2 # genes which removal leads to cell death
length_batch <- 3 # number genes in the batches of contiguous genes for insertion/deletion events
nb_generations <- 200 # number of simulation steps
# Custom parameters #
nb_bacteria <- 10000 # one chromosome each
nb_genes <- 100 # 
nb_persistent_genes <- 2 # genes which removal leads to cell death
length_batch <- 5 # number genes in the batches of contiguous genes for insertion/deletion events
nb_generations <- 1000 # number of simulation steps

is_essential <- TRUE # TRUE remove cells which loose an essential/persistent gene
is_control <- FALSE # TRUE transposition but no deletion/insertion (DNA part moved is conserved)
DEBUG <- FALSE

#### Mono-simulation ####
nb_iterations <- 1
for(it in 1:nb_iterations){
  ### INITIALIZATION ###
  # data.table string chromosomes #
  chromosomes <- as.data.table(matrix(ncol=1, nrow=nb_bacteria))
  setnames(chromosomes, old="V1", new="chromosome")
  chromosomes[, chromosome := as.character(chromosome)]
  chromosomes[, chromosome := random.chromosome(nb_genes, nb_persistent_genes, nb_bacteria)] # uniform istribution persistent genes
  class(chromosomes)
  
  # Results storage #
  col_names <- c("time","mean_dist","sd_dist","mean_pgene")
  data_results <- as.data.table(matrix(ncol=length(col_names), nrow=nb_generations+1))
  colnames(data_results) <- col_names
  data_results[, (col_names) := lapply(.SD, as.numeric), .SDcols=col_names]
  sapply(data_results, class)
  rm(col_names)
  
  ## Initial state t0 ##
  chromosomes[, nb_p := sapply(chromosome, function(x) str_count(x,"P"))]
  chromosomes[, dgenes := sapply(chromosome, function(x) distance.2.pgenes(x))]
  data_results[1,"time" := 0]
  data_results[1,"mean_dist" := mean(chromosomes$dgenes, na.rm=TRUE)]
  data_results[1,"sd_dist" := sd(chromosomes$dgenes)]
  data_results[1,"mean_pgene" := mean(chromosomes$nb_p)]
  
  ### SIMULATION ###
  for(i in 1:nb_generations){
    print(paste0("Iteration ",it," generation ",i,"/",nb_generations))
    if(is_control==FALSE){
      # Gene deletion #
      if(DEBUG){
        print("deletion")
      }
      chromosomes[, chromosome := sapply(chromosome, function(x) random.deletion(x,length_batch))]
      # Gene insertion #
      if(DEBUG){
        print("insertion")
      }
      chromosomes[, chromosome := sapply(chromosome, function(x) random.insertion(x,length_batch))]
    }else{
      # Gene transposition #
      if(DEBUG){
        print("transposition")
      }
      chromosomes[, chromosome := sapply(chromosome, function(x) random.transposition(x,length_batch))]
    }
    # Cell death #
    chromosomes[, nb_p := sapply(chromosome, function(x) str_count(x,"P"))]
    if(is_essential==TRUE){
      if(DEBUG){
        print("steady state")
      }
      # remove cells which lost an essential/persistent gene and repalce the lost cells with a random survivor
      nb_replacements <- nrow(chromosomes[nb_p<nb_persistent_genes]) # nb chromosomes which have deleted a persistent gene
      chromosomes[nb_p<nb_persistent_genes, chromosome:=sample(chromosomes[nb_p>=nb_persistent_genes, ]$chromosome, nb_replacements, replace=TRUE)]
      chromosomes[, nb_p := sapply(chromosome, function(x) str_count(x,"P"))] # check persistent genes
    }
    # Inter persistent gene distance / clustering indicator #
    # Kuiper's test, alpha = 0.05 if >2penes
    if(DEBUG){
      print("distance")
    }
    chromosomes[, dgenes := sapply(chromosome, function(x) distance.2.pgenes(x))]
    
    # Store mean pgene distance #
    if(DEBUG){
      print("data storage")
    }
    data_results[i+1,"time" := i]
    data_results[i+1,"mean_dist" := mean(chromosomes$dgenes, na.rm=TRUE)]
    data_results[i+1,"sd_dist" := sd(chromosomes$dgenes)]
    data_results[i+1,"mean_pgene" := mean(chromosomes$nb_p)]
  }
  data_results[, iteration := it]
  if(it==1){
    data_results_all <- data_results
  }else{
    data_results_all <- rbind(data_results_all, data_results)
  }
}

#### Parallel simulation ####
nb_iterations <- 10
data_results_all <- foreach (it=1:nb_iterations, .combine=rbind, .packages=c("stringr","data.table"), .verbose=TRUE) %dopar% {
  ### INITIALIZATION ###
  # data.table string chromosomes #
  chromosomes <- as.data.table(matrix(ncol=1, nrow=nb_bacteria))
  setnames(chromosomes, old="V1", new="chromosome")
  chromosomes[, chromosome := as.character(chromosome)]
  chromosomes[, chromosome := random.chromosome(nb_genes, nb_persistent_genes, nb_bacteria)] # uniform distribution persistent genes
  class(chromosomes)
  
  # Results storage #
  col_names <- c("time","mean_dist","sd_dist","mean_pgene")
  data_results <- as.data.table(matrix(ncol=length(col_names), nrow=nb_generations+1))
  colnames(data_results) <- col_names
  data_results[, (col_names) := lapply(.SD, as.numeric), .SDcols=col_names]
  sapply(data_results, class)
  rm(col_names)
  
  ## Initial state t0 ##
  chromosomes[, nb_p := sapply(chromosome, function(x) str_count(x,"P"))]
  chromosomes[, dgenes := sapply(chromosome, function(x) distance.2.pgenes(x))]
  data_results[1,"time" := 0]
  data_results[1,"mean_dist" := mean(chromosomes$dgenes, na.rm=TRUE)]
  data_results[1,"sd_dist" := sd(chromosomes$dgenes)]
  data_results[1,"mean_pgene" := mean(chromosomes$nb_p)]
  data_results[1,"x11" := sum(chromosomes$nb_p==2)]
  data_results[1,"x01" := sum(chromosomes$nb_p==1)]
  data_results[1,"x00" := sum(chromosomes$nb_p==0)]
  
  ### SIMULATION ###
  for(i in 1:nb_generations){
    print(paste0("Iteration ",it," generation ",i,"/",nb_generations))
    if(is_control==FALSE){
      # Gene deletion #
      if(DEBUG){
        print("deletion")
      }
      chromosomes[, chromosome := sapply(chromosome, function(x) random.deletion(x,length_batch))]
      # Gene insertion #
      if(DEBUG){
        print("insertion")
      }
      chromosomes[, chromosome := sapply(chromosome, function(x) random.insertion(x,length_batch))]
    }else{
      # Gene transposition #
      if(DEBUG){
        print("transposition")
      }
      chromosomes[, chromosome := sapply(chromosome, function(x) random.transposition(x,length_batch))]
    }
    # Cell death #
    chromosomes[, nb_p := sapply(chromosome, function(x) str_count(x,"P"))]
    if(is_essential==TRUE){
      if(DEBUG){
        print("steady state")
      }
      # remove cells which lost an essential/persistent gene and repalce the lost cells with a random survivor
      nb_replacements <- nrow(chromosomes[nb_p<nb_persistent_genes]) # nb chromosomes which have deleted a persistent gene
      chromosomes[nb_p<nb_persistent_genes, chromosome:=sample(chromosomes[nb_p>=nb_persistent_genes, ]$chromosome, nb_replacements, replace=TRUE)]
      chromosomes[, nb_p := sapply(chromosome, function(x) str_count(x,"P"))] # check persistent genes
    }
    # Inter persistent gene distance / clustering indicator #
    # Kuiper's test, alpha = 0.05 if >2penes
    if(DEBUG){
      print("distance")
    }
    chromosomes[, dgenes := sapply(chromosome, function(x) distance.2.pgenes(x))]
    
    # Store mean pgene distance #
    if(DEBUG){
      print("data storage")
    }
    data_results[i+1,"time" := i]
    data_results[i+1,"mean_dist" := mean(chromosomes$dgenes, na.rm=TRUE)]
    data_results[i+1,"sd_dist" := sd(chromosomes$dgenes, na.rm=TRUE)]
    data_results[i+1,"mean_pgene" := mean(chromosomes$nb_p)]
    data_results[i+1,"x11" := sum(chromosomes$nb_p==2)]
    data_results[i+1,"x01" := sum(chromosomes$nb_p==1)]
    data_results[i+1,"x00" := sum(chromosomes$nb_p==0)]
  } # end simulation
  data_results[,"iteration" := it]
  return(data_results)
} # end foreach


#### Plots ####
figure_path <- "./Figures"
fig_specs <- list(
  "width"=9,
  "heigth"=6,
  "dpi"=150
)

# d=f(t) #
P <- ggplot(data_results_all[x11>10])
P <- P + geom_line(aes(x=time, y=mean_dist, color=as.factor(iteration)))
P <- P + scale_x_continuous(name="generations", limits=c(0,100))
P <- P + scale_x_continuous(name="generations")
P <- P + scale_y_continuous(name="mean inter-gene distance")
P <- P + labs(color="iteration")
P <- P + theme_light()
plot(P)
ggsave("dist=f(t)_essential.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)
ggsave("dist=f(t)_unsessential.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)
ggsave("dist=f(t)_control.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)

# cell=f(t) #
color_cells <- c("green2","blue2","magenta2")
P <- ggplot(data_results_all)
P <- P + geom_line(aes(x=time, y=x00, color="x00"))
P <- P + geom_line(aes(x=time, y=x01, color="x01"))
P <- P + geom_line(aes(x=time, y=x11, color="x11"))
P <- P + scale_x_continuous(name="generations", limits=c(0,100))
P <- P + scale_x_continuous(name="generations")
P <- P + scale_y_continuous(name="cells", breaks=seq(0,1E4,2E3), labels=c(0,"2E4","4E4","6E4","8E4","1E5"))
P <- P + scale_color_manual(values=color_cells)#, breaks=c("x00","x01","x11"))
P <- P + labs(color="cell type")
P <- P + theme_light()
plot(P)
ggsave("cell=f(t)_essential.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)
ggsave("cell=f(t)_unsessential.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)
ggsave("cell=f(t)_control.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)


#### Functions ####

get.deletion <- function(nb_genes, length_batch){
  # return random consecutive deletion positions
  start <- sample(1:nb_genes, 1)
  if( (start+length_batch-1)>nb_genes ){
    del <- c(start:nb_genes, 1:(length_batch-length(start:nb_genes)))
  }else{
    del <- start:(start+length_batch-1)
  }
  if(length(del)!=length_batch){
    stop("length deletion different thant length batch")
  }
  return(del)
}
get.deletion(10,3)
class(get.deletion(10,3))
class(1:3)

random.chromosome <- function(nb_genes, nb_persistent_genes, nb_bacteria){
  # generate a random string chromosome
  all_chrom <- c()
  for(b in 1:nb_bacteria){
    chrom <- paste(rep("A", nb_genes), collapse="")
    for(i in sample(1:nb_genes, nb_persistent_genes)){
      str_sub(chrom, start=i, end=i) <- "P"
    }
    all_chrom <- c(all_chrom, chrom)
  }
  return(all_chrom)
}
# random.chromosome(20,3,5)

random.chromosome.bis <- function(nb_genes, nb_persistent_genes, nb_bact){
  # generate a random string chromosome
  all_chrom <- as.data.table(matrix(ncol=1, nrow=nb_bact))
  colnames(all_chrom) <- "chromo"
  all_chrom[, chromo:=as.character(chromo)]
  for(b in 1:nb_bact){
    chrom <- strrep("A",nb_genes)
    for(i in sample(1:nb_genes, nb_persistent_genes)){
      str_sub(chrom, start=i, end=i) <- "P"
    }
    all_chrom[b, chromo := chrom]
  }
  return(all_chrom$chromo)
}
# random.chromosome.bis(20,3,5)

random.deletion <- function(chromosome, length_batch){
  # delete a batch of genes in a string chromosome
  length_chrom <- str_length(chromosome)
  start <- sample(1:length_chrom, 1)
  if( (start+length_batch-1)>length_chrom ){ # circular chromosome
    str_sub(chromosome,start=start,end=length_chrom) <- ""
    str_sub(chromosome,start=1,end=length_batch-length(start:length_chrom)) <- ""
  }else{
    del <- start:(start+length_batch-1)
    str_sub(chromosome,start=start,end=start+length_batch-1) <- ""
  }
  return(chromosome)
}
# random.deletion("ABCDEFGHIJ",3)

random.insertion <- function(chromosome, length_batch){
  # add a batch of accessory genes "A" in a string chromosome
  length_chrom <- str_length(chromosome)
  start <- sample(1:length_chrom, 1)
  str_sub(chromosome,start=start,end=start) <- str_c(str_sub(chromosome,start=start,end=start), paste(rep("A",length_batch),collapse=""))
  return(chromosome)
}
# random.insertion("DEFGHIJ",3)

distance.2.pgenes <- function(chromosome){
  # return the number of accessory genes between 2 persistent genes
  if(nrow(str_locate_all(chromosome,"P")[[1]])==2){
    length_chrom <- str_length(chromosome)
    position1 <- str_locate_all(chromosome,"P")[[1]][1,"start"]
    position2 <- str_locate_all(chromosome,"P")[[1]][2,"start"]
    dist <- min(position2-position1-1, length_chrom-position2+position1-1)
  }else{
    dist <- NA
  }
  return(dist)
}
# distance.2.pgenes("APAPAAAAA")
# distance.2.pgenes("APAAAAPAA")
# distance.2.pgenes("APPAAPPAA")
# distance.2.pgenes("PAAAAAAAA")

random.transposition <- function(chromosome, length_batch){
  length_chrom <- str_length(chromosome)
  # deletion
  start <- sample(1:length_chrom, 1)
  if( (start+length_batch-1)>length_chrom ){ # circular chromosome
    a <- str_sub(chromosome,start=start,end=length_chrom)
    b <- str_sub(chromosome,start=1,end=length_batch-length(start:length_chrom))
    seq <- paste0(a,b)
    str_sub(chromosome,start=start,end=length_chrom) <- ""
    str_sub(chromosome,start=1,end=length_batch-length(start:length_chrom)) <- ""
  }else{
    seq <- str_sub(chromosome,start=start,end=start+length_batch-1)
    str_sub(chromosome,start=start,end=start+length_batch-1) <- ""
  }
  # insertion
  length_chrom <- str_length(chromosome) # length after deletion
  start <- sample(1:length_chrom, 1)
  str_sub(chromosome,start=start,end=start) <- str_c(str_sub(chromosome,start=start,end=start), seq)
  
  return(chromosome)
}
# random.transposition("ABCDEFGHIJ",3)

system.time({ random.chromosome(100,2,1E5) }) # 26s
system.time({ random.chromosome.bis(100,2,1E5) }) # 29s





#### Notes ####

# Kuiper test for clustering indicator (if nb_persistent_genes > 2)
# https://www.rdocumentation.org/packages/circular/versions/0.4-93/topics/kuiper.test






