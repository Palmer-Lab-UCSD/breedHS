## UTILITY FUNCTIONS FOR PALMER LAB AIL BREEDER SELECTION CODE
## written by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

source("kinship.R")
source("find_mates.R")

encode.sex <- function(ped) {
  ped[,4] = as.character(ped[,4])
  ped[ped[,4]=="F",4] = 0  # Female = 0
  ped[ped[,4]=="M",4] = 1  # Male = 1
  ped[,4] = as.numeric(ped[,4])
  ped = data.matrix(ped)
  return(ped)
}

# identify pedigree errors 
# (for now, only identifies missing parents)
find.ped.errors <- function(first_gen,  # the desired starting generation
                            last_gen,   # the desired last generation
                            data_dir,   # the directory housing pedigree files
                            file_stem,  # the stem name shared by all pedigree files
                            write_file=TRUE,
                            return_ids=TRUE)
{

  # empty element to hold the eventual pedigree 
  ped <- c() 
    
  # empty elements to hold vectors of missing parents  
  missing_parents <- c()
  missing_gens <- c()
  missing_parents_files <- c()
  missing_kin_files <- c()

  first_nchar <- nchar(first_gen)
  last_nchar <- nchar(last_gen)
  
  ## loop across generations
  for (i in first_gen:last_gen){

    # save the file name for the previous generation
    if (i > first_gen){
        prev.file <- file.name
    }

    # get the file name for the current generation
    file.name <- file.path(data_dir, paste0(file_stem, i, ".csv"))

    if (!file.exists(file.name)) {
        gen_nchar <- nchar(i)
        n_zeros <- last_nchar - gen_nchar
        i_str <- paste0(rep('0', n_zeros), i)
        file.name <- file.path(data_dir, paste0(file_stem, i_str, ".csv"))
    }

    if (!file.exists(file.name)) {
        print(paste('Cannot find a pedigree file for generation', i, 
                   'at', file.name))
    }
    
    ## read in the pedigree for a given generation
    ped.tmp = read.table(file=file.name, sep=",", header=T,
                         as.is=TRUE, na.strings="?")
    
    ## format desired data columns
    ped.ids <- ped.tmp[,1]
    ped.cols <- c('sire', 'dam', 'sex')
    ped.tmp <- cbind(ped.ids, ped.tmp[ped.cols])
    colnames(ped.tmp) <- c('id', ped.cols)

    # remove parental data from the first generation
    if (i == first_gen){
      ped.tmp[,c(2,3)] <- 0}
    
    ## check that parents are in the previous generation
    if (i != first_gen){
      tmp <- check.ped(ped.tmp , ped.prev)   
      if (length(tmp) > 0){
        missing_parents <- c(missing_parents, tmp)
        missing_gens <- c(missing_gens, rep(i,length(tmp)))
        missing_parents_files <- c(missing_parents_files, rep(prev.file, length(tmp)))
        missing_kin_files <- c(missing_kin_files, rep(file.name, length(tmp)))
        print(paste('Generation', paste0(i, ':'), length(tmp), 'parent IDs were not found in generation', i-1))
      }
    }
    ped.prev <- ped.tmp	
    ped <- rbind(ped,ped.tmp)
  }


  if (is.null(missing_parents)) {
      print('No errors found in the pedigree')
  } else {
      ## write missing parents to a file for investigation
      if (write_file) {
          missing_df <- data.frame(
          id = missing_parents,
          birth_gen = missing_gens-1,
          kin_gen = missing_gens,
          self_ped = missing_parents_files,
          kin_ped = missing_kin_files)

          timestamp <- format(Sys.time(),'%Y%m%d-%H:%M:%S')
          outstr <- paste0(file_stem, 
                           rep('0',(last_nchar-first_nchar)),
                           first_gen, '-', last_gen,
                           '_missing_parents_', timestamp, '.csv')
          outfile <- file.path(data_dir, outstr)
          print(paste('See', outfile, 'for details'))
          write.csv(missing_df, outfile, row.names=F,quote=F)
      }
      if (return_ids) {
          return(missing_parents)        
      }
  }
}


# format multiple pedigree files (from one population) into a single pedigree
# (assumes that utils.R is in the same directory as kinship.R and find_mates.R)
format.pedigree <- function(first_gen,  # the desired starting generation
                            last_gen,   # the desired last generation
                            data_dir,   # the directory housing pedigree files
                            file_stem){ # the stem name shared by all pedigree files

  # empty element to hold the eventual pedigree 
  ped <- c() 
    
  # empty elements to hold vectors of missing parents  
  missing_parents <- c()
  missing_gens <- c()

  first_nchar <- nchar(first_gen)
  last_nchar <- nchar(last_gen)
  
  ## loop across generations
  for (i in first_gen:last_gen){

    # get the file name for the current generation
    file.name <- file.path(data_dir, paste0(file_stem, i, ".csv"))

    if (!file.exists(file.name)) {
        gen_nchar <- nchar(i)
        n_zeros <- last_nchar - gen_nchar
        i_str <- paste0(rep('0', n_zeros), i)
        file.name <- file.path(data_dir, paste0(file_stem, i_str, ".csv"))
    }

    if (!file.exists(file.name)) {
        print(paste('Cannot find a pedigree file for generation', i, 
                   'at', file.name))
    }
    
    ## read in the pedigree for a given generation
    ped.tmp = read.table(file=file.name, sep=",", header=T,
                         as.is=TRUE, na.strings="?")
    
    ## format desired data columns
    ped.ids <- ped.tmp[,1]
    ped.cols <- c('sire', 'dam', 'sex')
    ped.tmp <- cbind(ped.ids, ped.tmp[ped.cols])
    colnames(ped.tmp) <- c('id', ped.cols)

    # remove parental data from the first generation
    if (i == first_gen){
      ped.tmp[,c(2,3)] <- 0}
    
    ## check that parents are in the previous generation
    if (i != first_gen){
      tmp <- check.ped(ped.tmp , ped.prev)   
      if (length(tmp) > 0){
        print(paste('Generation', paste0(i, ':'), length(tmp), 'parent IDs were not found in generation', i-1))
      }
    }
    
    ped.prev <- ped.tmp	
    ped <- rbind(ped,ped.tmp)
      
  } # end of generation loop

  ## change sex coding from alphabetical to numeric, where F:0, M:1
  ped.tmp <- encode.sex(ped.tmp)
  ped <- encode.sex(ped)
  
  out <- list(ped.tmp = ped.tmp, ped = ped)
  return(out)
  
}


# combine pedigrees
get.ped.comb <- function(prev_pedigree, breed_results){
  first.os <- max(prev_pedigree[,1]) + 1000
  next_pedigree <- cbind(first.os:(first.os+(dim(breed_results)[1]-1)),
                         breed_results[,1:2], rep(1,dim(breed_results)[1]))
  prev_pedigree <- as.data.frame(prev_pedigree, )
  next_pedigree <- as.data.frame(next_pedigree, )
  names(next_pedigree)[1:4] <- c("ID","Father","Mother","Sex")
  names(prev_pedigree) <- c("ID","Father","Mother","Sex")
  ped.comb <- rbind(prev_pedigree, next_pedigree)
  out <- list(ped.comb=ped.comb, ped.prev=prev_pedigree, ped.next=next_pedigree)
  return(out)
}


# wrapper function to conduct all breedail breeder selection
select.breeders <- function(first_gen,              # first generation of the pedigree
                            last_gen,               # final generation of the pedigree
                            data_dir,               # directory housing pedigree files
                            out_dir,                # the output directory in which to save results
                            file_stem,              # stem name for pedigree files
                            n_rounds,               # number of rounds of breeder selection desired
                            one_per_sibship = TRUE) # whether one (T) or multiple breeders per sibship
{
  
  # format pedigree files into one pedigree
  ped.out <- format.pedigree(first_gen, last_gen, data_dir, file_stem)
  ped.tmp <- ped.out$ped.tmp
  ped <- ped.out$ped
  
  ## build the covariance matrix
  k <- kinship(ped)
  
  ## get dimensions from kinship matrix
  N <- dim(k)[1]
  no <- dim(ped.tmp)[1]
  idx <- (N-no+1):N
  k.prev <- k[idx,idx]
  ped.prev <- ped[idx,]
  
  # first mate pairing
  if (one_per_sibship){
    
    sibs <- 'one_per_sibship'
    first.breeders <- find.mates.res(ped.prev,k.prev)
    colnames(first.breeders) <- c('sire', 'dam', 'kinship')
  } else {
    
    sibs <- 'multi_per_sibship'
    first.breeders <- find.mates(ped.prev,k.prev)
  }
  
  round.1 <- rep(1, nrow(first.breeders))
  
  if (n_rounds==1){
    first.breeders <- cbind(round.1, first.breeders)
    write.table(first.breeders, 
                file.path(out_dir, paste0('breedpairs_F', last_gen, '_n1_', 
                                          sibs, '.txt')),
                quote=F, row.names=F, sep='\t')
    return(first.breeders)
  }
  
  ped.comb.out <- get.ped.comb(ped.prev, first.breeders)
  ped.comb <- ped.comb.out$ped.comb
  ped.prev <- ped.comb.out$ped.prev
  ped.next <- ped.comb.out$ped.next
  
  ## update kinship matrix ##
  k.comb <- kinship.update(ped.prev , k.prev, ped.next)
  pos.prev <- 1:dim(ped.prev)[1]
  pos.next <- dim(ped.prev)[1]+(1:dim(ped.next)[1])
  not.avail.idx <- ped.prev[,1] %in% c(ped.next[,2],ped.next[,3])
  not.avail.pos <- c(1:dim(ped.prev)[1])[not.avail.idx]
  
  # second mate pairing
  new.mates <- find.mates.given.pop(ped.comb, k.comb, pos.prev, pos.next, not.avail.pos, max.pairs=1)
  breeders <- rbind(first.breeders, new.mates[,1:3])

  if (!is.null(new.mates)){
    round.2 <- rep(2, nrow(new.mates))
  } else(round.2 <- NULL)
  
  if (n_rounds==2){
    
    round <- c(round.1, round.2)
    breeders <- cbind(round, breeders)  
    write.table(breeders,
                file.path(out_dir, paste0('breedpairs_F', last_gen, '_n2_', 
                                          sibs, '.txt')),
                quote=F, row.names=F, sep='\t')
    return(breeders)
  }
  
  else if (n_rounds > 2){
    
    round.i <- list()
    all.mates <- list(first.breeders, new.mates[,1:3])
    
    for (i in 3:n_rounds){
      
      # dataframe of all breeder pairs so far
      all.breeders <- do.call(rbind, all.mates[1:i-1])

      ped.comb.out <- get.ped.comb(ped.prev, all.breeders)
      ped.comb <- ped.comb.out$ped.comb
      ped.prev <- ped.comb.out$ped.prev
      ped.next <- ped.comb.out$ped.next
      
      ## update kinship matrix 
      k.comb <- kinship.update(ped.prev , k.prev, ped.next)
      pos.prev <- 1:dim(ped.prev)[1]
      pos.next <- dim(ped.prev)[1]+(1:dim(ped.next)[1])
      not.avail.idx <- ped.prev[,1] %in% c(ped.next[,2],ped.next[,3])
      not.avail.pos <- c(1:dim(ped.prev)[1])[not.avail.idx]
      
      # next round of mate pairing
      new.mates <- find.mates.given.pop(ped.comb, k.comb, 
                                        pos.prev, pos.next, not.avail.pos)
      

      all.mates[[i]] <- new.mates[,1:3]
      round.i[[i]] <- rep(i, nrow(new.mates))
      
    }
    
    all.breeders <- do.call(rbind, all.mates)
    round <- c(round.1, round.2, unlist(round.i))
    all.breeders <- cbind(round, all.breeders)
    write.table(all.breeders,
                file.path(out_dir, paste0('breedpairs_F', last_gen, '_n', n_rounds, 
                                          '_', sibs, '.txt')),
                quote=F, row.names=F, sep='\t')
    
    return(all.breeders)
    
    
  }
  
  
}

# simulate breeding between mate pairs
mate.breeders <- function(breedpairs,   # matrix output by select.breeders, or path to file
                          n_pairs=NULL, # number of pairs to breed, default is all pairs
                          n_sibs,       # number of offspring per pair
                          outfile=NULL)      # file to save new gen pedigree
{     
  
  if (class(breedpairs)[1]=='matrix'){
    breedpairs <- as.data.frame(breedpairs)
  } else if (class(breedpairs)[1]=='character'){
    breedpairs <- read.table(breedpairs, sep='\t', header=T)
  }
  
  if (!is.null(n_pairs)){
    sample_idx <- sort(sample(1:nrow(breedpairs), n_pairs))
    dam <- breedpairs$dam[sample_idx]
    dam <- rep(dam, each=n_sibs)
    sire <- breedpairs$sire[sample_idx]
    sire <- rep(sire, each=n_sibs)
  } else {
    dam <- rep(breedpairs$dam, each=n_sibs)
    sire <- rep(breedpairs$sire, each=n_sibs)
  }
  
  sex <- rep(c('F','M'), length.out=length(dam))
  first_id <- max(c(breedpairs$dam, breedpairs$sire)) + 1
  last_id <- max(c(breedpairs$dam, breedpairs$sire)) + length(dam)
  id <- seq.int(first_id,last_id)
  
  new_gen <- data.frame(id, sex, dam, sire)
  
  if (!is.null(outfile)){
    write.csv(new_gen, outfile, row.names=F, quote=F)
  }
  return(new_gen)
}

# format the desired columns from a pedigree dataframe
format.ped.cols <- function(df){
  
  ## Format desired data columns
  ped.ids <- df[,1]
  ped.cols <- c('sire', 'dam', 'sex')
  ped <- cbind(ped.ids, df[ped.cols])
  colnames(ped) <- c('id', ped.cols)
  return(ped)
}


# get breeding pairs for a breeder exchange between two colonies
# population 1 must be the original founder of pop 2
exchange.breeders <- function(
        dir_1,                # directory of ped files for population 1
        stem_1,               # file name stem for pop 1 pedigree files
        first_gen_1,          # number of the first generation analyzed from pop 1
        last_gen_1,           # number of the final (exchange) generation from pop 1
        dir_2,                # directory of ped files for population 2
        stem_2,               # file name stem for pop 2 pedigree files
        first_gen_2,          # number of the first generation analyzed from pop 2
        last_gen_2,           # number of the final (exchange) generation from pop 2
        out_dir,              # directory in which to store output files
        out_stem,             # file name stem of output files
        one_per_sibship=TRUE) # boolean: one or multiple breeders per sibship
{
    
    ped1_all <- format.pedigree(first_gen_1, last_gen_1-1, dir_1, stem_1)$ped      
    ped2_all <- format.pedigree(first_gen_2, last_gen_2-1, dir_2, stem_2)$ped      
    ped1_now <- format.pedigree(first_gen_1, last_gen_1, dir_1, stem_1)$ped.tmp   
    ped2_now <- format.pedigree(first_gen_2, last_gen_2, dir_2, stem_2)$ped.tmp 
    
    # remove duplicates from hsw
    ped2_all <- ped2_all[!ped2_all[,1] %in% ped1_all[,1],]
    merged_ped <- rbind(ped1_all, ped2_all)
    merged_ped <- merged_ped[order(merged_ped[,1]),]
    
    # artificially create a final generation by merging potential breeders
    # from either population
    m1 <- ped1_now[ped1_now[,4]==1,]
    f1 <- ped1_now[ped1_now[,4]==0,]
    m2 <- ped2_now[ped2_now[,4]==1,]
    f2 <- ped2_now[ped2_now[,4]==0,]
    
    m1f2_now <- rbind(m1, f2) # wfu M, hsw F final gen
    m2f1_now <- rbind(m2, f1) # hsw M, wfu F final gen
    
    ped1 <- rbind(merged_ped, m1f2_now)
    ped2 <- rbind(merged_ped, m2f1_now)
    
    # estimate kinship and get dimensions from kinship matrix
    k1 <- kinship(ped1)
    k2 <- kinship(ped2)
    
    N1 <- dim(k1)[1]           # total number of rats in the pop1 pedigree
    no1 <- dim(m1f2_now)[1]    # number of rats in final generation
    idx1 <- (N1-no1+1):N1      # index positions for final gen rats
    k.prev1 <- k1[idx1,idx1]   # kinship matrix of previous generations
    ped.prev1 <- ped1[idx1,]   # pedigree of just final gen
    
    N2 <- dim(k2)[1]
    no2 <- dim(m2f1_now)[1]
    idx2 <- (N2-no2+1):N2
    k.prev2 <- k2[idx2,idx2]
    ped.prev2 <- ped2[idx2,]
    
    # first mate pairing
    if (one_per_sibship){
        round1.breeders1 <- find.mates.res(ped.prev1,k.prev1)
        round1.breeders2 <- find.mates.res(ped.prev2,k.prev2)
    } else {
        round1.breeders1 <- find.mates(ped.prev1,k.prev1)
        round1.breeders2 <- find.mates(ped.prev2,k.prev2)
    }
    
    round1.idx1 <- rep(1, nrow(round1.breeders1))
    round1.idx2 <- rep(1, nrow(round1.breeders2))
    
    all.breeders1 <- cbind(round1.idx1, round1.breeders1)
    all.breeders2 <- cbind(round1.idx2, round1.breeders2)
    
    colnames(all.breeders1) <- c('round', 'sire', 'dam', 'kinship')
    colnames(all.breeders2) <- c('round', 'sire', 'dam', 'kinship')
    
    # update pedigrees
    first.os1 <- max(ped.prev1[,1])+9999
    ped.next1 <- cbind(first.os1:(first.os1+(dim(round1.breeders1)[1]-1)),
                       round1.breeders1[,1:2], rep(1,dim(round1.breeders1)[1]))
    ped.prev1 <- as.data.frame(ped.prev1, )
    ped.next1 <- as.data.frame(ped.next1, )
    ped.cols <- c('ID', 'Father', 'Mother', 'Sex')
    names(ped.next1) <- ped.cols
    names(ped.prev1) <- ped.cols
    ped.comb1 <- rbind(ped.prev1, ped.next1)
    
    first.os2 <- max(ped.prev2[,1])+1000
    ped.next2 <- cbind(first.os2:(first.os2+(dim(round1.breeders2)[1]-1)),
                       round1.breeders2[,1:2], rep(1,dim(round1.breeders2)[1]))
    ped.prev2 <- as.data.frame(ped.prev2, )
    ped.next2 <- as.data.frame(ped.next2, )
    names(ped.next2)[1:4] <- ped.cols
    names(ped.prev2) <- ped.cols
    ped.comb2 <- rbind(ped.prev2, ped.next2)
    
    # update kinship matrices
    k.comb1 <- kinship.update(ped.prev1 , k.prev1, ped.next1)
    pos.prev1 <- 1:dim(ped.prev1)[1]
    pos.next1 <- dim(ped.prev1)[1]+(1:dim(ped.next1)[1])
    not.avail.idx1 <- ped.prev1[,1] %in% c(ped.next1[,2],ped.next1[,3])
    not.avail.pos1 <- c(1:dim(ped.prev1)[1])[not.avail.idx1]
    
    k.comb2 <- kinship.update(ped.prev2 , k.prev2, ped.next2)
    pos.prev2 <- 1:dim(ped.prev2)[1]
    pos.next2 <- dim(ped.prev2)[1]+(1:dim(ped.next2)[1]) 
    not.avail.idx2 <- ped.prev2[,1] %in% c(ped.next2[,2],ped.next2[,3])
    not.avail.pos2 <- c(1:dim(ped.prev2)[1])[not.avail.idx2]
    
    # second mate pairing 
    round2.breeders1 <- find.mates.given.pop(ped.comb1, k.comb1, 
                                             pos.prev1, pos.next1, not.avail.pos1)
    round2.breeders2 <- find.mates.given.pop(ped.comb2, k.comb2, 
                                             pos.prev2, pos.next2, not.avail.pos2)
    
    if (!is.null(round2.breeders1)){
        
        round2.idx1 <- rep(2, nrow(round2.breeders1))
        round2.breeders1 <- cbind(round2.idx1, round2.breeders1[,1:3])
        all.breeders1 <- rbind(all.breeders1, round2.breeders1)
    } 
    
    if (!is.null(round2.breeders2)){
        
        round2.idx2 <- rep(2, nrow(round2.breeders2))
        round2.breeders2 <- cbind(round2.idx2, round2.breeders2[,1:3])
        all.breeders2 <- rbind(all.breeders2, round2.breeders2)
    } 
    
    # update pedigrees
    breed.rslts1 <- rbind(round1.breeders1 , round2.breeders1[,2:4])
    ped.next1 <- cbind(first.os1:(first.os1+(dim(breed.rslts1)[1]-1)),
                       breed.rslts1[,1:2], rep(1,dim(breed.rslts1)[1]))
    ped.prev1 <- as.data.frame(ped.prev1, )
    ped.next1 <- as.data.frame(ped.next1, )
    names(ped.next1)[1:4] <- ped.cols
    names(ped.prev1) <- ped.cols
    ped.comb1 <- rbind(ped.prev1, ped.next1)
    
    breed.rslts2 <- rbind(round1.breeders2 , round2.breeders2[,2:4])
    ped.next2= cbind(first.os2:(first.os2+(dim(breed.rslts2)[1]-1)),
                     breed.rslts2[,1:2], rep(1,dim(breed.rslts2)[1]))
    ped.prev2 <- as.data.frame(ped.prev2, )
    ped.next2 <- as.data.frame(ped.next2, )
    names(ped.next2)[1:4] <- ped.cols
    names(ped.prev2) <- ped.cols
    ped.comb2 <- rbind(ped.prev2, ped.next2)
    
    # update kinship matrices
    k.comb1 <- kinship.update(ped.prev1 , k.prev1, ped.next1)
    pos.prev1 <- 1:dim(ped.prev1)[1]
    pos.next1 <- dim(ped.prev1)[1]+(1:dim(ped.next1)[1])
    not.avail.idx1 <- ped.prev1[,1] %in% c(ped.next1[,2],ped.next1[,3])
    not.avail.pos1 <- c(1:dim(ped.prev1)[1])[not.avail.idx1]
    
    k.comb2 <- kinship.update(ped.prev2 , k.prev2, ped.next2)
    pos.prev2 <- 1:dim(ped.prev2)[1]
    pos.next2 <- dim(ped.prev2)[1]+(1:dim(ped.next2)[1])
    not.avail.idx2 <- ped.prev2[,1] %in% c(ped.next2[,2],ped.next2[,3])
    not.avail.pos2 <- c(1:dim(ped.prev2)[1])[not.avail.idx2]
    
    # third mate pairing
    round3.breeders1 <- find.mates.given.pop(ped.comb1, k.comb1, 
                                             pos.prev1, pos.next1, not.avail.pos1)
    round3.breeders2 <- find.mates.given.pop(ped.comb2, k.comb2, 
                                             pos.prev2, pos.next2, not.avail.pos2)
    
    if (!is.null(round3.breeders1)){
        
        round3.idx1 <- rep(3, nrow(round3.breeders1))
        round3.breeders1 <- cbind(round3.idx1, round3.breeders1[,1:3])
        all.breeders3 <- rbind(all.breeders1, round3.breeders1)
    } 
    
    if (!is.null(round3.breeders2)){
        
        round3.idx2 <- rep(3, nrow(round3.breeders2))
        round3.breeders2 <- cbind(round3.idx2, round3.breeders2[,1:3])
        all.breeders2 <- rbind(all.breeders2, round3.breeders2)
    } 
    
    
    # update pedigrees
    breed.rslts1 <- rbind(round1.breeders1 , round2.breeders1[,2:4], round3.breeders1[,2:4])
    ped.next1 <- cbind(first.os1:(first.os1+(dim(breed.rslts1)[1]-1)),
                       breed.rslts1[,1:2], rep(1,dim(breed.rslts1)[1]))
    ped.prev1 <- as.data.frame(ped.prev1, )
    ped.next1 <- as.data.frame(ped.next1, )
    names(ped.next1)[1:4] <- ped.cols
    names(ped.prev1) <- ped.cols
    ped.comb1 <- rbind(ped.prev1, ped.next1)
    
    breed.rslts2 <- rbind(round1.breeders2 , round2.breeders2[,2:4], round3.breeders2[,2:4])
    ped.next2= cbind(first.os2:(first.os2+(dim(breed.rslts2)[1]-1)),
                     breed.rslts2[,1:2], rep(1,dim(breed.rslts2)[1]))
    ped.prev2 <- as.data.frame(ped.prev2, )
    ped.next2 <- as.data.frame(ped.next2, )
    names(ped.next2)[1:4] <- ped.cols
    names(ped.prev2) <- ped.cols
    ped.comb2 <- rbind(ped.prev2, ped.next2)
    
    # update kinship matrices
    k.comb1 <- kinship.update(ped.prev1 , k.prev1, ped.next1)
    pos.prev1 <- 1:dim(ped.prev1)[1]
    pos.next1 <- dim(ped.prev1)[1]+(1:dim(ped.next1)[1])
    not.avail.idx1 <- ped.prev1[,1] %in% c(ped.next1[,2],ped.next1[,3])
    not.avail.pos1 <- c(1:dim(ped.prev1)[1])[not.avail.idx1]
    
    k.comb2 <- kinship.update(ped.prev2 , k.prev2, ped.next2)
    pos.prev2 <- 1:dim(ped.prev2)[1]
    pos.next2 <- dim(ped.prev2)[1]+(1:dim(ped.next2)[1])
    not.avail.idx2 <- ped.prev2[,1] %in% c(ped.next2[,2],ped.next2[,3])
    not.avail.pos2 <- c(1:dim(ped.prev2)[1])[not.avail.idx2]
    
    # fourth mate pairing
    round4.breeders1 <- find.mates.given.pop(ped.comb1, k.comb1, 
                                             pos.prev1, pos.next1, not.avail.pos1)
    
    round4.breeders2 <- find.mates.given.pop(ped.comb2, k.comb2, 
                                             pos.prev2, pos.next2, not.avail.pos2)
    
    if (!is.null(round4.breeders1)){
        
        round4.idx1 <- rep(4, nrow(round4.breeders1))
        round4.breeders1 <- cbind(round4.idx1, round4.breeders1[,1:3])
        all.breeders1 <- rbind(all.breeders1, round4.breeders1)
    } 
    
    if (!is.null(round4.breeders2)){
        
        round4.idx2 <- rep(4, nrow(round4.breeders2))
        round4.breeders2 <- cbind(round4.idx2, round4.breeders2[,1:3])
        all.breeders2 <- rbind(all.breeders2, round4.breeders2)
    } 
    

    if (one_per_sibship){
        
        write.table(all.breeders1, file.path(out_dir, paste0(out_stem, '_M1F2_breedpairs_one_per_sibship.txt')),
                    sep='\t', row.names=F, quote=F)
        write.table(all.breeders2, file.path(out_dir, paste0(out_stem, '_M2F1_breedpairs_one_per_shibship.txt')),
                    sep='\t', row.names=F, quote=F)
        
    } else {
        
        write.table(all.breeders1, file.path(out_dir, paste0(out_stem, '_M1F2_breedpairs_multi_per_sibship.txt')),
                    sep='\t', row.names=F, quote=F)
        write.table(all.breeders2, file.path(out_dir, paste0(out_stem, '_M2F1_breedpairs_multi_per_sibship.txt')),
                    sep='\t', row.names=F, quote=F)
        
    }
    
    return(list(M1F2 = all.breeders1, M2F1 = all.breeders2))
    
}

# create a named vector of all generations (from a single population) involved in sample exchanges
name.exchanges <- function(sent,     # vector of generations sent to the other population
                           received) # vector of generations that received from the other population
{

    all_exchanges <- sort(unique(c(sent, received)))
    events <- c()
    
    for (event in all_exchanges) {
        if ((event %in% sent) & !(event %in% received)) {
            event_type <- 'sent'
        }
        else if (!(event %in% sent) & (event %in% received)) {
            event_type <- 'received'
        }
        else if ((event %in% sent) & (event %in% received)) {
            event_type <- 'exchanged'
        }
        events <- c(events, event_type)
    }
    
    names(all_exchanges) <- events
    return(all_exchanges)

}


# format multiple pedigree files (from one population) to prep them 
# for merging with pedigree files from another population
format.pedigree.for.merge <- function(
                            first_gen,  # the desired starting generation
                            last_gen,   # the desired last generation
                            data_dir,   # the directory housing pedigree files
                            file_stem)  # the stem name shared by all pedigree files
{
    
    # empty element to hold the eventual pedigree 
    ped <- c() 
      
    first_nchar <- nchar(first_gen)
    last_nchar <- nchar(last_gen)
    
    ## loop across generations
    for (i in first_gen:last_gen){
        
  
        # get the file name for the current generation
        file.name <- file.path(data_dir, paste0(file_stem, i, ".csv"))
  
        if (!file.exists(file.name)) {
            gen_nchar <- nchar(i)
            n_zeros <- last_nchar - gen_nchar
            i_str <- paste0(rep('0', n_zeros), i)
            file.name <- file.path(data_dir, paste0(file_stem, i_str, ".csv"))
        }
    
        if (!file.exists(file.name)) {
            print(paste('Cannot find a pedigree file for generation', i, 
                       'at', file.name))
        }
        
        ## read in the pedigree for a given generation
        ped.tmp = read.table(file=file.name, sep=",", header=T,
                             as.is=TRUE, na.strings="?")
        
        ## format desired data columns
        ped.ids <- ped.tmp[,1]
        ped.cols <- c('sire', 'dam', 'sex', 'generation')
        ped.tmp <- cbind(ped.ids, ped.tmp[ped.cols])
        colnames(ped.tmp) <- c('id', ped.cols)

        # order by ID
        ped.tmp <- as.data.frame(ped.tmp)
        ped.tmp <- as.matrix(ped.tmp[order(ped.tmp[,1]),])

    
        # remove parental data from the first generation
        if (i == first_gen){
            ped.tmp[,c(2,3)] <- 0}
        
        ped <- rbind(ped,ped.tmp)  
        
    } # end of generation loop
  
    # change sex coding from alphabetical to numeric, where F:0, M:1
    ped <- encode.sex(ped)
    ped <- as.data.frame(ped)
    
    # split pedigree by generations 
    for (col in colnames(ped)) {
        ped[[col]] <- as.numeric(ped[[col]])
    }
    ped$generation <- as.numeric(ped$generation)
    ped_gens <- split(ped, ped$generation)
  
      
    out <- ped_gens
    return(out)
  
}
