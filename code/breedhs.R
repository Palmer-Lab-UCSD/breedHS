## UTILITY FUNCTIONS FOR PALMER LAB AIL BREEDER SELECTION CODE
## written by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

library(tools)
source("kinship.R")
source("find_mates.R")
source("hsw_utils.R")
source("wfu_utils.R")

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
                            return_ids=TRUE,
                            print_ids=TRUE)
{

  # empty element to hold the eventual pedigree 
  ped <- c() 
    
  # empty elements to hold vectors of missing parents  
  missing_parents <- c()
  missing_gens <- c()
  missing_parents_files <- c()
  missing_kin_files <- c()
  empty_parents <- list()
  kin_no_parents <- list()

  first_nchar <- nchar(first_gen)
  last_nchar <- nchar(last_gen)
  
  ## loop across generations
  for (i in first_gen:last_gen){

    # save the file name for the previous generation
    if (i > first_gen){
        prev.file <- file.name
    }

    # get the file name for the current generation
    file.name <- file.path(data_dir, paste0(file_stem, i, '.csv'))

    if (!file.exists(file.name)) {
        gen_nchar <- nchar(i)
        n_zeros <- last_nchar - gen_nchar
        i_str <- paste0(rep('0', n_zeros), i)
        file.name <- file.path(data_dir, paste0(file_stem, i_str, '.csv'))
    }

    if (!file.exists(file.name)) {
        cat('Cannot find a pedigree file for generation', i, 'at', file.name, '\n')
    }
    
    ## read in the pedigree for a given generation
    ped.tmp = read.table(file=file.name, sep=',', header=T,
                         as.is=TRUE, na.strings=c('','?','NA'))
    
    ## format desired data columns
    ped.ids <- ped.tmp[,1]
    ped.cols <- c('sire', 'dam', 'sex')
    ped.tmp <- cbind(ped.ids, ped.tmp[ped.cols])
    colnames(ped.tmp) <- c('id', ped.cols)

    # remove parental data from the first generation
    if (i == first_gen){
      ped.tmp[,c(2,3)] <- 0}
    
    if (i != first_gen){
      ## check that parents are in the previous generation
      tmp <- check.ped(ped.tmp , ped.prev)   
      if (length(tmp) > 0){
        missing_parents <- c(missing_parents, tmp)
        missing_gens <- c(missing_gens, rep(i,length(tmp)))
        missing_parents_files <- c(missing_parents_files, rep(prev.file, length(tmp)))
        missing_kin_files <- c(missing_kin_files, rep(file.name, length(tmp)))

        # get IDs of offspring whose parents are missing
        dams <- unique(ped.tmp[,2])
        sires <- unique(ped.tmp[,3])
        dams <- dams[which(dams != 0)]
        sires <- sires[which(sires != 0)]
        missing_dams_idx <- which(!dams %in% ped.prev[,1])
        missing_sires_idx <- which(!sires %in% ped.prev[,1])
        missing_idx <- union(missing_dams_idx, missing_sires_idx)
        missing_dams <- dams[missing_dams_idx]
        missing_sires <- sires[missing_sires_idx]

        if (length(missing_idx) > 0) {

          kin_missing_parents <- data.frame(
            id = ped.tmp[,1][missing_idx],
            gen = i,
            parent_gen = i-1,
            missing_dam = NA,
            missing_sire = NA)
        
          for (r in 1:nrow(kin_missing_parents)) {
            id <- kin_missing_parents$id[r]
            dam <- ped.tmp[,2][ped.tmp[,1] == id]
            sire <- ped.tmp[,3][ped.tmp[,1] == id]
            if (dam %in% missing_dams) {
              kin_missing_parents$missing_dam[r] <- dam
            }
            if (sire %in% missing_sires) {
              kin_missing_parents$missing_sire[r] <- sire
            }
          }
          kin_no_parents[[i]] <- kin_missing_parents
        }

        cat('Generation', paste0(i, ':'), length(tmp), 'parent IDs were not found in generation', i-1, '\n')
        if (print_ids) {
          cat('\t', sort(tmp), '\n')
        }
      }

      ## check for missing values
      parent_nans <- ped.tmp[(is.na(ped.tmp$dam)) | (is.na(ped.tmp$sire)),]
      if(nrow(parent_nans) > 0) {
        parent_nans$generation <- i
        parent_nans <- parent_nans[,c('generation','id','dam','sire','sex')]
        empty_parents[[i]] <- parent_nans
      }
    }
    ped.prev <- ped.tmp	
    ped <- rbind(ped,ped.tmp)
  }
  
  empty_parents <- do.call(rbind, empty_parents)
  kin_no_parents <- do.call(rbind, kin_no_parents)

  if ((is.null(missing_parents)) & (length(empty_parents == 0))){
      cat('No errors found in the pedigree \n')
  } else {
      ## write missing parents to a file for investigation
      if (write_file) {
          
          cat('See file(s) for details: \n')

          if (!is.null(missing_parents)) {

            missing_df <- data.frame(
            id = missing_parents,
            gen = missing_gens-1,
            kin_gen = missing_gens,
            ped = missing_parents_files,
            kin_ped = missing_kin_files)

            datestamp <- format(Sys.time(),'%Y%m%d')
            outstr1 <- paste0(file_stem, '_',
                            rep('0',(last_nchar-first_nchar)),
                            first_gen, '_', last_gen,
                            '_ped_error_missing_parents_', datestamp, '.csv')
            outfile1 <- file.path(data_dir, outstr1)
            write.csv(missing_df, outfile1, row.names=F, quote=F)
            cat(outfile1, '\n')

            outstr2 <- paste0(file_stem, '_',
                            rep('0',(last_nchar-first_nchar)),
                            first_gen, '_', last_gen,
                            '_ped_error_kin_w_missing_parents_', datestamp, '.csv')
            outfile2 <- file.path(data_dir, outstr2)
            write.csv(kin_no_parents, outfile2, row.names=F, quote=F)
            cat(outfile2, '\n')
          }

          if (length(empty_parents == 0)) {

            outstr3 <- paste0(file_stem, '_',
                            rep('0',(last_nchar-first_nchar)),
                            first_gen, '_', last_gen,
                            '_ped_error_parent_NAs_', datestamp, '.csv')
            outfile3 <- file.path(data_dir, outstr3)
            write.csv(empty_parents, outfile3, row.names=F, quote=F)
            cat(outfile3, '\n')
          }
      }
      if (return_ids) {
          return(list(missing_ids = missing_ids,
                      empty_parents = empty_parents))
          return(empty_parents)
      }
  }
}


# format multiple pedigree files (from one population) into a single pedigree
format.pedigree <- function(first_gen,  # the desired starting generation
                            last_gen,   # the desired last generation
                            data_dir,   # the directory housing pedigree files
                            file_stem, # the stem name shared by all pedigree files
                            print_errors = TRUE)
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
        cat('Cannot find a pedigree file for generation', i, 'at', file.name, '\n')
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
    
    if (i != first_gen){
      ## drop rows with NA values for parents
      # ped.tmp <- ped.tmp[complete.cases(ped.tmp),]
      if (print_errors){
        # check that parents are in the previous generation
        tmp <- check.ped(ped.tmp , ped.prev)   
        if (length(tmp) > 0){
          cat('Generation', paste0(i, ':'), length(tmp), 'parent IDs were not found in generation', i-1, '\n')
          cat('\t', sort(tmp), '\n')
        }
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
                            one_per_sibship = TRUE) # whether one (T) or multiple breeders per sibship
{
  MAX_ROUNDS <- 10  # Maximum number of rounds to attempt
  
  # format pedigree files into one pedigree
  ped.out <- format.pedigree(first_gen, last_gen, data_dir, file_stem, print_errors=T)
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
  } else {
    sibs <- 'multi_per_sibship'
    first.breeders <- find.mates(ped.prev,k.prev)
  }

  # note: breedail pedigree columns are ordered sire/dam
  # but find.mates outputs are ordered dam/sire
  colnames(first.breeders) <- c('dam', 'sire', 'kinship')
  round.1 <- rep(1, nrow(first.breeders))
  round <- round.1
  
  ped.comb.out <- get.ped.comb(ped.prev, first.breeders)
  ped.comb <- ped.comb.out$ped.comb
  ped.prev <- ped.comb.out$ped.prev
  ped.next <- ped.comb.out$ped.next
  
  ## update kinship matrix ##
  k.comb <- kinship.update(ped.prev, k.prev, ped.next)
  pos.prev <- 1:dim(ped.prev)[1]
  pos.next <- dim(ped.prev)[1]+(1:dim(ped.next)[1])
  not.avail.idx <- ped.prev[,1] %in% c(ped.next[,2],ped.next[,3])
  not.avail.pos <- c(1:dim(ped.prev)[1])[not.avail.idx]
  
  # second mate pairing (now using default max_pairs)
  new.mates <- find.mates.given.pop(ped.comb, k.comb, pos.prev, pos.next, not.avail.pos)
  
  # Initialize storage for all rounds
  all.mates <- list()
  round.numbers <- list()
  
  all.mates[[1]] <- first.breeders
  round.numbers[[1]] <- round.1
  
  if (!is.null(new.mates)) {
    all.mates[[2]] <- new.mates[,1:3]
    round.numbers[[2]] <- rep(2, nrow(new.mates))
  }
  
  # Try additional rounds (3 onwards)
  for (i in 3:MAX_ROUNDS) {
    tryCatch({
      # Get combined pedigree for all previous rounds
      all.breeders <- do.call(rbind, all.mates[1:(i-1)])
      
      ped.comb.out <- get.ped.comb(ped.prev, all.breeders)
      ped.comb <- ped.comb.out$ped.comb
      ped.prev <- ped.comb.out$ped.prev
      ped.next <- ped.comb.out$ped.next
      
      ## update kinship matrix
      k.comb <- kinship.update(ped.prev, k.prev, ped.next)
      pos.prev <- 1:dim(ped.prev)[1]
      pos.next <- dim(ped.prev)[1]+(1:dim(ped.next)[1])
      not.avail.idx <- ped.prev[,1] %in% c(ped.next[,2],ped.next[,3])
      not.avail.pos <- c(1:dim(ped.prev)[1])[not.avail.idx]
      
      # next round of mate pairing
      new.mates <- find.mates.given.pop(ped.comb, k.comb, 
                                      pos.prev, pos.next, not.avail.pos)
      
      if (!is.null(new.mates) && nrow(new.mates) > 0) {
        all.mates[[i]] <- new.mates[,1:3]
        round.numbers[[i]] <- rep(i, nrow(new.mates))
      } else {
        break
      }
      
    }, error = function(e) {
      cat(sprintf("Round %d: %s\n", i, conditionMessage(e)))
    })
    
    # If no new mates were found or an error occurred, stop trying new rounds
    if (is.null(all.mates[[i]])) {
      break
    }
  }
  
  # Process all successful rounds
  all.breeders <- do.call(rbind, all.mates)
  round <- unlist(round.numbers)
  all.breeders <- as.data.frame(cbind(round, all.breeders))
  all.breeders <- all.breeders[order(all.breeders$kinship),]
  
  timestamp <- format(Sys.time(),'%Y%m%d-%H:%M:%S')
  outfile <- file.path(out_dir, paste0('breedpairs_F', last_gen, '_n', length(all.mates), 
                                      '_', sibs, '_', timestamp, '.csv'))
  write.csv(all.breeders, outfile, quote=F, row.names=F)
  cat('Successfully paired', nrow(all.breeders), 'breeder pairs \n')
  cat('Pairing file written to', outfile, '\n')
  
  return(list(pairs = all.breeders, file = outfile))
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
    
    ped1_all <- format.pedigree(first_gen_1, last_gen_1-1, dir_1, stem_1, print_errors = F)$ped      
    ped2_all <- format.pedigree(first_gen_2, last_gen_2-1, dir_2, stem_2, print_errors = F)$ped      
    ped1_now <- format.pedigree(first_gen_1, last_gen_1, dir_1, stem_1, print_errors = F)$ped.tmp   
    ped2_now <- format.pedigree(first_gen_2, last_gen_2, dir_2, stem_2, print_errors = F)$ped.tmp 
    
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

    all_exchanges <- as.character(sort(unique(as.numeric(c(sent, received)))))
    events <- c()
    
    for (event in all_exchanges) {
        if ((event %in% sent) & !(event %in% received)) {
            event_type <- 'sent'
            if (as.numeric(event) == min(as.numeric(all_exchanges))){
              event_type <- 'sent founders'
            }
        }
        else if (!(event %in% sent) & (event %in% received)) {
            event_type <- 'received'
            if (as.numeric(event) == min(as.numeric(all_exchanges))){
              event_type <- 'received founders'
            }
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
            cat('Cannot find a pedigree file for generation', i, 'at', file.name, '\n')
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


# merge the pedigrees of two exchanging populations into a single pedigree
# with a new set of IDs for all individuals
merge.pedigrees <- function(
    ped_map,     # pedigree map: path to csv with per-population gen numbers for all shared generations
    ex_1_2,      # exchange history from pop1 to pop2: path to csv w/ cols pop1_from, pop2_to
    ex_2_1,      # exchange history from pop2 to pop1: path csv w/ cols pop2_from, pop1_to
    dir_1,       # path to pop1 pedigrees
    stem_1,      # filename stem for pop1 pedigrees
    first_gen_1, # number of the first generation to include from the pop1 pedigree
    last_gen_1,  # number of the final generation to include from the pop1 pedigree
    dir_2,       # path to pop2 pedigrees
    stem_2,      # filename stem for pop2 pedigrees
    first_gen_2, # number of the first generation to include from the pop2 pedigree
    last_gen_2,  # number of the final generation to include from the pop2 pedigree
    merge_into,  # set the merge direction: 1 = from pop2 into pop1, 2 = from pop1 into pop2
    as_df=TRUE,  # set the output type: T = dataframe, F = list of per-gen dataframes
    out_dir,     # the path and base file stem for all output pedigree files
    out_stem)    # the stem name for all merged pedigree files to be saved
{

    # read in generation data and classify pedigree generations for each population
    ped_map <- read.csv(ped_map)
    ex_1_2 <- read.csv(ex_1_2)
    ex_2_1 <- read.csv(ex_2_1)

    pop1_all_gens <- as.character(ex_1_2[,1])
    pop2_all_gens <- as.character(ex_2_1[,1])
    pop1_shared_gens <- as.character(ped_map[,1])
    pop2_shared_gens <- as.character(ped_map[,2]) 
    pop1_separate_gens <- setdiff(pop1_all_gens, pop1_shared_gens)
    pop2_separate_gens <- setdiff(pop2_all_gens, pop2_shared_gens)
    
    pop1_sent <- ex_1_2[complete.cases(ex_1_2),]
    pop2_sent <- ex_2_1[complete.cases(ex_2_1),]
    pop1_received <- as.character(pop2_sent[,2])
    pop2_received <- as.character(pop1_sent[,2])
    pop1_sent <- as.character(pop1_sent[,1])
    pop2_sent <- as.character(pop2_sent[,1])

    # construct named vectors of exchanges
    pop1_exchanges <- name.exchanges(pop1_sent, pop1_received)
    pop2_exchanges <- name.exchanges(pop2_sent, pop2_received)
    
    # format pedigrees for each population
    pop1_ped <- format.pedigree.for.merge(first_gen_1, last_gen_1, dir_1, stem_1)
    pop2_ped <- format.pedigree.for.merge(first_gen_2, last_gen_2, dir_2, stem_2)

    # set up generations to merge: any that recieved samples from the other population,
    # plus the generation prior to each receiving generation
    # merge by 'perspective' of the receiving population: merge the sending pop into the receiving pop
    if (merge_into == 1) {
        receiving_pop <- 'pop1'
        sending_pop <- 'pop2'
        receiving_ped <- pop1_ped
        sending_ped <- pop2_ped
        receiving_exchanges <- pop1_exchanges
        sending_exchanges <- pop2_exchanges
        receiving_all_gens <- pop1_all_gens    
        receiving_history <- ex_2_1 
    } else if (merge_into == 2) {
        receiving_pop <- 'pop2'
        sending_pop <- 'pop1'
        receiving_ped <- pop2_ped
        sending_ped <- pop1_ped
        receiving_exchanges <- pop2_exchanges
        sending_exchanges <- pop1_exchanges
        receiving_all_gens <- pop2_all_gens    
        receiving_history <- ex_1_2
    }
    
    received_gens <- receiving_exchanges[names(receiving_exchanges) %in% c('received','exchanged')]
    received_prior_gens <- as.character(as.numeric(received_gens) - 1)
    receiving_gens_to_merge <- sort(c(received_gens, received_prior_gens))
    
    sent_gens <- as.character(receiving_history[receiving_history[,2] %in% received_gens, 1])
    sent_prior_gens <- as.character(as.numeric(sent_gens) - 1)

    # incorporate simple generations (that don't need merging) into the merged pedigree
    merged_ped <- list()
    for (gen in receiving_all_gens) {
        if (!gen %in% receiving_gens_to_merge){
            ped <- receiving_ped[[gen]]
            ped$alt_generation <- NA
            merged_ped[[gen]] <- ped
        } 
    }
    
    # merge sent/received generations and incorporate into the pedigree
    i <- 0 # generation counter
    for (gen in received_gens) {
        i <- i + 1
        sent_gen <- sent_gens[i]
        received_ped <- receiving_ped[[gen]]
        received_ped$alt_generation <- NA
        sent_ped <- sending_ped[[sent_gen]]
        sent_ped$alt_generation <- sent_gen
        sent_ped$generation <- gen
        ped <- rbind(received_ped, sent_ped)
        ped <- ped[!duplicated(ped[,1]),]
        merged_ped[[gen]] <- ped
    }
    
    # merge generations prior to shipments/receipts and incorporate into the pedigree
    if (length(received_gens) != length(sent_gens)) {
        print(paste('Error: number of generations sent from', sending_pop, '(', length(sent_gens), ') \\
            is not the same as the number of generations received by', receiving_pop, '). Check your \\
            pedigree map and shipment history input files.'))
    }

    i <- 0
    for (gen in received_prior_gens) {
        i <- i + 1
        receipt_ped <- receiving_ped[[gen]]
        receipt_ped$alt_generation <- NA
        if (merge_into == 2) {
            sending_founder_gen <- sending_exchanges[names(sending_exchanges)=='sent founders']     
        } else if (merge_into == 1) {
            sending_founder_gen <- sending_exchanges[names(sending_exchanges)=='received founders']     
        }
       
        sending_first_gen_post_founding <- as.numeric(sending_founder_gen) + 1

        # sending gens: all generations between the last shipment until the current (ith) shipment
        # (all gens following the previous shipment through and including the ith 'prior' generation)
        sent_gen <- as.numeric(sent_gens[i])
        sent_prior_gen <- as.numeric(sent_prior_gens[i])
    
        if (i == 1) {
            sending_gens <- as.character(as.numeric(sending_founder_gen) : as.numeric(sent_prior_gen))
        } else if (i > 1) {
            previous_sent_gen <- sent_gens[i-1]
            sending_gens <- as.character(as.numeric(previous_sent_gen) : as.numeric(sent_prior_gen))
        }
    
        # list of pedigree generations to merge into the 'prior' generation
        sending_peds <- lapply(sending_gens, function(gen) sending_ped[[gen]])
        names(sending_peds) <- sending_gens
        sending_combined_ped <- do.call(rbind, sending_peds)
        sending_combined_ped$alt_generation <- sending_combined_ped$generation
        sending_combined_ped$generation <- gen
    
        # combine all previous generations from the sending population
        # into the receiving 'prior' generation
        ped <- rbind(receipt_ped, sending_combined_ped)
        ped <- ped[!duplicated(ped[,1]),]
        merged_ped[[gen]] <- ped

    } # end of received_prior_gens
    
    # reorder generations in the complete pedigree
    sorted_names <- sort(as.numeric(names(merged_ped)))
    merged_ped <- merged_ped[as.character(sorted_names)]
                               
    # ensure that all generations are present
    check <- sum(names(merged_ped) != receiving_all_gens)
                               
    if (check > 0) {
        cat('Error: Not all', receiving_pop, 'generations were included in the merged pedigree! \n')
    } else {

        # order each generation by ID, then assign new IDs to enable kinship matrix construction
        # by breedail in the correct order
        for (ped in merged_ped) {
            ped$true_id <- ped$id
            ped$sire_id <- ped$sire
            ped$dam_id <- ped$dam
            gen <- as.character(ped$generation[1])
            ped <- ped[order(as.numeric(ped$id)),]
            
            for (i in 1:nrow(ped)){
                
                i_chr <- nchar(i)
                n_zeros <- 4 - i_chr
                zero_str <- paste0(rep('0', n_zeros),collapse='')
                new_id <- as.numeric(paste0(gen, zero_str, i))
                ped$id[i] <- new_id
                id <- ped$id[i]
                rownames(ped) <- NULL
            }
            merged_ped[[gen]] <- ped
        } # end of ped (generation) loop
        
        # change all instances of a new ID throughout the pedigree
        all_ids <- unlist(lapply(merged_ped, function(df) df$id))
        all_true_ids <- unlist(lapply(merged_ped, function(df) df$true_id))
        names(all_ids) <- NULL
        names(all_true_ids) <- NULL
        dups <- duplicated(all_true_ids)
        all_ids <- all_ids[!dups]
        all_true_ids <- all_true_ids[!dups]
        names(all_true_ids) <- all_ids

        for (p in 1:length(merged_ped)) {
            ped <- merged_ped[[p]]
            ped$sire <- names(all_true_ids[match(ped$sire, all_true_ids)])
            ped$dam <- names(all_true_ids[match(ped$dam, all_true_ids)])
            if (p == 1) {
                ped$sire <- 0
                ped$dam <- 0
            }
            merged_ped[[p]] <- ped
        }

        # write pedigree generations to file
        if (dir.exists(out_dir) == FALSE) {
            dir.create(out_dir, showWarnings = TRUE)
        }

        min_gen <- min(as.numeric(receiving_all_gens))
        max_gen <- max(as.numeric(receiving_all_gens))
        max_nchar <- nchar(max_gen)
        for (gen in receiving_all_gens) {
            ped <- merged_ped[[gen]]
            gen_nchar <- nchar(gen)
            n_zeros <- max_nchar - gen_nchar
            zero_str <- paste0(rep('0', n_zeros),collapse='')
            gen_str <- paste0(zero_str, gen)
            file.name <- file.path(out_dir, paste0(out_stem, gen_str, '.csv'))
            write.csv(ped, file.name, row.names=F, quote=F, na='?')
        }
                                    
        ped_df <- do.call(rbind, merged_ped)
        ped_df$sex[ped_df$sex == 0] <- 'F'  # Female = 0
        ped_df$sex[ped_df$sex == 1] <- 'M'  # Male = 1
        outstr <- paste0(out_stem, '_0', min_gen, '_', max_gen, '_complete_ped.csv')
        outfile <- file.path(out_dir, outstr)
        write.csv(ped_df, outfile, row.names=F, quote=F, na='?')
        cat('Merged pedigree written to', outfile, '\n')
        cat('Individual merged pedigree files saved in', paste0(out_dir,'/'), '\n')
          
        # write an ID map for the merged pedigree
        id_map <- map.merged.ids(pop = merge_into)

        if (as_df) {
            merged_ped <- do.call(rbind, merged_ped)  
            rownames(merged_ped) <- NULL
        }
        return(list(pedigree = merged_ped, id.map = id_map, file = outfile))
    }
}


# write an ID map for a merged pedigree
map.merged.ids <- function(pop){
    if(pop==1) {
        id_map <- map_merged_ids_wfu(
            merged_ped = outfile,
            merged_stem = out_stem,
            dir_1 = dir_1,
            stem_1 = stem_1,
            first_gen_1 = first_gen_1,
            last_gen_1 = last_gen_1,
            dir_2 = dir_2,
            stem_2 = stem_2,
            first_gen_2 = first_gen_2,
            last_gen_2 = last_gen_2,
            out_dir = out_dir)
    } else if (pop==2) {
        id_map <- map_merged_ids_hsw(
            merged_ped = outfile,
            merged_stem = out_stem,
            dir_1 = dir_1,
            stem_1 = stem_1,
            first_gen_1 = first_gen_1,
            last_gen_1 = last_gen_1,
            dir_2 = dir_2,
            stem_2 = stem_2,
            first_gen_2 = first_gen_2,
            last_gen_2 = last_gen_2,
            out_dir = out_dir)
    }
    return(id_map)
}


# write multiple pedigree files to one complete pedigree file
write.complete.ped <- function(
    first_gen,
    last_gen,
    data_dir,
    file_stem,
    save_file = TRUE)
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
            cat('Cannot find a pedigree file for generation', i, 'at', file.name, '\n')
        }
        
        ## read in the pedigree for a given generation
        ped.tmp = read.table(file=file.name, sep=",", header=T,
                             as.is=TRUE, na.strings="?")
        # order by ID
        ped.tmp <- as.data.frame(ped.tmp)
        ped.tmp <- ped.tmp[order(ped.tmp[,1]),]
        
        ped <- rbind(ped, ped.tmp)
    }
    if (save_file) {
      outfile <- file.path(data_dir, paste0(file_stem, 
                          '_0', first_gen, '_', last_gen, '_complete_ped.csv'))
      write.csv(ped, outfile, row.names=F, quote=F, na='?')
      cat('Complete pedigree saved to', outfile, '\n')
    }
    return(ped)
}

# translate desired columns from a given ID type to another ID type
translate.merged.ids <- function(
    input,     # the file or list (output by select.breeders()) with IDs to convert
    id_map, # the ID map file to use for conversion
    cols,   # vector of column names from input that need translating
    from,   # vector of column ID types that need changing (a column name from id_map)
    to)     # vector of replacement ID type (a column name from id_map)
{
    if (class(input) == 'character') {
      df <- read.csv(input)
      basefile <- file_path_sans_ext(input)
    } else if (sum(!names(input) %in% c('pairs','file')) == 0) {
      df <- input$pairs
      basefile <- file_path_sans_ext(input$file)
    } else if (sum(!names(input) %in% c('pedigree','id.map','file')) == 0) {
      df <- input$pedigree
      basefile <- file_path_sans_ext(input$file)
    } 
    if (class(id_map) == 'character') {
        id_map <- read.csv(id_map)
    } else if (class(id_map) == 'data.frame') {
      id_map <- id_map
    } else if (sum(!names(input) %in% c('pedigree','id.map')) == 0) {
      id_map <- id_map$id.map
    } 
    
    for (i in 1:length(cols)) {
        col <- cols[i]
        from_id <- from[i]
        to_id <- to[i]
        use_map <- as.character(id_map[[to_id]])
        names(use_map) <- as.character(id_map[[from_id]])
        
        for (r in 1:nrow(df)) {
            id_to_change <- as.character(df[[col]][r])
            df[[col]][r] <- use_map[id_to_change]
        }
    }
    
    from_str <- paste0(unique(from), collapse='_')
    to_str <- paste0(unique(to), collapse='_')
    datestamp <- format(Sys.time(),'%Y%m%d')
    outfile <- paste0(basefile, '_translated_', from_str, '_to_', to_str, '_', datestamp, '.csv')
    write.csv(df, outfile, row.names=F, quote=F, na='')
    cat('Translated file saved to', outfile, '\n')

    return(list(pairs = df, file = outfile))
}


# get the kinship of each potential breeder pair in the current generation
current.kinship <- function(df, # colony df with ALL HSW rats
                           first_gen,
                           last_gen,
                           data_dir,
                           file_stem,
                           id_map = NULL)  # if using a merged pedigree
{
    df <- read.csv(df)
    if (as.numeric(df$generation[1]) != as.numeric(last_gen)){
        stop(paste0('Colony dataframe generation (', df$generation[1], 
                   ') and last_gen (', last_gen, ') must be identical'))
    }
    # format the pedigree for kinship estimation
    ped_for_kinship <- format.pedigree(
        first_gen = first_gen,
        last_gen = last_gen,
        data_dir = data_dir,
        file_stem = file_stem,
        print_errors = F)

    use_ped <- ped_for_kinship$ped
    
    # estimate kinship across the pedigree
    k_all <- kinship(use_ped)
    
    # subset the kinship matrix to only the final generation
    if (is.null(id_map)) {
        kinship_ids <- df$animalid
        kinship_ids <- sapply(kinship_ids, animalid_to_accessid)
        k_use <- k_all[kinship_ids, kinship_ids]
    } else {
        if (class(id_map) == 'character'){
            id_map <- read.csv(id_map)    
        } else if (class(id_map) == 'data.frame'){
          id_map <- id_map
        }
        kinship_map <- id_map[id_map$generation==last_gen,]
        kinship_ids <- as.character(kinship_map$merged_id)
        k_use <- k_all[kinship_ids, kinship_ids]
        for (i in 1:nrow(k_use)) {
            merged_id <- rownames(k_use)[i] 
            animalid <- id_map[id_map$merged_id==merged_id,]$animalid
            rownames(k_use)[i] <- animalid
            colnames(k_use)[i] <- animalid
        }
    }

    # create a file for all hypothetical pairings
    df <- df[df$animalid %in% rownames(k_use),]
    all_males <- df[df$sex=='M',]$animalid
    all_females <- df[df$sex=='F',]$animalid
    n_males <- length(all_males)
    n_females <- length(all_females)

    f_pairs <- rep(all_females, n_males)
    m_pairs <- rep(all_males, each = n_females)

    all_pairs <- data.frame(dam = f_pairs, sire = m_pairs)
    all_pairs$dam_fam <- sapply(all_pairs$dam, function(x) sub('.*_B(\\d+)_.*', '\\1', x))
    all_pairs$sire_fam <- sapply(all_pairs$sire, function(x) sub('.*_B(\\d+)_.*', '\\1', x))

    for (i in 1:nrow(all_pairs)) {
        dam <- all_pairs$dam[i]
        sire <- all_pairs$sire[i]
        all_pairs$kinship[i] <- k_use[dam, sire]
    }
    all_pairs <- all_pairs[order(all_pairs$dam, all_pairs$sire),]
    
    # write kinship and pairing files
    kstr <- paste0(file_stem, '_gen', last_gen, '_kinship_matrix.csv')
    kfile <- file.path(data_dir, kstr)
    write.csv(k_use, kfile, row.names=T, quote=F, na='')
    cat('Kinship matrix saved to', kfile, '\n')
    pairstr <- paste0(file_stem, '_gen', last_gen, '_kinship_all_pairings.csv')
    pairfile <- file.path(data_dir, pairstr)
    write.csv(all_pairs, pairfile, row.names=F, quote=F, na='')
    cat('Kinship pairings saved to', pairfile, '\n')

    return(list(k = k_use, pairs = all_pairs))
}

# print time-stamped messages to stdout
printout <- function(str) {
    cat(paste0('\n[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ']'), 
    str, '\n\n')
}
