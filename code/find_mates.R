find.mates <- function(ped, k, verbose=FALSE){

  ## k is a matrix of kinship coefficients
  ## sex is a matrix of 0/1, where 1 indicates that the
  ## two individuals are of opposite sex and 0 they are of the same sex

  ##k.orig = k
  
  ibc = diag(k)
  
  idx.m = ped[,4]=="M"
  idx.f = ped[,4]=="F"
  if (sum(c(idx.m,idx.f)) ==  dim(ped)[1]){
    ped[idx.m,4] = 1
    ped[idx.f,4] = 0
  }
  ped.orig = ped

  sex.vec = as.numeric(ped[,4])
  ids = as.numeric(ped[,1])
  N = length(ids)

  sex.1 = kronecker(as.numeric(ped[,4]),t(rep(1,N)))
  sex.2 = t(sex.1)
  sex = sex.1 != sex.2

  output = c()
  
  while (sum(sex)>0){
    
    min.k = min(k[sex])
    min.idx = k == min.k & upper.tri(k,diag=FALSE) & sex
    m = kronecker(1:N,t(rep(1,N)))
    id.1 = m[min.idx]
    id.2 = t(m)[min.idx]
    ## if more than one choice, choose those with smallest average inbreeding
    ## coefficient
    if (length(id.1) > 1){
      mean.ibc = apply(rbind(ibc[id.1],ibc[id.2]),2,mean)
      ord = order(mean.ibc)
      id.1 = id.1[ord[1]]
      id.2 = id.2[ord[1]]
    }
    
    ## remove appropriate rows and columns ##
    N = N-2
    sex = sex[-c(id.1,id.2),-c(id.1,id.2)]
    k = k[-c(id.1,id.2),-c(id.1,id.2)]
    ibc = ibc[-c(id.1,id.2)]
    output = rbind(output ,
      c(c(ids[id.1],ids[id.2])[sex.vec[c(id.1,id.2)]+1] , round(min.k,4)) )
    sex.vec = sex.vec[-c(id.1,id.2)]
    ids = ids[-c(id.1,id.2)]
  }
  
  return(output)
}


#################################################
## Same as above but only using one individual ##
## per sibship                                 ##
#################################################
find.mates.res <- function(ped, k, verbose=FALSE){

  ## k is a matrix of kinship coefficients
  ## sex is a matrix of 0/1, where 1 indicates that the
  ## two individuals are of opposite sex and 0 they are of the same sex

  ##k.orig = k
  
  ibc = diag(k)
  
  idx.m = ped[,4]=="M"
  idx.f = ped[,4]=="F"

  if (sum(c(idx.m,idx.f)) ==  dim(ped)[1]){
    ped[idx.m,4] = 1
    ped[idx.f,4] = 0
  }
  ##ped.orig = ped
  ped.tmp = ped
  
  sex.vec = as.numeric(ped[,4])
  ids = as.numeric(ped[,1])
  N = length(ids)

  sex.1 = kronecker(as.numeric(ped[,4]),t(rep(1,N)))
  sex.2 = t(sex.1)
  sex = sex.1 != sex.2

  output = c()
  
  while (sum(sex)>0){
    
    min.k = min(k[sex])
    min.idx = k == min.k & upper.tri(k,diag=FALSE) & sex
    m = kronecker(1:N,t(rep(1,N)))
    id.1 = m[min.idx]
    id.2 = t(m)[min.idx]
    ## if more than one choice, choose those with smallest average inbreeding
    ## coefficient
    if (length(id.1) > 1){
      mean.ibc = apply(rbind(ibc[id.1],ibc[id.2]),2,mean)
      ord = order(mean.ibc)
      id.1 = id.1[ord[1]]
      id.2 = id.2[ord[1]]
    }
    parents.1 = as.numeric(ped.tmp[ped.tmp[,1]==ids[id.1],2:3])
    parents.2 = as.numeric(ped.tmp[ped.tmp[,1]==ids[id.2],2:3])

    ## check if any of parents are shared between breeding pair ##
    problem = sum(parents.1 == parents.2)
    
    if (verbose) {
    print(paste('id.1:', id.1))
    print(paste('id.2:', id.2))
    print('ped.tmp:')
    print(ped.tmp)
    print(paste('parents.1:', parents.1))
    print(paste('parents.1:', parents.2))
    print('problem:')
    print(problem)
    }
    
    if (problem){
      cat("Remaining possible mating pairs share at least one parent.  Exiting. . .\n")
      return(output)
    }
    
    ## find siblings ##
    match1 = ped.tmp[,2] == parents.1[1] & ped.tmp[,3] == parents.1[2]
    match2 = ped.tmp[,2] == parents.2[1] & ped.tmp[,3] == parents.2[2]
    sibs.1 = ped.tmp[match1,1]
    sibs.2 = ped.tmp[match2,1]
    idx.rm = ids %in% c(sibs.1,sibs.2)
    idx.keep = idx.rm==FALSE
    ## remove appropriate rows and columns ##
    N = N-sum(idx.rm)
    sex = sex[idx.keep,idx.keep]
    k = k[idx.keep,idx.keep]
    ped.tmp = ped.tmp[idx.keep,]
    ibc = ibc[idx.keep]
    output = rbind(output ,
      c(c(ids[id.1],ids[id.2])[sex.vec[c(id.1,id.2)]+1] , round(min.k,4)) )
    ## sex.vec = sex.vec[-c(id.1,id.2)]
    sex.vec = sex.vec[idx.keep]
    ids = ids[idx.keep]
  }
  
  return(output)
}




####################################################
## Given a set of mates (or their offspring)      ##
## and a non-overlapping set of individuals from  ##
## which to select additional breeding pairs,     ##
## select the pairs the minimizes the mean        ##
## coefficient of the offspring                   ##
####################################################

find.mates.given.pop <- function(ped , k , parent.pos, children.pos, cant.breed.pos = c(), max.pairs = 50, verbose=FALSE) {

  ## Given children of generation n-1 (generation n), select 
  ## additional parents from generation n-1 to minimize the
  ## relationship between the current set of children and those
  ## from the selected parents

  ## cant.breed.pos = position of the parental population that are no long
  ## available for breeding, but who we would still like to have
  ## represented in the offspring. E.g. Use their kinship coefficients in
  ## calculating mean kinships, but don't select as parents.
  
  ## potential parent pedigree ##
  if (length(c(parent.pos,children.pos)) != dim(ped)[1]){
    cat("Warning: number of parents and children specified doesn't ",
        "match the number of individuals in the pedigree file.",
        
        "Was this intentional?\n")
  }

  ibc = diag(k)
  ped.p = ped[parent.pos,]
  ped.tmp = ped.p
  
  ids = as.numeric(ped.tmp[,1])

  sex.vec = as.numeric(ped.p[,4])
  N = length(ids)
  
  sex.1 = kronecker(as.numeric(ped.tmp[,4]),t(rep(1,N)))
  sex.2 = t(sex.1)
  sex = sex.1 != sex.2
  ## adjust sex so that individuals who can't be breed are always 0 ##
  if (length(cant.breed.pos) > 0){
    can.breed = rep(TRUE,N)
    can.breed[cant.breed.pos] = FALSE
    breed.1 = kronecker(can.breed , t(rep(1,N)))
    breed.2 = t(breed.1)
    breed = breed.1 & breed.2
    ## update sex to to exclude non-breeders ##
    sex = sex & breed
  }
  
  ## calculate mean kinship coef between each parent and the offspring
  k.mean = rep(0,N)
  for (i in 1:N){
    k.mean[i] = sum(k[i,children.pos])
  }
  k.tmp = k[parent.pos,parent.pos]
  k.funky = k.tmp +
    kronecker( k.mean , t(rep(1,N)) ) +
    t(kronecker( k.mean , t(rep(1,N))) )

  ## No longer need kinship coef or ids of the children ##
  output = c()
  count = 0
  while(sum(sex) > 0 & count < max.pairs) {

    # cat(count, ":", output,"\n")
    count = count+1
    
    min.k = min(k.funky[sex])
    min.idx = k.funky == min.k & upper.tri(k.funky,diag=FALSE) & sex
    m = kronecker(1:N,t(rep(1,N)))
    id.1 = m[min.idx]
    id.2 = t(m)[min.idx]

    ## if more than one choice, choose those with smallest average
    ## inbreeding coefficient
    if (length(id.1) > 1){
      mean.ibc = apply(rbind(ibc[id.1],ibc[id.2]),2,mean)
      ord = order(mean.ibc)
      id.1 = id.1[ord[1]]
      id.2 = id.2[ord[1]]
    }
    ## parents.1 = as.numeric(ped.tmp[ped.tmp[,1]==ids[id.1],2:3])
    ## parents.2 = as.numeric(ped.tmp[ped.tmp[,1]==ids[id.2],2:3])
    parents.1 = c(ped.tmp[ped.tmp[,1]==ids[id.1],2:3])
    parents.2 = c(ped.tmp[ped.tmp[,1]==ids[id.2],2:3])

    ## find siblings ##
    match1 = ped.tmp[,2] == parents.1[1] & ped.tmp[,3] == parents.1[2]
    match2 = ped.tmp[,2] == parents.2[1] & ped.tmp[,3] == parents.2[2]
    sibs.1 = ped.tmp[match1,1]
    sibs.2 = ped.tmp[match2,1]
    idx.rm = ids %in% c(sibs.1,sibs.2)
    idx.keep = idx.rm==FALSE

    k.pair = k.tmp[id.1,id.2]
    ## if more breedings are possible ##
    if (sum(idx.keep)>1){
    ## update k.funky matrix ##  
      for (i in 1:(sum(idx.keep)-1)) {
        for (j in (i+1):(sum(idx.keep)) ){
          id.i = c(1:N)[idx.keep][i]
          id.j = c(1:N)[idx.keep][j]
          k.funky[id.i,id.j] = k.funky[id.i,id.j] +
            .5 * ( k.tmp[id.i , id.1] + k.tmp[id.i , id.2] +
                  k.tmp[id.j , id.1] + k.tmp[id.j , id.2])
          k.funky[id.j,id.i] = k.funky[id.i,id.j]
        }
      }

    }
    
         
    ## remove appropriate rows and columns ##
    N = N-sum(idx.rm)
    sex = sex[idx.keep,idx.keep]
    k = k[idx.keep,idx.keep]
    ped.tmp = ped.tmp[idx.keep,]
    k.tmp = k.tmp[idx.keep,idx.keep]
    k.funky = k.funky[idx.keep , idx.keep]
    ibc = ibc[idx.keep]
    output = rbind(output ,
      c(c(ids[id.1],ids[id.2])[sex.vec[c(id.1,id.2)]+1] , round(k.pair,4),round(min.k,4)) )
    ## sex.vec = sex.vec[-c(id.1,id.2)]
    sex.vec = sex.vec[idx.keep]
    ids = ids[idx.keep]
  }
  
  return(output)    
}


################################################
## base pairing on mean kinship coefficients  ##
################################################

find.mates.mk <- function(ped, verbose=FALSE){

  k = k.orig
  ibc = diag(k)
  N = length(ibc)
  
  ped = ped.orig
  sex = as.numeric(ped[,4])
  
  idx.m = sex==1
  idx.f = sex==0

  code.m = c(1:N)[idx.m]
  code.f = c(1:N)[idx.f]
  id.m = ped[idx.m,1]
  id.f = ped[idx.f,1]
  
  mean.k = apply(k, 2, mean) - ibc/N 
  order.m = order(mean.k[idx.m])
  order.f = order(mean.k[idx.f])
  N.pair = min(sum(idx.m),sum(idx.f))

  k.pairs = rep(0,N.pair)
  for (i in 1:N.pair){
    k.pairs[i] = k[code.m[order.m[i]] , code.f[order.f[i]]]
  }

  out = cbind(id.f[order.f[N.pair]] , id.m[order.m[N.pair]],
    mean.k[code.f[order.f[N.pair]]], mean.k[code.m[order.m[N.pair]]] , k.pairs) 
  return(out)
}

######################################################
## CHECK TO INSURE THAT PARENTS ARE IN THE PEDIGREE ##
######################################################

check.ped <- function(ped.current, ped.prev){

  parents = c(ped.current[,2],ped.current[,3])
  parents = unique(parents)
  idx = parents == 0
  parents = parents[idx==F]

  idx = parents %in% ped.prev[,1] == F
  missing.parents = parents[idx]

  return(missing.parents)
}


### modifications to prevent reciprocal crosses (eg if B01_M x B02_F disallow B01_F x B02_M)
UNKNOWN.FAM <- "0_0"

make.fam.key <- function(p1, p2) paste(p1, p2, sep = "_")
fam.pair.key <- function(f1, f2) paste(sort(c(f1, f2)), collapse = "||")

## family key of one individual, looked up in a reference pedigree
fam.of <- function(id, ped.ref) {
  r <- which(ped.ref[, 1] == id)[1]
  make.fam.key(ped.ref[r, 2], ped.ref[r, 3])
}

## family-pair keys represented by a table of pairs (cols 1,2 = the two IDs)
banned.from.pairs <- function(pairs, ped.ref) {
  if (is.null(pairs) || nrow(pairs) == 0) return(character(0))
  keys <- apply(pairs, 1, function(row) {
    f1 <- fam.of(as.numeric(row[1]), ped.ref)
    f2 <- fam.of(as.numeric(row[2]), ped.ref)
    if (f1 == UNKNOWN.FAM || f2 == UNKNOWN.FAM) return(NA_character_)
    fam.pair.key(f1, f2)
  })
  unique(keys[!is.na(keys)])
}

## turn off any opposite-sex entry that mates the same family or an
## already-used family pair
apply.family.mask <- function(sex, fam, banned.pairs) {
  n <- length(fam)
  if (n < 2) return(sex)
  for (a in 1:(n - 1)) {
    fa <- fam[a]; if (fa == UNKNOWN.FAM) next
    for (b in (a + 1):n) {
      if (!sex[a, b]) next
      fb <- fam[b]; if (fb == UNKNOWN.FAM) next
      if (fa == fb || fam.pair.key(fa, fb) %in% banned.pairs) {
        sex[a, b] <- FALSE; sex[b, a] <- FALSE
      }
    }
  }
  sex
}

# this function is derived from find.mates.given.pop() to prevent reciprocal crosses
unique.mates.given.pop <- function(
  ped, 
  k, 
  parent.pos, 
  children.pos,
  cant.breed.pos = c(), 
  max.pairs = 50,
  banned.pairs = character(0), 
  verbose = FALSE) {

  if (length(c(parent.pos, children.pos)) != dim(ped)[1]) {
    cat("Warning: number of parents and children specified doesn't ",
        "match the number of individuals in the pedigree file.",
        "Was this intentional?\n")
  }

  ibc     <- diag(k)
  ped.p   <- ped[parent.pos, ]
  ped.tmp <- ped.p
  ids     <- as.numeric(ped.tmp[, 1])
  sex.vec <- as.numeric(ped.p[, 4])
  N       <- length(ids)

  fam <- make.fam.key(ped.tmp[, 2], ped.tmp[, 3])          # NEW: family key per parent

  sex.1 <- kronecker(as.numeric(ped.tmp[, 4]), t(rep(1, N)))
  sex.2 <- t(sex.1)
  sex   <- sex.1 != sex.2

  if (length(cant.breed.pos) > 0) {
    can.breed <- rep(TRUE, N); can.breed[cant.breed.pos] <- FALSE
    breed.1 <- kronecker(can.breed, t(rep(1, N))); breed.2 <- t(breed.1)
    sex <- sex & (breed.1 & breed.2)
  }

  sex <- apply.family.mask(sex, fam, banned.pairs)         # NEW: forbid same-family & used pairs

  k.mean <- rep(0, N)
  for (i in 1:N) k.mean[i] <- sum(k[i, children.pos])
  k.tmp   <- k[parent.pos, parent.pos]
  k.funky <- k.tmp + kronecker(k.mean, t(rep(1, N))) +
                     t(kronecker(k.mean, t(rep(1, N))))

  output <- c(); count <- 0
  while (sum(sex) > 0 & count < max.pairs) {
    count <- count + 1
    min.k   <- min(k.funky[sex])
    min.idx <- k.funky == min.k & upper.tri(k.funky, diag = FALSE) & sex
    m       <- kronecker(1:N, t(rep(1, N)))
    id.1    <- m[min.idx]; id.2 <- t(m)[min.idx]
    if (length(id.1) > 1) {
      mean.ibc <- apply(rbind(ibc[id.1], ibc[id.2]), 2, mean)
      ord <- order(mean.ibc); id.1 <- id.1[ord[1]]; id.2 <- id.2[ord[1]]
    }

    k.pair <- k.tmp[id.1, id.2]

    if (fam[id.1] != UNKNOWN.FAM && fam[id.2] != UNKNOWN.FAM)   # NEW: ban this combination
      banned.pairs <- union(banned.pairs, fam.pair.key(fam[id.1], fam[id.2]))

    idx.keep <- rep(TRUE, N); idx.keep[c(id.1, id.2)] <- FALSE  # NEW: remove ONLY the two chosen

    if (sum(idx.keep) > 1) {
      keep.pos <- which(idx.keep)
      for (a in 1:(length(keep.pos) - 1)) for (b in (a + 1):length(keep.pos)) {
        id.i <- keep.pos[a]; id.j <- keep.pos[b]
        k.funky[id.i, id.j] <- k.funky[id.i, id.j] +
          .5 * (k.tmp[id.i, id.1] + k.tmp[id.i, id.2] +
                k.tmp[id.j, id.1] + k.tmp[id.j, id.2])
        k.funky[id.j, id.i] <- k.funky[id.i, id.j]
      }
    }

    output <- rbind(output,
      c(c(ids[id.1], ids[id.2])[sex.vec[c(id.1, id.2)] + 1],
        round(k.pair, 4), round(min.k, 4)))

    N       <- sum(idx.keep)
    sex     <- sex[idx.keep, idx.keep, drop = FALSE]         # drop=FALSE guards the 1-left collapse
    k       <- k[idx.keep, idx.keep, drop = FALSE]
    ped.tmp <- ped.tmp[idx.keep, , drop = FALSE]
    k.tmp   <- k.tmp[idx.keep, idx.keep, drop = FALSE]
    k.funky <- k.funky[idx.keep, idx.keep, drop = FALSE]
    ibc     <- ibc[idx.keep]; sex.vec <- sex.vec[idx.keep]
    ids     <- ids[idx.keep]; fam <- fam[idx.keep]

    sex <- apply.family.mask(sex, fam, banned.pairs)        # NEW: re-mask with the new ban
  }
  return(output)
}