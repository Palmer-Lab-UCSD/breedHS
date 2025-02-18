## PEDIGREE-RELATED FUNCTIONS SPECIFIC TO WAKE FOREST
## written by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

## These are used in conjunction with utils.R, which is more general
## These functions are designed to adhere specifically to data and organizational
## conventions used in-house at the WFU colony. They are likely not useful
## to outside users, though utils.R is of general utility once pedigrees have already
## been properly formatted for use in breedail

library(readxl)


# convert WFU SW.ID to access ID
swid_to_accessid <- function(id, wfu_map) {

    accessid <- wfu_map[wfu_map[['swid']] == id,][['accessid']]
    return(accessid)
}

# convert access ID to WFU SW.ID
accessid_to_swid <- function(id, wfu_map) {

    swid <- wfu_map[wfu_map[['accessid']] == id,][['swid']]
    return(swid)
}

# find a WFU animal ID (SW.ID) given a WFU access ID
get_wfu_swid <- function(id, wfu_df){
    wfu_df <- read.csv(wfu_df, na.str=c('','NA','NaN','nan'))
    sw_id <- wfu_df[wfu_df[,1] == id,]$SW.ID
    return(sw_id)
}


get_wfu_accessid <- function(
    id, 
   shipping_sheet) # formatted shipping sheet, read into R
{
    wfu <- shipping[wfu$animalid == id,]
    accessid <- wfu$id[1]
    accessid <- gsub('_', '', accessid)
    return(accessid)
}

# function to read in an unformatted WFU pedigree file, return an R dataframe
read_wfu_raw_ped <- function(
    ped) # path to csv or xlsx file
{

    # read in the pedigree file
    if (file_ext(ped) == 'xlsx') {
        # suppress warnings temporarily - excel formatting can produce a lot
        oldw <- getOption('warn')
        options(warn = -1)
        wfu <- as.data.frame(read_excel(ped))
        options(warn = oldw) # allow warnings again
    } else if (file_ext(ped) == 'csv') {
        wfu <- read.csv(ped, na.str=c('','NA','NaN','nan'))
    }

    # format generation - remove trailing double zeros
    # wfu$Generation <- gsub('00$','',wfu$Generation)
    
    # keep only real data columns
    keep_cols <- c()
    for (col in colnames(wfu)) {
        if (!startsWith(col,'...')){
            keep_cols <- c(keep_cols, col)
        }
    }
    wfu <- wfu[,keep_cols]

    # format NA's
    for (col in 1:ncol(wfu)){
        wfu[,col] <- as.character(wfu[,col])
        wfu[,col][wfu[,col]=='NA'] <- NA
    }

    return(wfu)
}

# split a raw WFU pedigree into single-generation csv files, still in raw format
split_wfu_raw_ped <- function(
    ped,    # the complete pedigree, either xlsx or csv format
    outdir) # the desired directory in which to save per-gen raw pedigree files
{
    
    wfu <- read_wfu_raw_ped(ped)

    # split the pedigree by generation
    ped_gens <- split(wfu, wfu$Generation)
        
    # save generations to separate files
    if (dir.exists(outdir) == FALSE) {
        dir.create(outdir, showWarnings = TRUE)
    }

    for (i in 1:length(ped_gens)){
        
        ped <- ped_gens[[i]]
        gen <- names(ped_gens)[[i]]
        gen <- gsub('00$','', gen)
        if (nchar(gen) == 1) {
            gen <- paste0('0', gen)
        }
        write.csv(ped, file.path(outdir, paste0('wfu_raw_gen', gen, '.csv')),
                  row.names=F, quote=F, na='')
    }
    cat('Split pedigree generations saved to', outdir, '\n')
}


# function to format a raw complete pedigree from WFU for use in breedail
format_wfu_raw_ped <- function(
    ped,    # path to any raw pedigree file (single- or multi-gen, xlsx or csv)
    outdir) # desired output directory path
{

    wfu <- read_wfu_raw_ped(ped)
    
    # subset only the rows needed to identify breeders using breedail
    keep_cols <- c('ID.F51','Dam.ID','Sire.ID','Sex','Generation','Transpondernumber','SW.ID','Dam.SW.ID','Sire.SW.ID')
    wfu_keepcols <- setdiff(keep_cols, colnames(wfu))
    if (length(wfu_keepcols)>0) {
        keep_cols <- c('IDF51','DamID','SireID','Sex','Generation','Transpondernumber','SWID','DamSWID','SireSWID')
        wfu_keepcols <- setdiff(keep_cols, colnames(wfu))
        if (length(wfu_keepcols)>0) {
            cat('Make sure the pedigree has the following columns:', keep_cols, '\n')
            cat('File:', ped)
        }
    }
    wfu <- wfu[,keep_cols]

    # rename columns for compatibility with breedail
    new_colnames <- c('id','dam','sire','sex','generation','rfid','swid','dam_swid','sire_swid')
    colnames(wfu) <- new_colnames

    # format NA's for compatibility with breedail
    for (col in 1:ncol(wfu)){
        wfu[,col][which(wfu[,col]=='NA')] <- '?'
        wfu[,col][which(is.na(wfu[,col]))] <- '?'
    }

    # format generation
    wfu$generation <- gsub('00$','',wfu$generation)
    
    # split the pedigree by generation
    ped_gens <- split(wfu, wfu$generation)
        
    # save generations to separate files
    if (dir.exists(outdir) == FALSE) {
        dir.create(outdir, showWarnings = TRUE)
    }
    
    for (i in 1:length(ped_gens)){
        
        ped <- ped_gens[[i]]
        ped <- ped[order(ped$id),]
        gen <- names(ped_gens)[[i]]
        if (nchar(gen) == 1) {
            gen <- paste0('0', gen)
        }

        outfile <- file.path(outdir, paste0('wfu_gen', gen, '.csv'))
        write.csv(ped, outfile, row.names=F, quote=F, na='?')
    }
}
