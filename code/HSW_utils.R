## PEDIGREE-RELATED FUNCTIONS SPECIFIC TO HS WEST
## written by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

## These are used in conjunction with utils.R, which is more general
## These functions are designed to adhere specifically to data and organizational
## conventions used in-house at the HS West colony. They are likely not useful
## to outside users, though utils.R is of general utility once pedigrees have already
## been properly formatted for use in breedail

source('utils.R')

# convert HSW animal IDs to HSW access IDS
animalid_to_accessid <- function(id){
    
    # remove letters and underscores
    clean_id <- gsub('[A-Za-z_]', '', id)
    
    # append '11' to the start of the new ID
    accessid <- paste0('11', clean_id)
    
    return(accessid)
}

# convert HSW access IDs to HSW animal IDS
accessid_to_animalid <- function(id){
    
    # remove letters and underscores
    id <- as.character(id)
    clean_id <- unlist(strsplit(id,'11'))[2]
    
    rat_no <- substr(clean_id, start = nchar(clean_id)-1, stop = nchar(clean_id))
    pair_no <- substr(clean_id, start = nchar(clean_id)-3, stop = nchar(clean_id)-2)
    # regular expression pattern to match both rat and pair numbers, replace with an empty string
    gen_no <- gsub(paste0(rat_no, '|', pair_no), '', clean_id)
    
    # construct the animal ID
    animalid <- paste0('G', gen_no, '_', 'B', pair_no, '_', rat_no)
    
    return(animalid)
}


# split a raw WFU pedigree into single-generation csv files, still in raw format
split_wfu_raw_ped <- function(ped, outdir) {
    
    # suppress warnings temporarily - excel formatting can produce a lot
    oldw <- getOption("warn")
    options(warn = -1)

    # read in the pedigree excel file
    wfu <- as.data.frame(read_excel(ped))

    # allow warnings again
    options(warn = oldw)

    # format generation
    wfu$Generation <- gsub('00$','',wfu$Generation)
    
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

    # split the pedigree by generation
    ped_gens <- split(wfu, wfu$Generation)
        
    # save generations to separate files
    for (i in 1:length(ped_gens)){
        
        ped <- ped_gens[[i]]
        gen <- names(ped_gens)[[i]]
        if (nchar(gen) == 1) {
            gen <- paste0('0', gen)
        }
        
        write.csv(ped, file.path(outdir, paste0('wfu_raw_gen', gen, '.csv')),
                  row.names=F, quote=F, na='')
    }
}


# function to format a raw complete pedigree from WFU for use in breedail
format_wfu_raw_ped <- function(ped,    # path to raw pedigree excel file
                               outdir) # desired output directory path
{
    # suppress warnings temporarily - excel formatting can produce a lot
    oldw <- getOption("warn")
    options(warn = -1)

    # read in the pedigree excel file
    wfu_in <- as.data.frame(read_excel(ped))

    # allow warnings again
    options(warn = oldw)

    # subset only the rows needed to identify breeders using breedail
    keep_cols <- c('ID.F51', 'Dam.ID', 'Sire.ID', 'Sex', 'Generation','Transpondernumber',
               'SW.ID','Dam.SW.ID','Sire.SW.ID')
    wfu <- wfu_in[,keep_cols]

    # rename columns for compatibility with breedail
    new_colnames <- c('id', 'dam', 'sire', 'sex', 'generation',
                      'rfid','animalid','dam_animalid','sire_animalid')
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
    for (i in 1:length(ped_gens)){
        
        ped <- ped_gens[[i]]
        ped <- ped[order(ped$id),]
        gen <- names(ped_gens)[[i]]
        if (nchar(gen) == 1) {
            gen <- paste0('0', gen)
        }
        
        write.csv(ped, file.path(outdir, paste0('wfu_gen', gen, '.csv')),
                  row.names=F, quote=F, na='?')
    }

}


# function to format a raw WFU pedigree into HSW gens 0-41 (54-95) for breedail
wfu_raw_to_hsw <- function(ped, outdir) {
    
    # suppress warnings temporarily - excel formatting can produce a lot
    oldw <- getOption("warn")
    options(warn = -1)

    # read in the pedigree excel file
    wfu_in <- as.data.frame(read_excel(ped))

    # allow warnings again
    options(warn = oldw)

    # subset only the rows needed to identify breeders using breedail
    keep_cols <- c('ID.F51', 'Dam.ID', 'Sire.ID', 'Sex', 'Generation','Transpondernumber',
               'SW.ID','Dam.SW.ID','Sire.SW.ID')
    wfu <- wfu_in[,keep_cols]

    # rename columns for compatibility with breedail
    new_colnames <- c('id', 'dam', 'sire', 'sex', 'generation',
                      'rfid','animalid','dam_animalid','sire_animalid')
    colnames(wfu) <- new_colnames

    # format NA's for compatibility with breedail
    for (col in 1:ncol(wfu)){
        wfu[,col][which(wfu[,col]=='NA')] <- '?'
        wfu[,col][which(is.na(wfu[,col]))] <- '?'
    }
    
    # remove extra 0's from generation numbers
    wfu$generation <- as.character(wfu$generation)
    wfu$generation <- gsub('00$','',wfu$generation)

    # drop all generations after 40
    # (HS West was founded with rats sent from gen042)
    # (gen041 ped is only the parents of those rats sent in gen042)
    wfu$generation <- as.numeric(wfu$generation)
    hsw <- wfu[wfu$generation < 41,]

    # add 54 to all generation numbers to equal HSW numbering convention
    # (HSW generation counts are in TOTAL HS rat generations, not only WFU generations)
    hsw$wfu_generation <- as.numeric(hsw$generation)
    hsw$generation <- hsw$wfu_generation + 54

    # reorder columns
    new_col_order <- c('id', 'dam', 'sire', 'sex', 'generation', 'wfu_generation',
                      'rfid','animalid','dam_animalid','sire_animalid')
    hsw <- hsw[,new_col_order]

    # split the pedigree by generation
    ped_gens <- split(hsw, hsw$generation)
    print(str(ped_gens))
    # save generations to separate files
    for (i in 1:length(ped_gens)){
        
        ped <- ped_gens[[i]]
        ped <- ped[order(ped$id),]
        gen <- names(ped_gens)[[i]]
        if (nchar(gen) == 1) {
            gen <- paste0('00', gen)
        }
        if (nchar(gen) == 2) {
            gen <- paste0('0', gen)
        }
        
        write.csv(ped, file.path(outdir, paste0('hsw_gen', gen, '.csv')),
                  row.names=F, quote=F, na='?')
    }

}


# function to produce a breedail pedigree from breeders assignments
# df is the path to either a complete assignment file (for all colony rats) or a breeders-specific assignment file
# wfu = path to WFU shipping sheet (xlsx), wfu_gen = WFU generation number provided
hsw_df_to_ped <- function(df, return_df=FALSE, wfu=NULL, wfu_sheet=NULL, wfu_gen=NULL, outdir=NULL) {

    library(readxl)
    df <- read.csv(df)
    hsw_gen <- df$generation[1]
    
    # subset the HSW df to only breeders, depending on df format
    if ('assignment' %in% colnames(df)) {
        df <- df[df$assignment == 'hsw_breeders',]
    } else if ('breeder' %in% colnames(df)) {
        df <- df[df$breeder == 1,]
    } else if ('hsw_breeders' %in% colnames(df)) {
        df <- df[df$hsw_breeders == 1,]
    } else {
        print('Cannot identify an assignment or breeders column to subset')
    }

    # format the HSW df for breedail
    hsw_keep_cols <- c('rfid','generation','animalid','accessid','sex','dam','sire')
    df <- df[,hsw_keep_cols]
    colnames(df) <- c('rfid','generation','animalid','id','sex','dam_animalid','sire_animalid')
    df$wfu_generation <- NA
    df$dam <- sapply(df$dam_animalid, animalid_to_accessid)
    df$sire <- sapply(df$sire_animalid, animalid_to_accessid)
    col_order <- c('id','dam','sire','sex','generation','wfu_generation','rfid',
                   'animalid','dam_animalid','sire_animalid')
    df <- df[,col_order]

    # format the WFU shipping sheet if provided
    if (!is.null(wfu)){
        if (is.null(wfu_gen)){
            print('Please provide a WFU generation number')
            stop()
        }
        if (is.null(wfu_sheet)){
            print('Please indicate which excel sheet to read from the WFU shipping sheet')
            stop()
        }
        
        wfu <- as.data.frame(read_excel(wfu, sheet=wfu_sheet))
        wfu_keep_cols <- c('Dam','Sire','Access ID','Animal ID','Transponder ID','Sex')
        wfu <- wfu[,wfu_keep_cols]
        colnames(wfu) <- c('dam','sire','id','animalid','rfid','sex')
        wfu <- wfu[!is.na(wfu$id),]
        wfu$wfu_generation <- wfu_gen
        wfu$dam <- gsub('_', '', wfu$dam)
        wfu$sire <- gsub('_', '', wfu$sire)
        wfu$id <- gsub('_', '', wfu$id)
        wfu$dam_animalid <- NA
        wfu$sire_animalid <- NA
        wfu$generation <- hsw_gen
        wfu <- wfu[,col_order]

        # concatenate HSW and WFU rats
        df <- rbind(df, wfu)
    }

    df <- df[order(as.numeric(df$id)),]

    if (!is.null(outdir)){
        
        if (nchar(hsw_gen) == 1) {
            gen <- paste0('00', hsw_gen)
        } else if (nchar(hsw_gen) == 2) {
            gen <- paste0('0', hsw_gen)
        } else if (nchar(hsw_gen) == 3) {
            gen <- hsw_gen
        }
        
        write.csv(df, file.path(outdir, paste0('hsw_gen', gen, '.csv')), row.names=F, quote=F, na='?')
    }

    if (!is.logical(return_df)) {
        stop("return_df should be a logical value")
    }

    if (return_df) {
        return(df)    
    } 

}
