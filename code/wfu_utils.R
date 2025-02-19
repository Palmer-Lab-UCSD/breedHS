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
    new_colnames <- c('id','dam','sire','sex','generation','rfid','animalid','dam_animalid','sire_animalid')
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

# format a breeder file to assist with pairing in the HSW colony
## INCOMPLETE
create_wfu_breeder_file <- function(
    pairs,          # R dataframe or path to csv, as output by select.breeders, dam/sire must be SW IDs
    raw_ped,        # the raw WFU pedigree for the current generation
    hsw_ss=NULL,    # HSW shipping sheet, if pairing with shipped HSW rats
    outdir=NULL)
{
    if (class(pairs) == 'data.frame') {
        pairs <- pairs
    } else if (class(pairs) == 'character') {
        pairs <- read.csv(pairs, na.str=c('','NA','NaN','nan'))
    } else if (sum(!names(pairs) %in% c('pairs','file')) == 0) {
        pairs <- pairs$pairs
    }
    if (!is.null(hsw_ss)) {    
        hsw <- read.csv(hsw_ss)
    }
    ped <- read_wfu_raw_ped(raw_ped)
    gen <- as.numeric(ped$Generation[1])
    pairs$generation <- gen
    pairs$breederpair <- c(
        paste0('WFU', gen + 1, '_B0', seq.int(1,9)),
        paste0('WFU', gen + 1, '_B', seq.int(10,nrow(pairs)))
    )

    pairs$sire_swid <- sapply(pairs$sire, function(id) {
        if (id %in% ped$IDF51) {
            sire_id <- ped$SWID[which(ped$IDF51==id)]
        } else if (id %in% ped$SWID) {
            sire_id <- id
        }
        return(sire_id)
    })
    pairs$sire_swid <- sapply(pairs$sire, function(id) {
        if (id %in% ped$IDF51) {
            sire_id <- ped$SWID[which(ped$IDF51==id)]
        } else if (id %in% ped$SWID) {
            sire_id <- id
        }
        return(sire_id)
    })

    pairs$sire_animalid <- pairs$sire
    pairs$dam_animalid <- pairs$dam

    all_cols <- c('generation','breederpair','kinship','dam_rfid','dam_animalid','dam_earpunch',
                       'dam_coatcolor','dam_rack_num','dam_rack_pos','dam_wfu_cage',
                       'sire_rfid','sire_animalid','sire_earpunch','sire_coatcolor',
                       'sire_rack_num','sire_rack_pos','sire_wfu_cage')    
    
    for (col in setdiff(all_cols, names(pairs))) {
        pairs[[col]] <- NA
    }
    for (i in 1:nrow(pairs)) {
    
        dam <- pairs$dam[i]
        
        if (substr(dam,1,1)=='G') {
            dam_df <- df[df$animalid==dam,]
            pairs$dam_rfid[i] <- as.character(dam_df$rfid)
            pairs$dam_earpunch[i] <- dam_df$earpunch
            pairs$dam_coatcolor[i] <- dam_df$coatcolor
            pairs$dam_rack_num[i] <- dam_df$rack_number
            pairs$dam_rack_pos[i] <- dam_df$rack_position            
        } else if (substr(dam,1,1)=='H') {
            dam_df <- wfu[wfu[['Animal ID']]==dam,]
            pairs$dam_rfid[i] <- dam_df[['Transponder ID']]
            pairs$dam_earpunch[i] <- dam_df[['Ear Punch']]
            pairs$dam_coatcolor[i] <- dam_df[['Coat Color']]
            pairs$dam_rack_num[i] <- NA
            pairs$dam_rack_pos[i] <- NA
            pairs$dam_wfu_cage[i] <- dam_df[['Ship Box']]
        }

        sire <- pairs$sire[i]
        
        if (substr(sire,1,1)=='G') {
            sire_df <- df[df$animalid==sire,]
            pairs$sire_rfid[i] <- as.character(sire_df$rfid)
            pairs$sire_earpunch[i] <- sire_df$earpunch
            pairs$sire_coatcolor[i] <- sire_df$coatcolor
            pairs$sire_rack_num[i] <- sire_df$rack_number
            pairs$sire_rack_pos[i] <- sire_df$rack_position            
        } else if (substr(sire,1,1)=='H') {
            sire_df <- wfu[wfu[['Animal ID']]==sire,]
            pairs$sire_rfid[i] <- sire_df[['Transponder ID']]
            pairs$sire_earpunch[i] <- sire_df[['Ear Punch']]
            pairs$sire_coatcolor[i] <- sire_df[['Coat Color']]
            pairs$sire_rack_num[i] <- NA
            pairs$sire_rack_pos[i] <- NA
            pairs$sire_wfu_cage[i] <- sire_df[['Ship Box']]
        }
    }

    if (is.null(wfu_ss)) {
        col_order <- c('generation','breederpair','kinship','dam_rfid','dam_animalid','dam_earpunch',
                       'dam_coatcolor','dam_rack_num','dam_rack_pos','sire_rfid',
                       'sire_animalid','sire_earpunch','sire_coatcolor','sire_rack_num','sire_rack_pos')

    } else {
        col_order <- c('generation','breederpair','kinship','dam_rfid','dam_animalid','dam_earpunch',
                       'dam_coatcolor','dam_rack_num','dam_rack_pos','dam_wfu_cage',
                       'sire_rfid','sire_animalid','sire_earpunch','sire_coatcolor',
                       'sire_rack_num','sire_rack_pos','sire_wfu_cage')

    }
    
    # final formatting
    pairs <- pairs[,col_order]
    pairs$kinship <-format(round(pairs$kinship, 4), nsmall=4)
    pairs$kinship <- sapply(pairs$kinship, function(x) paste0(x, paste0(rep(0,6-nchar(x)), collapse='')))

    # final columns: two to fill/log as physical pairing happens, one with pairing notes
    pairs$paired_dam <- 'NONE'
    pairs$paired_sire <- 'NONE'
    pairs$comments <- 'breederpair assigned using breedHS'
    
    if (!is.null(outdir)) {
        datestamp <- format(Sys.time(),'%Y%m%d')
        outfile <- paste0('hsw_gen', gen, '_', gen+1, '_breeders_proposed_', datestamp, '.csv')
        outfile <- file.path(outdir, outfile)
        write.csv(pairs, outfile, row.names=F, quote=F, na='')
        cat('Breeder file written to', outfile, '\n')
    }
    return(list(pairs = pairs, file = outfile))
}


# incorporate HSW IDs into the WFU pedigree prior to an animal transfer
hsw_into_wfu_raw <- function(
    wfu_raw,   # raw WFU pedigree for the current generation
    hsw_assignments,    # a current HSW assignment sheet (csv) with IDs to add to the current WFU generation
    outdir=NULL,    # directory in which to save the formatted pedigree file
    return_df=FALSE) # whether to return the final df to the R console
{
    wfu <- read.csv(wfu_raw, na.str=c('','NA','NaN','nan'))
    hsw <- read.csv(hsw_assignments, na.str=c('','NA','NaN','nan'))

    wfu_gen <- wfu$Generation[1]
    hsw <- hsw[hsw$assignment=='hsw_breeders' & hsw$sex=='M',]

    # format HSW breeders to append to the WFU pedigree
    hsw$sire_animalid <- hsw$sire
    hsw$dam_animalid <- hsw$dam
    hsw$sire_accessid <- sapply(hsw$sire_animalid, animalid_to_accessid)
    hsw$dam_accessid <- sapply(hsw$dam_animalid, animalid_to_accessid)
    hsw$homecage <- NA
    hsw$generation <- wfu_gen

    hsw_col_order <- c('accessid','coatcolor','sex','breederpair','dob','homecage','sire_accessid','dam_accessid',
    'rfid','animalid','generation','earpunch','sire_animalid','dam_animalid')
    wfu_colnames <- c('IDF51','CC','Sex','FamNo','DOB','HomeCage','SireID','DamID',
    'Transpondernumber','SWID','Generation','EarPunch','SireSWID','DamSWID')
    hsw <- hsw[,hsw_col_order]
    colnames(hsw) <- wfu_colnames
    
    # concatenate HSW and WFU rats
    out <- rbind(wfu, hsw)

    if (!is.null(outdir)){
        
        wfu_gen <- gsub('00$','', wfu$Generation[1])
        if (nchar(wfu_gen) == 1) {
            gen <- paste0('0', wfu_gen)
        } else if (nchar(wfu_gen) == 2) {
            gen <- wfu_gen
        }
        
        if (dir.exists(outdir) == FALSE) {
            dir.create(outdir, showWarnings = TRUE)
        }

        outfile <- file.path(outdir, paste0('wfu_plus_hsw_raw_gen', gen, '.csv'))
        write.csv(out, wfu_raw, row.names=F, quote=F, na='')
        cat('HSW IDs incorporated into WFU pedigree file', wfu_raw, '\n')
    }

    if (!is.logical(return_df)) {
        stop("return_df should be a logical value")
    }

    if (return_df) {
        return(out)    
    } 
    
}


# create a map for HSW IDs used in a merged pedigree
map_merged_ids_wfu <- function(
    merged_ped, # path to a complete merged pedigree
    merged_stem,
    dir_1,
    stem_1,
    first_gen_1,
    last_gen_1,
    dir_2,
    stem_2,
    first_gen_2,
    last_gen_2,
    out_dir=NULL)
{

    merged_ped <- read.csv(merged_ped)
    
    ped1 <- write.complete.ped(
        first_gen = first_gen_1,
        last_gen = last_gen_1,
        data_dir = dir_1,
        file_stem = stem_1,
        save_file = FALSE)

    ped2 <- write.complete.ped(
        first_gen = first_gen_2,
        last_gen = last_gen_2,
        data_dir = dir_2,
        file_stem = stem_2,
        save_file = FALSE)

    use_cols <- c(1,6,7)
    ped1 <- ped1[,use_cols]
    ped2 <- ped2[,use_cols]
    wfu_cols <- c('id','rfid','animalid')
    colnames(ped2) <- wfu_cols
    ped_all <- rbind(ped1, ped2)

    id_map <- data.frame(
        generation = merged_ped$generation,
        merged_id = merged_ped$id,
        accessid = merged_ped$true_id)

    id_map <- merge(id_map, ped_all, by.x='accessid', by.y='id')
    id_map <- id_map[,c('generation','merged_id','accessid','animalid','rfid')]
    id_map <- id_map[order(as.numeric(id_map$merged_id)),]
    id_map <- id_map[!duplicated(id_map$merged_id),]
    
    filename <- paste0(merged_stem, '_id_map.csv')
    outfile <- file.path(out_dir, filename)
    write.csv(id_map, outfile, row.names=F, quote=F, na='')
    cat('ID map written to', outfile, '\n')
    return(id_map)
}
