## PEDIGREE-RELATED FUNCTIONS SPECIFIC TO HS WEST
## written by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

## These are used in conjunction with utils.R, which is more general
## These functions are designed to adhere specifically to data and organizational
## conventions used in-house at the HS West colony. They are likely not useful
## to outside users, though utils.R is of general utility once pedigrees have already
## been properly formatted for use in breedail

source('utils.R')
library(readxl)


# print time-stamped messages to stdout
printout <- function(str) {
    cat(paste0('\n[', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ']'), 
    str, '\n\n')
}

# convert HSW animal IDs to HSW access IDS
animalid_to_accessid <- function(id){
            
    # remove letters & underscores, append '11' to the start of the new ID
    clean_id <- gsub('[A-Za-z_]', '', id)
    accessid <- paste0('11', clean_id)
    return(accessid)
}

# animalid_to_accessid <- function(id, wfu_sl=NULL){
    
#     # check for WFU IDs
#     if (grepl('^(HS|WHS)', id)) {
#         if (is.null(wfu_sl)) {
#             cat('This is a WFU ID. Please provide the WFU shipping list to enable conversion of WFU IDs \n')
#         } else {
#             wfu <- read.csv(wfu_sl)
#             id_row <- wfu[wfu[['sw.id']] == id,]
#             accessid <- id_row[['accessid']]
#         }
#     } else {
#         # remove letters & underscores, append '11' to the start of the new ID
#         clean_id <- gsub('[A-Za-z_]', '', id)
#         accessid <- paste0('11', clean_id)
#     }
#     return(accessid)
# }

# convert HSW access IDs to HSW animal IDS
accessid_to_animalid <- function(id){
    
    # remove letters and underscores
    id <- as.character(id)
    clean_id <- substr(id, start = 3, stop = nchar(id))
    
    rat_no <- substr(clean_id, start = nchar(clean_id)-1, stop = nchar(clean_id))
    pair_no <- substr(clean_id, start = nchar(clean_id)-3, stop = nchar(clean_id)-2)
    gen_no <- unlist(strsplit(clean_id, pair_no))[1]
    
    # construct the animal ID
    animalid <- paste0('G', gen_no, '_', 'B', pair_no, '_', rat_no)
    
    return(animalid)
}

# convert WFU SW.ID to access ID
swid_to_accessid <- function(id, wfu_map) {

    accessid <- wfu_map[wfu_map[['swid']] == id,][['accessid']]
    return(accessid)
}

# convert access ID to WFU SW.ID
swid_to_accessid <- function(id, wfu_map) {

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

read_wfu_raw_ped <- function(ped)
{
    cat('Function: read_wfu_raw_ped() \n')

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

    return(wfu)
}

# split a raw WFU pedigree into single-generation csv files, still in raw format
split_wfu_raw_ped <- function(
    ped,    # the complete pedigree, either xlsx or csv format
    outdir) # the desired directory in which to save per-gen raw pedigree files
{
    
    cat('Function: split_wfu_raw_ped() \n')
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
        if (nchar(gen) == 1) {
            gen <- paste0('0', gen)
        }
        write.csv(ped, file.path(outdir, paste0('wfu_raw_gen', gen, '.csv')),
                  row.names=F, quote=F, na='')
    }
}


# function to format a raw complete pedigree from WFU for use in breedail
format_wfu_raw_ped <- function(
    ped,    # path to raw pedigree file (xlsx or csv)
    keep_cols = NULL,    # vector of column names to keep
    new_colnames = NULL, # vector of desired new column names, in same order as keep_cols
    outdir) # desired output directory path
{

    cat('Function: format_wfu_raw_ped() \n')
    wfu <- read_wfu_raw_ped(ped)

    # subset only the rows needed to identify breeders using breedail
    if (!is.null(keep_cols)) {
        wfu <- wfu[,keep_cols]
    }
    # rename columns for compatibility with breedail
    if (!is.null(new_colnames)) {
        colnames(wfu) <- new_colnames
    }

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

        write.csv(ped, file.path(outdir, paste0('wfu_gen', gen, '.csv')),
                  row.names=F, quote=F, na='?')
    }
}


# function to format a raw WFU pedigree into HSW gens 0-41 (54-95) for breedail
wfu_raw_to_hsw <- function(
    keep_cols = NULL,
    new_colnames = NULL,
    ped, 
    outdir) 
{
    cat('Function: wfu_raw_to_hsw() \n')
    wfu <- read_wfu_raw_ped(ped)

    # subset only the rows needed to identify breeders using breedail
    if (!is.null(keep_cols)) {
        wfu <- wfu[,keep_cols]
    }
    # rename columns for compatibility with breedail
    if (!is.null(new_colnames)) {
        colnames(wfu) <- new_colnames
    }

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
    gen_col <- which(colnames(hsw) == 'generation')
    new_col_order <- c(1:gen_col, ncol(hsw), (gen_col+1):ncol(wfu))
    hsw <- hsw[,new_col_order]

    # split the pedigree by generation
    ped_gens <- split(hsw, hsw$generation)

    # save generations to separate files
    if (dir.exists(outdir) == FALSE) {
        dir.create(outdir, showWarnings = TRUE)
    }

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
format_hsw_raw_ped <- function(
    df,                # an HSW assignment sheet (csv) w/ breeder assignments
    wfu_map=NULL,      # the WFU shipment list (csv) with WFU access IDs and animal IDs (SW IDs)
    outdir=NULL,       # directory in which to save the formatted pedigree file
    return_df=FALSE)   # whether to return the final df to the R console
{
    if (class(df) == 'character') {
        df <- read.csv(df, na.str=c('','NA','NaN','nan'))
    } else if (class(df) == 'data.frame') {
        df <- df
    }
    hsw_gen <- df$generation[1]

    # subset the HSW hsw to only breeders, depending on hsw format
    if ('assignment' %in% colnames(df)) {
        df <- df[df$assignment == 'hsw_breeders' & !is.na(df$assignment),]
    } else if ('breeder' %in% colnames(df)) {
        df <- df[df$breeder == 1,]
    } else if ('hsw_breeders' %in% colnames(df)) {
        df <- df[df$hsw_breeders == 1,]
    } else {
        print('Cannot identify an assignment or breeders column to subset')
    }

    # format the HSW pedigree for breedail
    hsw_keep_cols <- c('rfid','generation','animalid','accessid','sex','dam','sire')
    df <- df[,hsw_keep_cols]
    colnames(df) <- c('rfid','generation','animalid','id','sex','dam_animalid','sire_animalid')
    df$wfu_generation <- NA

    # convert all dam animal IDs to access IDs
    df$dam <- sapply(df$dam_animalid, animalid_to_accessid)

    # convert male animal IDs based on HSW vs WFU format
    for (i in 1:nrow(df)) {

        animal_id <- df$sire_animalid[i]
        if (grepl('^G', animal_id)) {
            # process HSW animal IDs as normal
            df$sire[i] <- animalid_to_accessid(animal_id)
        } else if (grepl('^(HS|WHS)', animal_id)) {
            # use wfu_map for WFU IDs
            if (!is.null(wfu_map)) {
                df$sire[i] <- animalid_to_accessid(animal_id, wfu_sl = wfu_map)
            }
        } else if (grepl('^[0-9]+$', animal_id)) {
            # access IDs remain unchanged
            df$sire[i] <- animal_id
        } else if (is.na(animal_id)) {
            df$sire[i] <- '?'
        } else {
            if (is.null(wfu_map)) {
                cat('Error: Cannot read animal ID', animal_id, 'from row', i, '\n')
            } else {
                df$sire[i] <- wfu_map[wfu_map[['sw.id']]==animal_id,]$id
            }
        }
    }
    col_order <- c('id','dam','sire','sex','generation','wfu_generation','rfid',
                   'animalid','dam_animalid','sire_animalid')
    df <- df[,col_order]
    df <- df[order(as.numeric(df$id)),]

    if (!is.null(outdir)){
        
        if (nchar(hsw_gen) == 1) {
            gen <- paste0('00', hsw_gen)
        } else if (nchar(hsw_gen) == 2) {
            gen <- paste0('0', hsw_gen)
        } else if (nchar(hsw_gen) == 3) {
            gen <- hsw_gen
        }
        
        if (dir.exists(outdir) == FALSE) {
            dir.create(outdir, showWarnings = TRUE)
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


# incorporate WFU IDs into the HSW pedigree following an animal transfer
wfu_into_hsw_gen <- function(hsw_raw,   # an HSW assignment sheet (csv) w/ breeder assignments
                             wfu_ss,    # a WFU shipping sheet (xlsx) with IDs to add to the current HSW generation
                             wfu_sheet, # the name of the excel sheet to use from the shipping sheet file
                             ss_gen,    # the generation number being shipped from WFU
                             wfu_prev,  # the raw WFU pedigree from the previous generation (prior to shipment)
                             outdir=NULL,    # directory in which to save the formatted pedigree file
                             return_df=FALSE) # whether to return the final df to the R console
{
    cat('Function: wfu_into_hsw_gen() \n')
    hsw <- read.csv(hsw_raw, na.str=c('','NA','NaN','nan'))
    wfu <- as.data.frame(read_excel(wfu_ss, sheet=wfu_sheet, na=c('NA','','?')))
    wfu <- wfu[!is.na(wfu[['Animal ID']]),]
    # wfu_prev <- read.csv(wfu_prev, na.str=c('','NA','NaN','nan'))
    hsw_gen <- hsw$generation[1]

    # subset the HSW hsw to only breeders, depending on hsw format
    if ('assignment' %in% colnames(hsw)) {
        hsw <- hsw[hsw$assignment == 'hsw_breeders',]
    } else if ('breeder' %in% colnames(hsw)) {
        hsw <- hsw[hsw$breeder == 1,]
    } else if ('hsw_breeders' %in% colnames(hsw)) {
        hsw <- hsw[hsw$hsw_breeders == 1,]
    } else {
        print('Cannot identify an assignment or breeders column to subset')
    }

    # format the HSW hsw for breedail
    hsw_keep_cols <- c('rfid','generation','animalid','accessid','sex','dam','sire')
    hsw <- hsw[,hsw_keep_cols]
    colnames(hsw) <- c('rfid','generation','animalid','id','sex','dam_animalid','sire_animalid')
    hsw$wfu_generation <- NA
    hsw$dam <- sapply(hsw$dam_animalid, animalid_to_accessid)
    hsw$sire <- sapply(hsw$sire_animalid, animalid_to_accessid)
    col_order <- c('id','dam','sire','sex','generation','wfu_generation','rfid',
                   'animalid','dam_animalid','sire_animalid')
    hsw <- hsw[,col_order]
    hsw <- hsw[!is.na(hsw$id),]

    # format the WFU shipping sheet
    wfu_keep_cols <- c('Dam','Sire','Access ID','Animal ID','Transponder ID','Sex')
    wfu <- wfu[,wfu_keep_cols]
    colnames(wfu) <- c('dam','sire','id','animalid','rfid','sex')
    wfu <- wfu[!is.na(wfu$id),]
    wfu$wfu_generation <- ss_gen
    wfu$dam <- gsub('_', '', wfu$dam)
    wfu$sire <- gsub('_', '', wfu$sire)
    wfu$id <- gsub('_', '', wfu$id)
    wfu$dam_animalid <- sapply(wfu$dam, get_wfu_swid, wfu_df = wfu_prev)
    wfu$sire_animalid <- sapply(wfu$sire, get_wfu_swid, wfu_df = wfu_prev)
    wfu$generation <- hsw_gen
    wfu <- wfu[,col_order]
    wfu <- wfu[!is.na(wfu$id),]

    # concatenate HSW and WFU pedigrees
    out <- rbind(hsw, wfu)
    out <- out[order(as.numeric(out$id)),]

        if (!is.null(outdir)){
        
        if (nchar(hsw_gen) == 1) {
            gen <- paste0('00', hsw_gen)
        } else if (nchar(hsw_gen) == 2) {
            gen <- paste0('0', hsw_gen)
        } else if (nchar(hsw_gen) == 3) {
            gen <- hsw_gen
        }
        
        if (dir.exists(outdir) == FALSE) {
            dir.create(outdir, showWarnings = TRUE)
        }

        write.csv(out, file.path(outdir, paste0('hsw_gen', gen, '.csv')), row.names=F, quote=F, na='?')
    }

    if (!is.logical(return_df)) {
        stop("return_df should be a logical value")
    }

    if (return_df) {
        return(out)    
    } 
    
}


# format a breeder file to assist with pairing in the HSW colony
create_breeder_file <- function(
    pairs,          # R dataframe or path to csv, as output by select.breeders, dam/sire must be animal IDs
    df,             # colony dataframe or assignment file
    wfu_ss=NULL,    # WFU shipping sheet, if pairing with shipped WFU rats
    outdir=NULL)
{
    if (class(pairs) == 'data.frame') {
        pairs <- pairs
    } else if (class(pairs) == 'character') {
        pairs <- read.csv(pairs, na.str=c('','NA','NaN','nan'))
    } else if (sum(!names(pairs) %in% c('pairs','file')) == 0) {
        pairs <- pairs$pairs
    }
    if (!is.null(wfu_ss)) {    
        wfu <- as.data.frame(read_excel(wfu_ss))
    }
    df <- read.csv(df, na.str=c('','NA','NaN','nan'))
    gen <- as.numeric(df$generation[1])
    pairs$generation <- gen
    pairs$breederpair <- c(
        paste0('HSW', gen + 1, '_B0', seq.int(1,9)),
        paste0('HSW', gen + 1, '_B', seq.int(10,nrow(pairs)))
    )
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
    pairs$kinship <- round(pairs$kinship, 4)
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
        print(paste('Breeder file written to', outfile))
    }
    return(list(pairs = pairs, file = outfile))
}


# count the parental breeder pairs that have been or currently need assignment/pairing
count_breeder_pairs <- function(
    assignments, # assignment sheet with ALL assignments so far
    breederpairs,   # 'breederpair file' w/ all breeder pairings so far, or the output from create_breeder_file()
    outdir)
{
    if (class(assignments) == 'character') {
        assignments <- read.csv(assignments, na.str=(c('','NA','NaN','nan')))
    } else if (class(assignments) == 'data.frame') {
        assignments <- assignments
    }
    if (class(breederpairs)=='character'){
        all_pairs <- read.csv(breederpair, na.str=c('','NA','NaN','nan'))    
    } else {
        breederpairs <- breederpairs
    }
    assigned_breeders <- assignments[assignments$assignment=='hsw_breeders',]
    paired_so_far <- breederpairs[(breederpairs$paired_dam != 'NONE') | (breederpairs$paired_sire != 'NONE'),]

    # all animal IDs that have been assigned and paired so far
    assigned <- assigned_breeders$animalid
    paired <- c(paired_so_far$dam_animalid, paired_so_far$sire_animalid)

    all_f <- as.matrix(table(assignments$breederpair, assignments$sex))[,1]
    all_m <- as.matrix(table(assignments$breederpair, assignments$sex))[,2]
    all_breederpairs <- data.frame(
        breederpair = names(all_f),
        females = all_f,
        males = all_m)

    no_females <- all_breederpairs[all_breederpairs$females==0,]$breederpair
    no_males <- all_breederpairs[all_breederpairs$males==0,]$breederpair

    paired.but.not.assigned <- setdiff(paired, assigned)
    assigned.but.not.paired <- setdiff(assigned, paired)

    # dataframe counting all breederpairs represented in assignments
    f_assigned <- as.matrix(table(assigned_breeders$breederpair, assigned_breeders$sex))[,1]
    m_assigned <- as.matrix(table(assigned_breeders$breederpair, assigned_breeders$sex))[,2]
    assigned.breederpairs <- data.frame(
        generation = assignments$generation[1]-1,
        breederpair = names(f_assigned),
        females_assigned = f_assigned,
        males_assigned = m_assigned)

    # dataframe counting all breederpairs represented in pairings
    paired_assignments <- assigned_breeders[assigned_breeders$animalid %in% paired,]
    if (nrow(paired_assignments) > 0) {
        f_paired <- as.matrix(table(paired_assignments$breederpair, paired_assignments$sex))[,1]
        m_paired <- as.matrix(table(paired_assignments$breederpair, paired_assignments$sex))[,2]
        paired.breederpairs <- data.frame(
            generation = assignments$generation[1]-1,
            breederpair = names(f_paired),
            females_paired = f_paired,
            males_paired = m_paired)
    } else {paired.breederpairs <- NULL}

    # dataframe counting bps lacking offspring assignments to hsw_breeders
    breederpairs.need.assignment <- 
        assigned.breederpairs[(assigned.breederpairs$females_assigned < 1) |
                                       (assigned.breederpairs$males_assigned < 1),]

    # remove IDs from needs.pairing/needs.assignment lists if there were no M or F to begin with
    needs_assignment <- c()
    for (i in 1:nrow(breederpairs.need.assignment)) {
        pair <- breederpairs.need.assignment$breederpair[i]
        total_males <- all_breederpairs[all_breederpairs$breederpair==pair,]$males
        total_females <- all_breederpairs[all_breederpairs$breederpair==pair,]$females
        males_assigned <- assigned.breederpairs[assigned.breederpairs$breederpair==pair,]$males_assigned
        females_assigned <- assigned.breederpairs[assigned.breederpairs$breederpair==pair,]$females_assigned
    
        if ((total_males > 0) & (males_assigned < 1)) {
            needs_assignment <- c(needs_assignment, pair)
        }
        if ((total_females > 0) & (females_assigned < 1)) {
            needs_assignment <- c(needs_assignment, pair)
        }
    }
    breederpairs.need.assignment <- breederpairs.need.assignment[needs_assignment,]

    # dataframe counting bps with offspring assigned to hsw_breeders that still need to be paired
    if (!is.null(paired.breederpairs)) {
        breederpairs.need.pairing <- 
            paired.breederpairs[(paired.breederpairs$females_paired < 1) |
                                        (paired.breederpairs$males_paired < 1),]
    } else {
        breederpairs.need.pairing <- assigned.breederpairs
    }

    needs_pairing <- c()
    for (i in 1:nrow(breederpairs.need.pairing)) {
        pair <- breederpairs.need.pairing$breederpair[i]
        total_males <- all_breederpairs[all_breederpairs$breederpair==pair,]$males
        total_females <- all_breederpairs[all_breederpairs$breederpair==pair,]$females
        males_paired <- paired.breederpairs[paired.breederpairs$breederpair==pair,]$males_paired
        females_paired <- paired.breederpairs[paired.breederpairs$breederpair==pair,]$females_paired
    
        if (!is.null(males_paired)) {
            if ((total_males > 0) & (males_paired < 1)) {
                needs_pairing <- c(needs_pairing, pair)
            }
        }
        if (!is.null(females_paired)) {
            if ((total_females > 0) & (!is.null(females_paired)) & (females_paired < 1)) {
                needs_pairing <- c(needs_pairing, pair)
            }
        }
    }
    if (!is.null(needs_pairing)) {
        breederpairs.need.pairing <- breederpairs.need.pairing[beederpairs.need.paring$breederpair %in% needs_pairing,]
    }

    still.avail.for.assignment <- assignments[assignments$assignment=='not_assigned',]
    
    if (!is.null(outdir)) {

        datestamp <- format(Sys.time(),'%Y%m%d')
        outdir <- file.path(outdir, 'bp_counts')
        outfile1 <- file.path(outdir, paste0('still_avail_for_assignment_', datestamp, '.csv'))
        dir.create(outdir, showWarnings=F)
        write.csv(still.avail.for.assignment, outfile1, row.names=F, quote=F, na='')

        still_avail_counts <- table(still.avail.for.assignment$breederpair, still.avail.for.assignment$sex)
        outfile2 <- file.path(outdir, paste0('still_avail_bp_counts_', datestamp, '.csv'))
        write.csv(still_avail_counts, outfile2, row.names=T, quote=F, na='')
    }
    out <- list(all.breederpairs = all_breederpairs,
                paired.but.not.assigned = paired.but.not.assigned, 
                assigned.but.not.paired = assigned.but.not.paired, 
                assigned.breederpairs = assigned.breederpairs,
                paired.breederpairs = paired.breederpairs,
                breederpairs.need.assignment = breederpairs.need.assignment,
                breederpairs.need.pairing = breederpairs.need.pairing,
                still.avail.for.assignment = still.avail.for.assignment$animalid)

    return(out)
}


# add proposed alternative breeder pairs to an existing breeder file
add_to_breeder_file <- function(
    orig_pairs,   # R dataframe or csv path: breeder file output by create_breeder_file()
    new_pairings, # R list output by get.alternative.pairings()
    df,           # colony dataframe or assignment file
    n_new,        # int, the number of desired new pairings to add
    wfu_ss=NULL,  # WFU shipping sheet, if pairing with shipped WFU rats
    outdir=NULL)
{
    if (class(orig_pairs) == 'data.frame') {
        orig_pairs <- orig_pairs
    } else if (class(orig_pairs) == 'character') {
        orig_pairs <- read.csv(orig_pairs, na.str=c('','NA','NaN','nan'))
    }
    if (!is.null(wfu_ss)) {    
        wfu <- as.data.frame(read_excel(wfu_ss))
    }
    new_pairings <- new_pairings$best.pairs  
    df <- read.csv(df, na.str=c('','NA','NaN','nan'))
    names(new_pairings)[1:2] <- c('dam_animalid','sire_animalid')
    gen <- as.numeric(df$generation[1])
    n_paired <- nrow(orig_pairs)
    total_pairs <- n_paired + n_new

    all_males <- unique(new_pairings$sire_animalid)
    all_females <- unique(new_pairings$dam_animalid)

    # subset the number of desired new pairings with minimum kinship
    for (i in 1:n_new) {
        if (i == 1) {
            pairs_to_add <- new_pairings[new_pairings$kinship == min(new_pairings$kinship),]
            unavail_ids <- c(pairs_to_add$dam_animalid, pairs_to_add$sire_animalid)
            available_f <- new_pairings[!new_pairings$dam_animalid %in% unavail_ids,] 
            available_m <- new_pairings[!new_pairings$sire_animalid %in% unavail_ids,]
            available_pairings <- rbind(available_f, available_m)
            available_pairings <- available_pairings[!duplicated(available_pairings),]
            available_pairings <- new_pairings[!new_pairings$dam_animalid %in% unavail_ids,]
            available_pairings <- available_pairings[!available_pairings$sire_animalid %in% unavail_ids,]
        }
        if (i > 1) {
            if (nrow(available_pairings) == 0) {
                stop(paste('No further combinations are available for pairing.', 
                           i-1, 'pairs successfully added'))
            }
            next_pair <- available_pairings[available_pairings$kinship == 
                                            min(available_pairings$kinship),]
            pairs_to_add <- rbind(pairs_to_add, next_pair)
            unavail_ids <- unique(c(unavail_ids, pairs_to_add$dam_animalid, pairs_to_add$sire_animalid))
            available_f <- available_pairings[!(new_pairings$dam_animalid %in% unavail_ids),] 
            available_m <- new_pairings[!(new_pairings$sire_animalid %in% unavail_ids),]
            available_pairings <- rbind(available_f, available_m)
        }
    }

    for (i in 1:nrow(pairs_to_add)) {
        
        dam <- pairs_to_add$dam_animalid[i]
        
        if (substr(dam,1,1)=='G') {
            dam_df <- df[df$animalid==dam,]
            pairs_to_add$dam_rfid[i] <- as.character(dam_df$rfid)
            pairs_to_add$dam_earpunch[i] <- dam_df$earpunch
            pairs_to_add$dam_coatcolor[i] <- dam_df$coatcolor
            pairs_to_add$dam_rack_num[i] <- dam_df$rack_number
            pairs_to_add$dam_rack_pos[i] <- dam_df$rack_position            
        } else if (substr(dam,1,1)=='H') {
            dam_df <- wfu[wfu[['Animal ID']]==dam,]
            pairs_to_add$dam_rfid[i] <- dam_df[['Transponder ID']]
            pairs_to_add$dam_earpunch[i] <- dam_df[['Ear Punch']]
            pairs_to_add$dam_coatcolor[i] <- dam_df[['Coat Color']]
            pairs_to_add$dam_rack_num[i] <- NA
            pairs_to_add$dam_rack_pos[i] <- NA
            pairs_to_add$dam_wfu_cage[i] <- dam_df[['Ship Box']]
        }

        sire <- pairs_to_add$sire_animalid[i]
        
        if (substr(sire,1,1)=='G') {
            sire_df <- df[df$animalid==sire,]
            pairs_to_add$sire_rfid[i] <- as.character(sire_df$rfid)
            pairs_to_add$sire_earpunch[i] <- sire_df$earpunch
            pairs_to_add$sire_coatcolor[i] <- sire_df$coatcolor
            pairs_to_add$sire_rack_num[i] <- sire_df$rack_number
            pairs_to_add$sire_rack_pos[i] <- sire_df$rack_position            
        } else if (substr(sire,1,1)=='H') {
            sire_df <- wfu[wfu[['Animal ID']]==sire,]
            pairs_to_add$sire_rfid[i] <- sire_df[['Transponder ID']]
            pairs_to_add$sire_earpunch[i] <- sire_df[['Ear Punch']]
            pairs_to_add$sire_coatcolor[i] <- sire_df[['Coat Color']]
            pairs_to_add$sire_rack_num[i] <- NA
            pairs_to_add$sire_rack_pos[i] <- NA
            pairs_to_add$sire_wfu_cage[i] <- sire_df[['Ship Box']]
        }
    }

    pairs_to_add$generation <- gen
    if (n_paired < 10) {
        pairs_to_add$breederpair <- c(
            paste0('HSW', gen + 1, '_B0', seq.int((n_paired+1),9)),
            paste0('HSW', gen + 1, '_B', seq.int(10,total_pairs))
        )
    } else if (n_paired >= 10) {
        pairs_to_add$breederpair <- paste0('HSW', gen + 1, '_B', seq.int((n_paired+1),total_pairs))
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
    
    pairs_to_add <- pairs_to_add[,col_order]
    pairs_to_add$kinship <- round(pairs_to_add$kinship, 4)
    all_pairs <- rbind(orig_pairs, pairs_to_add)
    all_pairs_slim <- all_pairs[,c('generation','breederpair','kinship','dam_rfid','sire_rfid',
                                    'dam_animalid','sire_animalid'),]

    if (!is.null(outdir)) {
        datestamp <- format(Sys.time(),'%Y%m%d')
        outfile <- paste0('hsw_gen', gen, '_', gen+1, '_parents_', datestamp, '.csv')
        slimfile <- gsub('.csv', '_ids.csv', outfile)
        colony_file <- gsub('.csv','_colony.csv',outfile)
        outfile <- file.path(outdir, outfile)
        slimfile <- file.path(outdir, slimfile)
        colony_file <- file.path(outdir, colony_file)
        write.csv(all_pairs, outfile, row.names=F, quote=F, na='')
        write.csv(all_pairs, colony_file, row.names=F, quote=F, na='')
        write.csv(all_pairs_slim, slimfile, row.names=F, quote=F, na='')
        print(paste('Updated breeder file written to', outfile))
    }

    return(all_pairs)
}

# get up to 6 next-best breeder pairings for a given animal ID
next_best_pairs <- function(
    id, # animalid of the breeder that needs a new mate pairing
    kinship, # pairwise kinship df or path to csv file
    available_ids=NULL, # vector of animalids, or txt file of animalids
    available_ids_file=NULL) # path to text file of animalids, one per line, no header
{
    if (class(kinship) == 'character') {
        assignments <- read.csv(assignments, na.str=(c('','NA','NaN','nan')))
    } else if (class(kinship) == 'data.frame') {
        kinship <- kinship
    }
    if (!is.null(available_ids)) {
        available_ids <- available_ids
    } 
    if (!is.null(available_ids_file)) {
        available_ids <- readLines(available_ids_file)
    }

    dam <- ifelse(id %in% kinship$dam, TRUE, FALSE)
    sire <- ifelse(id %in% kinship$sire, TRUE, FALSE)
    if (sum(c(dam, sire)) == 0){
        stop(paste('ID', id, 'is not in the kinship matrix'))
    }

    if (dam) {
        avail_pairs <- kinship[kinship$dam == id,]
    } else {
        avail_pairs <- kinship[kinship$sire == id,]
    }

    avail_pairs <- avail_pairs[order(avail_pairs$kinship),]
    
    if (!is.null(available_ids)) {
        avail_pairs <- avail_pairs[avail_pairs$dam %in% available_ids | avail_pairs$sire %in% available_ids,]
    }

    if (nrow(avail_pairs) > 6) {
        avail_pairs <- avail_pairs[1:6,]
    }
    return(avail_pairs)
}