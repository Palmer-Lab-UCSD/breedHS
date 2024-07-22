## PEDIGREE-RELATED FUNCTIONS SPECIFIC TO HS WEST
## written by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

## These are used in conjunction with utils.R, which is more general
## These functions are designed to adhere specifically to data and organizational
## conventions used in-house at the HS West colony. They are likely not useful
## to outside users, though utils.R is of general utility once pedigrees have already
## been properly formatted for use in breedail

source('utils.R')
library(readxl)


# convert HSW animal IDs to HSW access IDS
animalid_to_accessid <- function(id){
            
    # remove letters & underscores, append '11' to the start of the new ID
    clean_id <- gsub('[A-Za-z_]', '', id)
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

# find a WFU animal ID (SW.ID) given a WFU access ID
get_wfu_swid <- function(id, wfu_df){
    wfu_df <- read.csv(wfu_df)
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
hsw_df_to_ped <- function(df,              # an HSW assignment sheet (csv) w/ breeder assignments
                          prev_ped=NULL,   # the formatted HSW pedigree for the previous generation, used if WFU IDs are present
                          outdir=NULL,     # directory in which to save the formatted pedigree file
                          return_df=FALSE) # whether to return the final df to the R console
{
    df <- read.csv(df)
    prev_ped <- read.csv(prev_ped)
    hsw_gen <- df$generation[1]

    # subset the HSW hsw to only breeders, depending on hsw format
    if ('assignment' %in% colnames(df)) {
        df <- df[df$assignment == 'hsw_breeders',]
    } else if ('breeder' %in% colnames(df)) {
        df <- df[df$breeder == 1,]
    } else if ('hsw_breeders' %in% colnames(df)) {
        df <- df[df$hsw_breeders == 1,]
    } else {
        print('Cannot identify an assignment or breeders column to subset')
    }

    # format the HSW hsw for breedail
    hsw_keep_cols <- c('rfid','generation','animalid','accessid','sex','dam','sire')
    df <- df[,hsw_keep_cols]
    colnames(df) <- c('rfid','generation','animalid','id','sex','dam_animalid','sire_animalid')
    df$wfu_generation <- NA

    # convert all dam animal IDs to access IDs
    df$dam <- sapply(df$dam_animalid, animalid_to_accessid)

    # convert male animal IDs based on HSW vs WFU format
    for (i in 1:nrow(df)) {
        animal_id <- df$sire_animalid[i]
        if (substr(animal_id,1,1) == 'G') {
            df$sire[i] <- animalid_to_accessid(animal_id)
        } else {
            if (is.null(prev_ped)) {
                print('Please provide a previous pedigree file from which to extract WFU access IDs')
            }
            df$sire[i] <- prev_ped[prev_ped$animalid==animal_id,]$id
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
    hsw <- read.csv(hsw_raw)
    wfu <- as.data.frame(read_excel(wfu_ss, sheet=wfu_sheet))
    wfu_prev <- read.csv(wfu_prev)
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

    # format the WFU shipping sheet
    wfu_keep_cols <- c('Dam','Sire','Access ID','Animal ID','Transponder ID','Sex')
    wfu <- wfu[,wfu_keep_cols]
    colnames(wfu) <- c('dam','sire','id','animalid','rfid','sex')
    wfu <- wfu[!is.na(wfu$id),]
    wfu$wfu_generation <- ss_gen
    wfu$dam <- gsub('_', '', wfu$dam)
    wfu$sire <- gsub('_', '', wfu$sire)
    wfu$id <- gsub('_', '', wfu$id)
    wfu$dam_animalid <- sapply(wfu$dam, get_wfu_swid, wfu_df = wfu_46)
    wfu$sire_animalid <- sapply(wfu$sire, get_wfu_swid, wfu_df = wfu_46)
    wfu$generation <- hsw_gen
    wfu <- wfu[,col_order]

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
        
        write.csv(out, file.path(outdir, paste0('hsw_gen', gen, '.csv')), row.names=F, quote=F, na='?')
    }

    if (!is.logical(return_df)) {
        stop("return_df should be a logical value")
    }

    if (return_df) {
        return(out)    
    } 
    
}


create_breeder_file <- function(
    pairs,          # R dataframe or path to csv, as output by select.breeders, dam/sire must be animal IDs
    df,             # colony dataframe or assignment file
    wfu_ss=NULL,    # WFU shipping sheet, if pairing with shipped WFU rats
    outdir=NULL)
{
    if (class(pairs) == 'data.frame') {
        pairs <- pairs
    } else if (class(pairs) == 'character') {
        pairs <- read.csv(pairs)
    }
    if (!is.null(wfu_ss)) {    
        wfu <- as.data.frame(read_excel(wfu_ss))
    }
    df <- read.csv(df)
    gen <- as.numeric(df$generation[1])
    pairs$generation <- gen
    pairs$breederpair <- c(
        paste0('G', gen + 1, '_B0', seq.int(1,9)),
        paste0('G', gen + 1, '_B', seq.int(10,nrow(pairs)))
    )
    pairs$sire_animalid <- pairs$sire
    pairs$dam_animalid <- pairs$dam

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
        dam_cols <- c('generation','breederpair','kinship','dam_rfid','dam_animalid','dam_earpunch',
                       'dam_coatcolor','dam_rack_num','dam_rack_pos')
        sire_cols <- c('generation','breederpair','kinship','sire_rfid','sire_animalid','sire_earpunch',
                       'sire_coatcolor','sire_rack_num','sire_rack_pos')

    } else {
        col_order <- c('generation','breederpair','kinship','dam_rfid','dam_animalid','dam_earpunch',
                       'dam_coatcolor','dam_rack_num','dam_rack_pos','dam_wfu_cage',
                       'sire_rfid','sire_animalid','sire_earpunch','sire_coatcolor',
                       'sire_rack_num','sire_rack_pos','sire_wfu_cage')
        dam_cols <- c('generation','breederpair','kinship','dam_rfid','dam_animalid','dam_earpunch',
                       'dam_coatcolor','dam_rack_num','dam_rack_pos','dam_wfu_cage')
        sire_cols <- c('generation','breederpair','kinship','sire_rfid','sire_animalid','sire_earpunch',
                       'sire_coatcolor','sire_rack_num','sire_rack_pos','sire_wfu_cage')

    }
    
    pairs <- pairs[,col_order]
    pairs$kinship <- round(pairs$kinship, 4)
    if (!is.null(outdir)) {
        datestamp <- format(Sys.time(),'%Y%m%d')
        outfile <- paste0('hsw_gen', gen, '_', gen+1, '_parents_', datestamp, '.csv')
        write.csv(pairs, file.path(outdir, outfile), row.names=F, quote=F, na='')
    }
    return(pairs)
}


# count the parental breeder pairs that have been or currently need assignment/pairing
count_breeder_pairs <- function(
    assignments, # assignment sheet with ALL assignments so far
    paired_breeders)   # 'breeder file': all breeder pairings so far, or output from create_breeder_file()
{
    assignments <- read.csv(assignments, na.str=(''))
    assigned_breeders <- assignments[assignments$hsw_breeders==1,]
    if (class(paired_breeders)=='str'){
        paired_breeders <- read.csv(paired_breeders)    
    } else {
        paired_breeders <- paired_breeders
    }

    # all animal IDs that have been assigned and paired so far
    assigned <- assigned_breeders$animalid
    paired <- c(paired_breeders$dam_animalid, paired_breeders$sire_animalid)

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
    f_paired <- as.matrix(table(paired_assignments$breederpair, paired_assignments$sex))[,1]
    m_paired <- as.matrix(table(paired_assignments$breederpair, paired_assignments$sex))[,2]
    paired.breederpairs <- data.frame(
        generation = assignments$generation[1]-1,
        breederpair = names(f_paired),
        females_paired = f_paired,
        males_paired = m_paired)

    breederpairs.need.assignment <- 
        assigned.breederpairs[(assigned.breederpairs$females_assigned < 1) |
                                       (assigned.breederpairs$males_assigned < 1),]

    breederpairs.need.pairing <- 
        paired.breederpairs[(paired.breederpairs$females_paired < 1) |
                                       (paired.breederpairs$males_paired < 1),]

    all_f <- as.matrix(table(assignments$breederpair, assignments$sex))[,1]
    all_m <- as.matrix(table(assignments$breederpair, assignments$sex))[,2]
    all_breederpairs <- data.frame(
        breederpair = names(all_f),
        females = all_f,
        males = all_m
    )
    no_females <- all_breederpairs[all_breederpairs$females==0,]$breederpair
    no_males <- all_breederpairs[all_breederpairs$males==0,]$breederpair
    
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
    
    needs_pairing <- c()
    for (i in 1:nrow(breederpairs.need.pairing)) {
        pair <- breederpairs.need.pairing$breederpair[i]
        total_males <- all_breederpairs[all_breederpairs$breederpair==pair,]$males
        total_females <- all_breederpairs[all_breederpairs$breederpair==pair,]$females
        males_paired <- paired.breederpairs[paired.breederpairs$breederpair==pair,]$males_paired
        females_paired <- paired.breederpairs[paired.breederpairs$breederpair==pair,]$females_paired
    
        if ((total_males > 0) & (males_paired < 1)) {
            needs_pairing <- c(needs_pairing, pair)
        }
        if ((total_females > 0) & (females_paired < 1)) {
            needs_pairing <- c(needs_pairing, pair)
        }
    }
    breederpairs.need.pairing <- breederpairs.need.pairing[needs_pairing,]

    still.avail.for.assignment <- assignments[is.na(assignments$assignment),]$animalid
    
    out <- list(all.breederpairs = all_breederpairs,
                paired.but.not.assigned = paired.but.not.assigned, 
                assigned.but.not.paired = assigned.but.not.paired, 
                assigned.breederpairs = assigned.breederpairs,
                paired.breederpairs = paired.breederpairs,
                breederpairs.need.assignment = breederpairs.need.assignment,
                breederpairs.need.pairing = breederpairs.need.pairing,
                still.avail.for.assignment = still.avail.for.assignment)

    return(out)
}