## PEDIGREE-RELATED FUNCTIONS SPECIFIC TO HS WEST
## written by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

## These are used in conjunction with utils.R, which is more general
## These functions are designed to adhere specifically to data and organizational
## conventions used in-house at the HS West colony. They are likely not useful
## to outside users, though utils.R is of general utility once pedigrees have already
## been properly formatted for use in breedail

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
    clean_id <- substr(id, start = 3, stop = nchar(id))
    
    rat_no <- substr(clean_id, start = nchar(clean_id)-1, stop = nchar(clean_id))
    pair_no <- substr(clean_id, start = nchar(clean_id)-3, stop = nchar(clean_id)-2)
    gen_no <- unlist(strsplit(clean_id, pair_no))[1]
    
    # construct the animal ID
    animalid <- paste0('G', gen_no, '_', 'B', pair_no, '_', rat_no)
    
    return(animalid)
}



# function to format a raw WFU pedigree into HSW gens 0-41 (54-95) for breedail
wfu_raw_to_hsw <- function(
    ped, 
    outdir) 
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

    # # reorder columns
    # gen_col <- which(colnames(hsw) == 'generation')
    # new_col_order <- c(1:gen_col, ncol(hsw), (gen_col+1):ncol(wfu))
    # hsw <- hsw[,new_col_order]

    # split the pedigree by generation
    ped_gens <- split(hsw, hsw$generation)

    # save generations to separate files
    if (dir.exists(outdir) == FALSE) {
        dir.create(outdir, showWarnings = TRUE)
    }

    # exit the function for generations > 41 (whose length==0)
    if (length(ped_gens)==0) {
        return(NULL)
    } else {
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
}


# function to produce a raw pedigree file from an HSW assignment sheet
assignment_to_raw_ped <- function(
    assignments, # assignment file csv or dataframe
    outdir,
    wfu_ss = NULL, # wfu shipping sheet to incorporate into the pedigree, path to xlsx
    wfu_sheet = NULL
) 
{
    if (class(assignments)=='character') {
        assignments <- read.csv(assignments)
    }
    ped_cols <- c('generation','rfid','animalid','accessid','sex','coatcolor','earpunch','dam','sire','comments')
    gen <- assignments$generation[1]

    ped <- assignments[assignments$assignment=='hsw_breeders',]
    ped <- ped[,ped_cols]
    ped <- ped[order(ped$animalid),]

    if (!is.null(wfu_ss)) {
        
        if(class(wfu_ss)=='character'){
        wfu_ss <- as.data.frame(read_excel(wfu_ss, sheet=wfu_sheet))
        }

        wfu <- data.frame(
            generation = gen,
            rfid = wfu_ss['Transponder ID'],
            animalid = wfu_ss['Animal ID'],
            accessid = gsub('_', '', wfu_ss['Access ID']),
            sex = wfu_ss['Sex'],
            coatcolor = wfu_ss['Coat Color'],
            earpunch = wfu_ss['Ear Punch'],
            dam = gsub('_', '', wfu_ss['Dam']), 
            sire = gsub('_', '', wfu_ss['Sire']), 
            comments = NA
        )
        wfu <- wfu[order(wfu$animalid)]
        ped <- rbind(ped, wfu)
    }
    
    outfile <-  file.path(outdir, paste0('hsw_raw_gen', gen, '.csv')) 
    write.csv(ped, outfile, row.names=F, quote=F, na='')
    cat('Raw pedigree written to', outfile, '\n')
    
    return(outfile)
}

# function to produce a breedail pedigree from breeders assignments
format_hsw_raw_ped <- function(
    df,                 # an HSW raw pedigree (csv) w/ an assignment column including 'hsw_breeders' assignments
    wfu_map = NULL,     # the WFU ID map (csv) with WFU access IDs and animal IDs (SW IDs)
    hsw_map = NULL,     # the HSW ID map (csv) 
    outdir = NULL,      # directory in which to save the formatted pedigree file
    return_df = FALSE)  # whether to return the final df to the R console
{
    if (class(df) == 'character') {
        df <- read.csv(df, na.str=c('','NA','NaN','nan'))
    } else if (class(df) == 'data.frame') {
        df <- df
    }
    hsw_gen <- df$generation[1]
    if (!is.null(wfu_map)) {
        if(class(wfu_map) == 'character') {
            wfu_map <- read.csv(wfu_map)
        }
    }
    if (!is.null(hsw_map)) {
        if(class(hsw_map) == 'character') {
            hsw_map <- read.csv(hsw_map)
        }
    }
    
    # format the HSW pedigree for breedail
    hsw_keep_cols <- c('rfid','generation','animalid','accessid','sex','dam','sire')
    df <- df[,hsw_keep_cols]
    colnames(df) <- c('rfid','generation','animalid','id','sex','dam_animalid','sire_animalid')
    df$wfu_generation <- NA

    # convert dam animal IDs based on HSW vs WFU format
    df$dam <- unlist(sapply(df$dam_animalid, function(x) {
        if (is.na(x)) return(NA)
        if (grepl('^WHS', x)) {
            swid_to_accessid(x, wfu_map)
        } else if (grepl('^HSW', x)) {
            hsw_map[hsw_map[['animalid']] == x,][['accessid']]
        } else if (grepl('^G', x)) {
            if (grepl('^G01', x)) {
                hsw_map[hsw_map[['animalid']] == x,][['accessid']]
            } else {
                animalid_to_accessid(x)
            }
        } else {
            NA
        }
    }, USE.NAMES = FALSE))

    # convert sire IDs
    df$sire <- unlist(sapply(df$sire_animalid, function(x) {
        if (is.na(x)) return(NA)
        if (grepl('^WHS', x)) {
            swid_to_accessid(x, wfu_map)
        } else if (grepl('^HSW', x)) {
            hsw_map[hsw_map[['animalid']] == x,][['accessid']]
        } else if (grepl('^G', x)) {
            if (grepl('^G01', x)) {
                hsw_map[hsw_map[['animalid']] == x,][['accessid']]
            } else {
                animalid_to_accessid(x)
            }
        } else {
            NA
        }
    }, USE.NAMES = FALSE))

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
        cat('Cannot identify an assignment or breeders column to subset \n')
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
create_hsw_breeder_file <- function(
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


# count the parental breeder pairs that have been or currently need assignment/pairing
count_breeder_pairs <- function(
    assignments, # assignment sheet with ALL assignments so far
    breederpairs)   # 'breederpair file' w/ all breeder pairings so far, or the output from create_breeder_file()
{
    if (class(assignments) == 'character') {
        assignments <- read.csv(assignments, na.str=(c('','NA','NaN','nan')))
    } else if (class(assignments) == 'data.frame') {
        assignments <- assignments
    }
    if (class(breederpairs)=='character'){
        breederpairs <- read.csv(breederpairs, na.str=c('','NA','NaN','nan'))    
    } else {
        breederpairs <- breederpairs
    }
    assigned_breeders <- assignments[assignments$assignment=='hsw_breeders',]
    
    # all animal IDs that have been assigned and paired so far
    assigned <- assigned_breeders$animalid
    paired <- c(breederpairs$dam_animalid, breederpairs$sire_animalid)

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

    if (length(paired.but.not.assigned) == 0) {paired.but.not.assigned <- NULL}
    if (length(assigned.but.not.paired) == 0) {assigned.but.not.paired <- NULL}

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

    if (nrow(breederpairs.need.assignment) > 0) {
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
    } else {
        breederpairs.need.assignment <- NULL
    }

    # dataframe counting bps with offspring assigned to hsw_breeders that still need to be paired
    if (!is.null(paired.breederpairs)) {
        breederpairs.need.pairing <- 
            paired.breederpairs[(paired.breederpairs$females_paired < 1) |
                                        (paired.breederpairs$males_paired < 1),]
    } else {
        breederpairs.need.pairing <- assigned.breederpairs
    }

    if (nrow(breederpairs.need.pairing) > 0) {
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
            breederpairs.need.pairing <- breederpairs.need.pairing[breederpairs.need.pairing$breederpair %in% needs_pairing,]
        }
    } else {
        breederpairs.need.pairing <- NULL
    }
    
    out <- list(all.breederpairs = all_breederpairs,
                paired.but.not.assigned = paired.but.not.assigned, 
                assigned.but.not.paired = assigned.but.not.paired, 
                assigned.breederpairs = assigned.breederpairs,
                paired.breederpairs = paired.breederpairs,
                breederpairs.need.assignment = breederpairs.need.assignment,
                breederpairs.need.pairing = breederpairs.need.pairing)

    return(out)
}


# get n next-best replacements for a given set of animal IDs, or the n best pairings from all available IDs
best_alt_pairs <- function(
    pairs,  # breederpair file
    kinship, # kinship for all possible pairs, df or csv path
    avail_ids, # vector or path to animalids still available for pairing
    n_best = 4, # the number of next-best IDs to return per ID that needs replacing
    ids_to_replace=NULL, # vector or path
    ids_to_pair=NULL, # vector or path
    k_all=FALSE, # whether to return a kinship matrix for all available rats
    outdir
) {
    if (class(pairs)=='character') {
        pairs <- read.csv(pairs)
    }
    if (class(kinship)=='character') {
        kinship <- read.csv(kinship)
    } 
    if (length(avail_ids)==1 && file.exists(avail_ids)) {
        avail_ids <- readLines(avail_ids)
    } else { cat('Please provide a list of animal IDs available for pairing')}
    if (!is.null(ids_to_replace)) {
        if (length(ids_to_replace)==1 && file.exists(ids_to_replace)) {
            ids_to_replace <- readLines(ids_to_replace)
        } 
    }
    if (!is.null(ids_to_pair)) {
        if (length(ids_to_pair)==1 && file.exists(ids_to_pair)) {
            ids_to_pair <- readLines(ids_to_pair)
        }
    }
    datestamp <- format(Sys.time(),'%Y%m%d')
    out <- list()

    # get the breederpair combination for each possible pairing
    kinship$dam_fam <- sapply(kinship$dam_fam, function(x) {
        if (nchar(x)==1) return(paste0('0',x)) else return(x)})
    kinship$sire_fam <- sapply(kinship$sire_fam, function(x) {
        if (nchar(x)==1) return(paste0('0',x)) else return(x)})

    kinship$combo <- paste0(kinship$dam_fam, '-', kinship$sire_fam)

    # get breederpair for each paired dam or sire
    paired_dams <- pairs$paired_dam
    paired_sires <- pairs$paired_sire
    proposed_dams <- pairs$dam_animalid
    proposed_sires <- pairs$sire_animalid
    use_dam_ids <- ifelse(paired_dams=='NONE', proposed_dams, paired_dams)
    use_sire_ids <- ifelse(paired_sires=='NONE', proposed_sires, paired_sires)
    paired_dam_fam <- sapply(use_dam_ids, function(x) sub('.*_B(\\d+)_.*', '\\1', x))
    paired_sire_fam <- sapply(use_sire_ids, function(x) sub('.*_B(\\d+)_.*', '\\1', x))

    # breederpair combo for all pairs
    paired_combos <- paste0(paired_dam_fam, '-', paired_sire_fam)
    reverse_pairs <- sapply(strsplit(paired_combos, '-'), function(pair) paste(rev(pair), collapse='-'))
    paired_combos <- c(paired_combos, reverse_pairs)
                              
    # get all remaining novel pairing combos for all rats available for pairing
    avail_int <- sapply(strsplit(avail_ids, "_"), tail, 1)
    avail_sex <- sapply(avail_int, function(x) ifelse(as.integer(x) > 8, 'F', 'M'))

    avail_f <- avail_ids[which(avail_sex=='F')]
    avail_m <- avail_ids[which(avail_sex=='M')]

    avail_f_fam <- unique(sapply(avail_f, function(x) sub('.*_B(\\d+)_.*', '\\1', x)))
    avail_m_fam <- unique(sapply(avail_m, function(x) sub('.*_B(\\d+)_.*', '\\1', x)))

    avail_combos <- expand.grid(avail_f_fam, avail_m_fam)
    colnames(avail_combos) <- c('dam_fam','sire_fam')                        
    avail_combos$combo <- paste0(avail_combos[,1], '-', avail_combos[,2])
    avail_combos <- avail_combos[!avail_combos$combo %in% paired_combos,]

    # kinship of all possible combos still available to pair
    k_avail <- kinship[kinship$combo %in% avail_combos$combo,]
    k_avail <- k_avail[k_avail$dam_fam != k_avail$sire_fam,]
    k_avail <- k_avail[,3:(ncol(k_avail)-1)] # drop specific animal IDs

    if (k_all) {
        outfile <- file.path(outdir, paste0('kinship_all_avail_fams_',datestamp,'.csv'))
        write.csv(k_avail, outfile, row.names=F, quote=F, na='')
        cat('\nKinship for all available pairings written to', outfile, '\n')
    }

    # identify best replacements for each ID to be replaced
    if (!is.null(ids_to_replace)) {
        
        replacements_list <- list()
        ids_int <- sapply(ids_to_replace, function(x) tail(strsplit(x[1], "_")[[1]], 1))
        ids_sex <- sapply(ids_int, function(x) ifelse(as.integer(x) > 8, 'F', 'M'))
        f_idx <- which(ids_sex=='F')
        m_idx <- which(ids_sex=='M')
        f_to_replace <- ids_to_replace[f_idx]
        m_to_replace <- ids_to_replace[m_idx]
        
        # get the ids and families of mates proposed for the IDs that need replacement
        lonely_dams <- pairs[pairs$sire_animalid %in% m_to_replace,]$dam_animalid 
        lonely_sires <- pairs[pairs$dam_animalid %in% f_to_replace,]$sire_animalid 
        lonely_f_fams <- sapply(lonely_dams, function(x) sub('.*_B(\\d+)_.*', '\\1', x))
        lonely_m_fams <- sapply(lonely_sires, function(x) sub('.*_B(\\d+)_.*', '\\1', x))
        fams_to_replace <- sapply(ids_to_replace, function(x) sub('.*_B(\\d+)_.*', '\\1', x))
        f_fams_to_replace <- fams_to_replace[f_idx]
        m_fams_to_replace <- fams_to_replace[m_idx]

        # designate families to replace as 'available'
        repl_avail_f_fam <- unique(c(avail_f_fam, f_fams_to_replace, lonely_f_fams))
        repl_avail_m_fam <- unique(c(avail_m_fam, m_fams_to_replace, lonely_m_fams))

        repl_avail_combos <- expand.grid(repl_avail_f_fam, repl_avail_m_fam)
        colnames(repl_avail_combos) <- c('dam_fam','sire_fam')
        repl_avail_combos$combo <- paste0(repl_avail_combos$dam_fam, '-', repl_avail_combos$sire_fam)
        repl_k_avail <- kinship[kinship$combo %in% repl_avail_combos$combo,]
        repl_k_avail <- repl_k_avail[!repl_k_avail$combo %in% paired_combos,]
                                  
        for (id in ids_to_replace) {

            id_int <- tail(strsplit(id[1], "_")[[1]], 1)
            id_sex <- ifelse(as.integer(id_int) > 8, 'F', 'M')
            id_fam <- sub('.*_B(\\d+)_.*', '\\1', id)
            id_col <- ifelse(id_sex=='F', 'dam_animalid', 'sire_animalid')
            
            # identify the proposed mate for the ID that needs replacing
            mate_col <- setdiff(c('dam_animalid','sire_animalid'), id_col)
            mate_id <- pairs[pairs[[id_col]]==id,][[mate_col]]
            mate_sex <- setdiff(c('M','F'), id_sex)
            mate_fam <- sub('.*_B(\\d+)_.*', '\\1', mate_id)

            # get the list of IDs that can replace the input ID
            k_id_col <- ifelse(id_sex=='F', 'dam_fam', 'sire_fam')
            k_mate_col <- setdiff(c('sire_fam','dam_fam'), k_id_col)
            k_replacements <- repl_k_avail[repl_k_avail[[k_mate_col]] == mate_fam,]
            k_replacements <- k_replacements[order(k_replacements$kinship),]  
            n_rows <- nrow(k_replacements)

            if (n_rows > n_best) {
                k_replacements <- k_replacements[1:n_best,]
                n_rows <- n_best
            }
            # if no replacements are available without repeating an already-executed pairing,
            # append the ID of the mate to ids_to_pair
            if (n_rows == 0) {
                cat('Could not find a novel replacement for', id, 'to be paired with', mate_id, '\n')
                cat('Consider new pairings as an alternative \n')
            }

            replacements <- data.frame(
                id_to_pair = rep(mate_id, n_rows),
                id_to_replace = rep(id, n_rows),
                use_fam = k_replacements[[k_id_col]],
                kinship = k_replacements$kinship,
                paired_with = rep('NONE', n_rows),
                breederpair = rep('NONE', n_rows),
                comments = rep(NA, n_rows),
                dif_fam = rep(1, n_rows)) # column only used for sorting
            
            # add a first row using the same family as the replacement ID
            rep_same_fam <- data.frame(
                id_to_pair = mate_id,
                id_to_replace = id,
                use_fam = id_fam,
                kinship = kinship[kinship[k_id_col]==id_fam & kinship[k_mate_col]==mate_fam,]$kinship,
                paired_with = 'NONE',
                breederpair = 'NONE',
                comments = NA,
                dif_fam = 0) # column only used for sorting

            replacements <- rbind(rep_same_fam, replacements)
            replacements_list[[id]] <- replacements

        } # end of ID loop  

        # combine replacement lists for all IDs that need replacing 
        replacements <- do.call(rbind, replacements_list)
        replacements <- replacements[with(replacements, order(id_to_replace, dif_fam, kinship)),]
        replacements <- replacements[,1:(ncol(replacements)-1)] #drop dif_fam column
        out$replacements <- replacements
        outfile <- file.path(outdir, paste0('replacement_pairs_n',n_best,'_',datestamp,'.csv'))
        write.csv(replacements, outfile, row.names=F, quote=F, na='')
        cat('\nBest replacement pairings written to', outfile, '\n')
        
    } # end of ids_to_replace

    # identify best alternative pairings for IDs in need of a new mate pair
    if (!is.null(ids_to_pair)) {

        alt_pair_list <- list()
        ids_int <- sapply(ids_to_pair, function(x) tail(strsplit(x[1], "_")[[1]], 1))
        ids_sex <- sapply(ids_int, function(x) ifelse(as.integer(x) > 8, 'F', 'M'))
        f_idx <- which(ids_sex=='F')
        m_idx <- which(ids_sex=='M')

        fams_to_pair <- sapply(ids_to_pair, function(x) sub('.*_B(\\d+)_.*', '\\1', x))
        f_fams_to_pair <- fams_to_pair[f_idx]
        m_fams_to_pair <- fams_to_pair[m_idx]

        avail_f_fams <- sapply(avail_f, function(x) sub('.*_B(\\d+)_.*', '\\1', x))
        avail_m_fams <- sapply(avail_m, function(x) sub('.*_B(\\d+)_.*', '\\1', x))

        alt_avail_f_fams <- unique(c(avail_f_fams, f_fams_to_pair))
        alt_avail_m_fams <- unique(c(avail_m_fams, m_fams_to_pair))

        alt_avail_combos <- (expand.grid(alt_avail_f_fams, alt_avail_m_fams))
        colnames(alt_avail_combos) <- c('dam_fam','sire_fam')
        alt_avail_combos$combo <- paste0(alt_avail_combos[,1], '-', alt_avail_combos[,2])

        # kinship of all possible combos still available to pair
        alt_k_avail <- kinship[kinship$combo %in% alt_avail_combos$combo,]
        alt_k_avail <- alt_k_avail[!alt_k_avail$combo %in% paired_combos,]

        for (id in ids_to_pair) {

            id_int <- tail(strsplit(id[1], "_")[[1]], 1)
            id_sex <- ifelse(as.integer(id_int) > 8, 'F', 'M')
            id_fam <- sub('.*_B(\\d+)_.*', '\\1', id)
            id_col <- ifelse(id_sex=='F', 'dam_animalid', 'sire_animalid')
            mate_col <- setdiff(c('dam_animalid', 'sire_animalid'), id_col)
            mate_id <- pairs[pairs[id_col]==id,mate_col]
            mate_fam <- sub('.*_B(\\d+)_.*', '\\1', mate_id)

            # get the list of IDs that can be paired with the input ID
            k_id_col <- ifelse(id_sex=='F', 'dam_fam', 'sire_fam')
            k_mate_col <- setdiff(c('dam_fam','sire_fam'), k_id_col)
            k_alt_pairs <- alt_k_avail[alt_k_avail[[k_id_col]]==id_fam,]
            k_alt_pairs <- k_alt_pairs[order(k_alt_pairs$kinship),]
            n_rows <- nrow(k_alt_pairs)

            if (n_rows > n_best) {
                k_alt_pairs <- k_alt_pairs[1:n_best,]
                n_rows <- n_best
            }

            alt_pairs <- data.frame(
                animalid = rep(id, n_rows),
                use_fam = k_alt_pairs[[k_mate_col]],
                kinship = k_alt_pairs$kinship,
                paired_with = rep('NONE', n_rows),
                breederpair = rep('NONE', n_rows),
                comments = rep(NA, n_rows),
                dif_fam = rep(1, n_rows)) # column used only for sorting
   
            # add a first row using the same family as the replacement ID
            alt_same_fam <- data.frame(
                animalid = id,
                use_fam = mate_fam,
                kinship = kinship[kinship[k_id_col]==id_fam & kinship[k_mate_col]==mate_fam, 'kinship'],
                paired_with = 'NONE',
                breederpair = 'NONE',
                comments = NA,
                dif_fam = 0) # column used only for sorting


            alt_pairs <- rbind(alt_same_fam, alt_pairs)
            alt_pair_list[[id]] <- alt_pairs
            
        } # end of ID loop  

        # combine replacement lists for all IDs that need replacing 
        alt_pairs <- do.call(rbind, alt_pair_list)
        alt_pairs <- alt_pairs[with(alt_pairs, order(animalid, dif_fam, kinship)),]
        alt_pairs <- alt_pairs[,1:(ncol(alt_pairs)-1)] #drop dif_fam column
        out$alt_pairs <- alt_pairs
        outfile <- file.path(outdir, paste0('alternative_pairs_n',n_best,'_',datestamp,'.csv'))
        write.csv(alt_pairs, outfile, row.names=F, quote=F, na='')
        cat('\nBest alternative pairings written to', outfile, '\n')

    } # end of ids_to_pair
    
    # identify the best pairings among available IDs
    if (is.null(ids_to_replace) && is.null(ids_to_pair)) {
        
        # find the n_best lowest-kinship pairings that are available
        kin_best <- k_avail[order(k_avail$kinship),]
        new_pairs <- data.frame()
        used_dams <- c()
        used_sires <- c()
    
        for (i in 1:nrow(kin_best)) {
            row <- kin_best[i, , drop=F]
            if (!row$dam_fam %in% used_dams && !row$sire_fam %in% used_sires) {
                new_pairs <- rbind(new_pairs, row)
                used_dams <- c(used_dams, row$dam_fam)
                used_sires <- c(used_sires, row$sire_fam)
                if (nrow(new_pairs) == n_best) break
            }
        }
        new_pairs$paired_dam <- 'NONE'
        new_pairs$paired_sire <- 'NONE'
        new_pairs$breederpair <- 'NONE'
        new_pairs$comments <- NA

        # save best new pairs
        out$new_pairs <- new_pairs
        outfile <- file.path(outdir, paste0('new_pairs_n',n_best,'_',datestamp,'.csv'))
        write.csv(new_pairs, outfile, row.names=F, quote=F, na='')
        cat('\nTop', n_best, 'new pairings written to', outfile, '\n')
    }
}

# function to plot a histogram of pairwise kinship coefficients
plot_k_hist <- function(
    kinship, # path to the pairwise kinship or breederpair csv, or a corresponding R dataframe
    pop, # the population name, eg 'hsw' or 'wfu'
    gen, # the generation number
    sample = c('all','breederpairs'), # which sample of rats (ie which file type) to process
    out_dir # output directory path
) {
    library(viridis)

    if (class(kinship) == 'character') {
        kinship <- read.csv(kinship)
    } else if (class(kinship) == 'data.frame') {
        kinship <- kinship
    }
    
    kinship <- as.numeric(kinship$kinship)
    mean_k <- mean(kinship) 
    quantiles <- quantile(kinship, probs = seq(0, 1, 0.25))

    if (sample == 'all') {
        outfile <- file.path(out_dir, paste0(pop, '_gen', gen, '_k_hist_all.png'))
        binsize = 60
    } else if (sample == 'breederpairs') {
        outfile <- file.path(out_dir, paste0(pop, '_gen', gen, '_k_hist_n', length(kinship), '_breeders.png'))
        binsize = 30
    }

    # plot breederpairs kinship
    png(outfile, width=7, height=6, units='in', res=300)
    par(mar=c(4,4,3.6,1), oma=c(0.4,0.4,0,0))
    hist(kinship, binsize, col=viridis(1,0.6,0.4), xlab='kinship (K)', main='')
    if (sample == 'all') {
        title_str <- paste('(all', length(kinship), 'possible pairs)')
    } else if (sample == 'breederpairs') {
        title_str <- paste(paste0('(', length(kinship)), 'breeder pairs)')
    }
    title(main = paste(toupper(pop), paste0('gen', gen), 'pairwise kinship'), line=2.2)
    title(main = title_str, line=1)

    # add mean K and 25% quantile lines
    abline(v=mean_k, lwd=2, col='red')
    for (q in 2:4) {
        abline(v=quantiles[q], lty=3)
        if (q==3) {
            abline(v=quantiles[q], lty=1, lwd=2)
        }
    }

    # add quantiles legend
    legend('topright', legend=c(paste('max:', round(quantiles[5], 3)),
                                paste('75%:', round(quantiles[4], 3)),
                                paste('mean:', round(mean_k, 3)),
                                paste('median:', round(quantiles[3], 3)),
                                paste('25%:', round(quantiles[2], 3)),
                                paste('min:', round(quantiles[1], 3))
                                ), 
        lty=c(0,3,1,1,3,0), col=c(1,1,'red',1,1,1), lwd=c(1,1,2,2,1,1), cex=0.9, bty='n')
    dev.off()
    cat('Kinship histogram saved to', outfile, '\n')
}


plot_k_network <- function(
    kinship, # path to the pairwise kinship or breederpair csv, or a corresponding R dataframe
    pop, # the population name, eg 'hsw' or 'wfu'
    gen, # the generation number
    exchange = FALSE, # whether this is an exchange generation
    sample = c('all','breederpairs'), # which sample of rats (ie which file type) to process
    out_dir # output directory path
) {
    library(viridis)
    library(ggplot2)

    if (class(kinship) == 'character') {
        kinship <- read.csv(kinship)
    } else if (class(kinship) == 'data.frame') {
        kinship <- kinship
    }
    kinship$kinship <- as.numeric(kinship$kinship)
    
    if (sample == 'all') {
        outfile <- file.path(out_dir, paste0(pop, '_gen', gen, '_k_net_all.png'))
        title_str <- paste('(all', nrow(kinship), 'possible pairs)')
        kinship <- kinship[kinship$dam_fam != kinship$sire_fam,]
        lwd <- 0.2

    } else if (sample == 'breederpairs') {
        outfile <- file.path(out_dir, paste0(pop, '_gen', gen, '_k_net_n', nrow(kinship), '_breeders.png'))
        title_str <- paste(paste0('(', nrow(kinship)), 'breeder pairs)')
        lwd <- 2
        if (pop == 'hsw') {
            kinship$dam_fam <- gsub('B','',sapply(kinship$dam_animalid, function(x) {unlist(strsplit(x, '_'))[2]}))
            kinship$sire_fam <- gsub('B','',sapply(kinship$sire_animalid, function(x) {unlist(strsplit(x, '_'))[2]}))
            if (exchange) {
                kinship$dam_fam <- paste0('HSW_', kinship$dam_fam)
                kinship$sire_fam <- sapply(kinship$sire_animalid, function(x) {gsub(".*[0-9]([0-9]{2})$", "\\1", x)})
                kinship$sire_fam <- paste0('WFU_', kinship$sire_fam)
            }
        } else if (pop == 'wfu') {
            kinship$dam_fam <- sapply(kinship$dam_animalid, function(x) {gsub(".*[0-9]([0-9]{2})$", "\\1", x)})
            kinship$sire_fam <- sapply(kinship$sire_animalid, function(x) {gsub(".*[0-9]([0-9]{2})$", "\\1", x)})
            if (exchange) {
                kinship$dam_fam <- paste0('WFU_', kinship$dam_fam)
                kinship$sire_fam <- gsub('B','HSW_',sapply(kinship$sire_animalid, function(x) {unlist(strsplit(x, '_'))[2]}))                
            }
            
        }
    }

    main_title <- paste(toupper(pop), paste0('gen', gen), 'pairwise kinship')
    fam_ids <- sort(unique(c(kinship$dam_fam, kinship$sire_fam)), decreasing=T)
    n_fams <- length(fam_ids)
    min_k <- min(kinship$kinship)
    max_k <- max(kinship$kinship)

    # positions each family around the circle
    fam_pos <- data.frame(
        fam_id = fam_ids,
        # angle = seq(0, 2*pi, length.out = n_fams + 1)[1:n_fams],
        angle = seq(0, 360, length.out = n_fams + 1)[1:n_fams],
        x = cos(seq(0, 2*pi, length.out = n_fams + 1)[1:n_fams]),
        y = sin(seq(0, 2*pi, length.out = n_fams + 1)[1:n_fams])
    )
    fam_pos$angle <- ifelse(fam_pos$angle > 90 & fam_pos$angle <= 270, fam_pos$angle + 180, fam_pos$angle)

    # create connections between fams
    connections <- merge(kinship, fam_pos, by.x = 'dam_fam', by.y = 'fam_id')
    names(connections)[names(connections) %in% c('x','y')] <- c('x1','y1')
    connections <- merge(connections, fam_pos, by.x = 'sire_fam', by.y = 'fam_id')
    names(connections)[names(connections) %in% c('x','y')] <- c('x2','y2')    

    p <- ggplot() +
    # curves for family connections
    geom_curve(
      data = connections,
      aes(x = x1, y = y1, xend = x2, yend = y2, color = kinship, alpha = kinship), 
      curvature = 0, 
      linewidth = lwd) +  
    # points for family positions
    geom_point(
      data = fam_pos,
      aes(x = x, y = y),
      size = 4,
      color = "black") +
    # family labels
    geom_text(
      data = fam_pos,
      aes(x = x * 1.1, y = y * 1.1, label = fam_id, angle = angle), size = 3.5) +
    # color lines by K
    scale_color_gradient(
      low = viridis(1,1,0),
      high = viridis(1,1,0.5),
      name = 'Kinship\nCoefficient') +
    # scale size by K
    scale_size_continuous(range = c(0.5, 2)) +
    # scale opacity by K
    scale_alpha_continuous(range = c(3 * min_k, 3 * max_k), name = 'Kinship\nCoefficient') +
    guides(linewidth = 'none', size = 'none', alpha = 'none') +
    coord_fixed() +
    theme_void() +
    ggtitle(label = main_title, subtitle = title_str) +
    theme(
      legend.position = 'right',
      plot.margin = margin(10, 10, 10, 10),
      plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),  
      plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)) 
    ) 

    png(outfile, width=9, height=9, units='in', res=300)
    # par(mar=c(4,4,3.6,1), oma=c(0.4,0.4,0,0))
    print(p)
    dev.off()
    cat('Kinship network saved to', outfile, '\n')

}


merge_comments <- function(str1, str2) {
    # Initialize output vector with same length as inputs
    result <- vector("character", length(str1))
    
    # Remove commas from both vectors
    str1 <- gsub(',','', str1)
    str2 <- gsub(',','', str2)
    
    # Process each pair of strings
    for(i in seq_along(str1)) {
        s1 <- str1[i]
        s2 <- str2[i]
        
        if(is.na(s1) && is.na(s2)) {
            result[i] <- NA
        } else if(is.na(s1)) {
            result[i] <- s2
        } else if(is.na(s2)) {
            result[i] <- s1
        } else if(s1 == s2) {
            result[i] <- s1
        } else if(startsWith(s2, s1)) {
            # if s1 is prefix of s2, only show the additional part
            additional <- substring(s2, nchar(s1) + 1)
            result[i] <- paste(s1, additional, sep = ' | ')
        } else {
            result[i] <- paste(s1, s2, sep = ' | ')
        }
    }
    
    return(result)
}

# function to finalize a proposed breeder file with edits made in the colony
final_breeder_file <- function(
    proposed_pairs, # csv path or dataframe: the original proposed pairings file
    colony_pairs, # csv path or dataframe: the pairings file executed in the colony
    colony_df, # csv path or dataframe: the current generation's colony dataframe
    pop,
    gen,
    outdir=NULL, 
    alt_pairs = NULL, # csv path or dataframe: the 'alternative pairs' file executed in the colony
    rep_pairs = NULL, # csv path or dataframe: the 'replacement pairs' file executed in the colony
    new_pairs = NULL) # csv path or dataframe: the 'new pairs' file executed in the colony
{
    if (class(proposed_pairs)=='character') {
        proposed_pairs <- read.csv(proposed_pairs)
    }
    if (class(colony_pairs)=='character') {
        colony_pairs <- read.csv(colony_pairs)
    }
    if (class(colony_df)=='character') {
        colony_df <- read.csv(colony_df)
    }
    if (!is.null(alt_pairs)) {
        if (class(alt_pairs)=='character') {
            alt_pairs <- read.csv(alt_pairs)
        }
    }
    if (!is.null(rep_pairs)) {
        if (class(rep_pairs)=='character') {
            rep_pairs <- read.csv(rep_pairs)
        }
    }
    if (!is.null(new_pairs)) {
        if (class(new_pairs)=='character') {
            new_pairs <- read.csv(new_pairs)
        }
    }

    # check that IDs are identical between proposed and executed files
    cols_to_compare <- c('generation','breederpair','kinship','dam_animalid','dam_earpunch',
        'dam_coatcolor','sire_animalid','sire_earpunch','sire_coatcolor')
    keep_cols <- c('generation','breederpair','kinship','dam_rfid','sire_rfid','dam_animalid','sire_animalid')

    # eliminate extra rows from the colony file if they were already added manually
    colony_pairs <- colony_pairs[1:nrow(proposed_pairs),]

    # ensure that no IDs have been changed in the colony file
    crosscheck <- identical(proposed_pairs[,cols_to_compare], colony_pairs[,cols_to_compare])

    if (!crosscheck) {

        cat('Error: IDs in the proposed pairing file and the colony file are not identical! \n')
        cat('Check the following rows for mismatches before proceeding: \n')
        for (col in cols_to_compare) {
            if (!identical(proposed_pairs[[col]], colony_pairs[[col]])) {
                mismatches <- which(proposed_pairs[[col]]!= colony_pairs[[col]])
                cat(col, ':', mismatches, '\n')
                # print(cbind(proposed_pairs[[col]], (colony_pairs[[col]])))
                colname_proposed <- paste0('proposed_pairs$',col)
                colname_colony <- paste0('colony_pairs$',col)
                df <- data.frame(
                    prop = proposed_pairs[[col]],
                    col = colony_pairs[[col]]
                )
                df <- df[mismatches,]
                names(df) <- c(colname_proposed, colname_colony)
                print(df)
            }
        }
    } #else {

        pairs <- proposed_pairs[,keep_cols]
        pairs$dam_accessid <- sapply(pairs$dam_animalid, animalid_to_accessid)
        pairs$sire_accessid <- sapply(pairs$sire_animalid, animalid_to_accessid)

        # if only comments differ between files, incorporate new comments
        if (identical(proposed_pairs$comments, colony_pairs$comments)) {
            pairs$comments <- proposed_pairs$comments
        } else {
            pairs$comments <- merge_comments(proposed_pairs$comments, colony_pairs$comments)
        }
    #}

    # incorporate alternative pairings into the file
    if (!is.null(alt_pairs)) {

        alt_pairs <- alt_pairs[alt_pairs$paired_with != 'NONE',]
        alt_pairs$dam_id <- NA
        alt_pairs$sire_id <- NA
        alt_pairs$dam_rfid <- NA
        alt_pairs$sire_rfid <- NA

        add_alt_pairs <- data.frame(matrix(ncol = length(keep_cols), nrow = nrow(alt_pairs)))
        colnames(add_alt_pairs) <- keep_cols

        for (i in 1:nrow(alt_pairs)) {
            id <- alt_pairs$animalid[i]
            mate_id <- alt_pairs$paired_with[i]
            id_int <- tail(strsplit(id, "_")[[1]], 1)
            id_sex <- ifelse(as.integer(id_int) > 8, 'F', 'M')
            if (id_sex=='F') {
                alt_pairs$dam_id[i] <- id
                alt_pairs$sire_id[i] <- mate_id
            } else {
                alt_pairs$dam_id[i] <- mate_id
                alt_pairs$sire_id[i] <- id
            }
        }

        # save alternative pairings to the new df
        add_alt_pairs$generation <- proposed_pairs$generation[1]
        add_alt_pairs$breederpair <- alt_pairs$breederpair
        add_alt_pairs$kinship <- format(round(alt_pairs$kinship, 4), nsmall=4)
        add_alt_pairs$dam_animalid <- alt_pairs$dam_id
        add_alt_pairs$dam_rfid <- as.character(colony_df$rfid[match(alt_pairs$dam_id, colony_df$animalid)])
        add_alt_pairs$sire_animalid <- alt_pairs$sire_id
        add_alt_pairs$sire_rfid <- as.character(colony_df$rfid[match(alt_pairs$sire_id, colony_df$animalid)])
        add_alt_pairs$dam_accessid <- sapply(alt_pairs$dam_id, animalid_to_accessid)
        add_alt_pairs$sire_accessid <- sapply(alt_pairs$sire_id, animalid_to_accessid)
        default_comment <- 'alternative pairing assigned using breedHS'
        add_alt_pairs$comments <- ifelse(is.na(alt_pairs$comments),
                                        default_comment,
                                        paste(default_comment, alt_pairs$comments, sep = ' | '))
        
        # append alternative pairings to the breeder file
        pairs <- rbind(pairs, add_alt_pairs)
        
    } # end of alt_pairs

    # incorporate replacement pairings into the file
    if (!is.null(rep_pairs)) {

        rep_pairs <- rep_pairs[rep_pairs$paired_with != 'NONE',]
        rep_pairs$dam_id <- NA
        rep_pairs$sire_id <- NA
        rep_pairs$dam_rfid <- NA
        rep_pairs$sire_rfid <- NA
        rep_pairs$breedhs_comments <- NA
        
        add_rep_pairs <- data.frame(matrix(ncol = length(keep_cols), nrow = nrow(rep_pairs)))
        colnames(add_rep_pairs) <- keep_cols

        for (i in 1:nrow(rep_pairs)) {

            bp <- rep_pairs$breederpair[i]
            kinship <- format(round(rep_pairs$kinship[i], 4), nsmall=4)
            colony_comment <- rep_pairs$comments[i]
            id_to_replace <- rep_pairs$id_to_replace[i]
            id <- rep_pairs$id_to_pair[i]
            mate_id <- rep_pairs$paired_with[i]
            id_int <- tail(strsplit(id, "_")[[1]], 1)
            id_sex <- ifelse(as.integer(id_int) > 8, 'F', 'M')
            if (id_sex=='F') {
                dam_id <- id
                sire_id <- mate_id
                id_row <- which(pairs$dam_animalid == id)
                rep_id_row <- which(pairs$sire_animalid == id_to_replace)
                rep_breeder <- 'sire'
            } else {
                dam_id <- mate_id
                sire_id <- id
                id_row <- which(pairs$sire_animalid == id)
                rep_id_row <- which(pairs$dam_animalid == id_to_replace)
                rep_breeder <- 'dam'
            }

            # ensure that the correct replacement is being made - all indiv and bp IDs line up
            bp_row <- which(pairs$breederpair == bp)
            if (id_row == rep_id_row & id_row != bp_row) {
                
                cat('Check your breederpair number \n')
                cat('IDs', dam_id, 'and', sire_id, 'were originally proposed to form', pairs$breederpair[id_row], '\n')
                quit(status=1)

            # if IDs align, overwrite old IDs with the new pairing
            } else if (id_row == rep_id_row & id_row == bp_row) {

                dam_rfid <- colony_df$rfid[which(colony_df$animalid == dam_id)]
                dam_accessid <- animalid_to_accessid(dam_id)
                sire_rfid <- colony_df$rfid[which(colony_df$animalid == sire_id)]
                sire_accessid <- animalid_to_accessid(sire_id)
                
                pairs$dam_animalid[id_row] <- dam_id
                pairs$sire_animalid[id_row] <- sire_id
                pairs$dam_accessid[id_row] <- dam_accessid
                pairs$sire_accessid[id_row] <- sire_accessid
                pairs$dam_rfid[id_row] <- dam_rfid
                pairs$sire_rfid[id_row] <- sire_rfid
                pairs$kinship[id_row] <- kinship

                # log all comments made for this pairing
                pairs_comment <- colony_pairs$comments[id_row]
                pairs_comment <- ifelse(pairs_comment == 'breederpair assigned using breedHS', NA, pairs_comment)
                pairs_comment <- gsub('breederpair assigned using breedHS', '', pairs_comment)
                breedhs_comment <- paste('original', rep_breeder, id_to_replace, 'replaced with', mate_id, 'using breedHS')
                merged_comment <- merge_comments(pairs_comment, breedhs_comment)
                merged_comment <- merge_comments(merged_comment, colony_comment)
                merged_comment <- gsub('NA', '', merged_comment)
                merged_comment <- gsub(' \\| \\| ', ' | ', merged_comment)  # replace double separators
                merged_comment <- gsub('^\\s*\\|\\s*', '', merged_comment)  # reemove leading separator
                merged_comment <- gsub('\\s*\\|\\s*$', '', merged_comment)  # remove trailing separator
                pairs$comments[id_row] <-  merged_comment
                
            } else {
                cat('General error with replacement IDs Check your individual and breederpair IDs \n')
                quit(status=1)
            }
        
        } # end of rows loop                       

    } # end of rep_pairs

    # incorporate new pairings into the file
    if (!is.null(new_pairs)) {

        new_pairs <- new_pairs[new_pairs$breederpair != 'NONE',]
        new_pairs$dam_id <- NA
        new_pairs$sire_id <- NA
        new_pairs$dam_rfid <- NA
        new_pairs$sire_rfid <- NA

        add_new_pairs <- data.frame(matrix(ncol = length(keep_cols), nrow = nrow(new_pairs)))
        colnames(add_new_pairs) <- keep_cols

        # save novel pairings to the new df
        add_new_pairs$generation <- proposed_pairs$generation[1]
        add_new_pairs$breederpair <- new_pairs$breederpair
        add_new_pairs$kinship <- format(round(new_pairs$kinship, 4), nsmall=4)
        add_new_pairs$dam_animalid <- new_pairs$paired_dam
        add_new_pairs$dam_rfid <- as.character(colony_df$rfid[match(new_pairs$paired_dam, colony_df$animalid)])
        add_new_pairs$sire_animalid <- new_pairs$paired_sire
        add_new_pairs$sire_rfid <- as.character(colony_df$rfid[match(new_pairs$paired_sire, colony_df$animalid)])
        add_new_pairs$dam_accessid <- sapply(new_pairs$paired_dam, animalid_to_accessid)
        add_new_pairs$sire_accessid <- sapply(new_pairs$paired_sire, animalid_to_accessid)
        default_comment <- 'novel pairing assigned using breedHS'
        add_new_pairs$comments <- merge_comments(default_comment, new_pairs$comments)
        
        # append novel pairings to the breeder file
        pairs <- rbind(pairs, add_new_pairs)

    } # end of new_pairs

    if (!is.null(outdir)) {
        datestamp <- format(Sys.time(),'%Y%m%d')
        outfile <- paste0(pop, '_gen', gen, '_', gen+1, '_breeders_final_', datestamp, '.csv')
        write.csv(pairs, file.path(outdir, outfile), row.names=F, quote=F, na='')
        cat('Saved final breeders file to', file.path(outdir, outfile), '\n')
    }
    return(pairs)
}


concat_peds <- function(
    directory,
    stem,
    outstem=NULL,
    return_df=FALSE
) {
    ped_files <- list.files(directory, full.names=T, pattern=stem)
    all_peds <- list()
    gens <- c()

    # for (file in ped_files) {
    for (i in 1:length(ped_files)) {
        ped <- read.csv(ped_files[[i]])

        if ('Generation' %in% colnames(ped)) {
            min_gen <- gsub('00$','', min(ped['Generation']))
            max_gen <- gsub('00$','', max(ped['Generation']))
        } else {
            min_gen <- min(ped['generation'])
            max_gen <- max(ped['generation'])
        }
        gens <- c(gens, min_gen, max_gen)
        all_peds[[i]] <- ped
    }

    min_gen <- min(as.numeric(gens))
    max_gen <- max(as.numeric(gens))
    min_nchar <- nchar(min_gen)
    max_nchar <- nchar(max_gen)    

    full_ped <- do.call(rbind, all_peds)

    if ('SWID' %in% colnames(full_ped)) {
        full_ped <- full_ped[order(full_ped$Generation, full_ped$SWID),]
        min_gen <- ifelse(min_nchar==1, paste0('0',min_gen), min_gen)
        max_gen <- ifelse(max_nchar==1, paste0('0',max_gen), max_gen)

    } else {
        full_ped <- full_ped[order(full_ped$generation, full_ped$animalid),]
        min_gen <- ifelse(min_nchar==2, paste0('0',min_gen), min_gen)
        max_gen <- ifelse(max_nchar==2, paste0('0',max_gen), max_gen)

    }

    if (!is.null(outstem)) {
        outfile <- file.path(directory, paste0(outstem, min_gen, '_', max_gen, '.csv'))
        write.csv(full_ped, outfile, row.names=F, quote=F, na='')
        cat('Concatenated pedigree saved to', outfile, '\n')
    }
    if (return_df) {
        return(full_ped)
    }
}


# create a map for HSW IDs used in a merged pedigree
map_merged_ids_hsw <- function(
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

    use_cols <- c('id','rfid','animalid')
    ped1 <- ped1[,use_cols]
    ped2 <- ped2[,use_cols]
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


create_hsw_shipping_sheet <- function(
    pairs,        # R dataframe or csv path, as output by select.breeders, dam/sire must be animal IDs
    assignments,  # HSW colony assignments for the current generation of rats to be sent
    n_ship,       # the total number of rats to be shipped (must be <= nrow(pairs)
    outdir
) {
    if (class(pairs) == 'data.frame') {
        pairs <- pairs
    } else if (class(pairs) == 'character') {
        pairs <- read.csv(pairs, na.str=c('','NA','NaN','nan'))
    } else if (sum(!names(pairs) %in% c('pairs','file')) == 0) {
        pairs <- pairs$pairs
    }
    if (file.exists(assignments)) {    
        assignments <- read.csv(assignments)
    }
    gen <- assignments$generation[1]
    breeders <- assignments[assignments$assignment == 'hsw_breeders' & assignments$sex=='M',]
    breeders$cage <- NA
    breeders$usage <- 'breeder'
    colnames(breeders)[which(colnames(breeders)=='generation')] <- 'hsw_generation'

    ss <- data.frame(
        cage = NA,
        usage = 'breeder',
        hsw_generation = gen,
        animalid = pairs$sire
    )

    # get metadata for paired IDs
    ss$accessid <- sapply(ss$animalid, function(x) {animalid_to_accessid(x)} )
    ss$rfid <- as.character(sapply(ss$animalid, function(x) { breeders$rfid[breeders$animalid == x] }))
    ss$sex <- sapply(ss$animalid, function(x) { breeders$sex[breeders$animalid == x] })
    ss$coatcolor <- sapply(ss$animalid, function(x) { breeders$coatcolor[breeders$animalid == x] })
    ss$earpunch <- sapply(ss$animalid, function(x) { breeders$earpunch[breeders$animalid == x] })
    ss$dob <- sapply(ss$animalid, function(x) { breeders$dob[breeders$animalid == x] })
    ss$dow <- sapply(ss$animalid, function(x) { breeders$dob[breeders$animalid == x] })
    ss$breederpair <- sapply(ss$animalid, function(x) { breeders$breederpair[breeders$animalid == x] })
    ss$dam <- sapply(ss$animalid, function(x) { breeders$dam[breeders$animalid == x] })
    ss$sire <- sapply(ss$animalid, function(x) { breeders$sire[breeders$animalid == x] })

    # get all leftover breeders available to send as extras
    n_extras <- n_ship - nrow(pairs)
    extra_rats <- breeders[!breeders$animalid %in% ss$animalid, colnames(ss)]
    extra_rats$usage <- 'extra'

    # sample n_extras rats (one per breederpair) to include in the shipping sheet
    keep_fams <- sample(unique(extra_rats$breederpair), n_extras)
    send_extras <- data.frame()
    for (fam in keep_fams) {
        fam_males <- extra_rats[extra_rats$breederpair == fam, ]      
        send_male <- fam_males[sample(nrow(fam_males), 1), ]
        send_extras <- rbind(send_extras, send_male)
    }
    ss <- rbind(ss, send_extras)

    outfile <- file.path(outdir, paste0('hsw_gen', gen, '_breeders_for_wfu.csv'))
    write.csv(ss, outfile, row.names=F, quote=F, na='')
    cat('HSW shipping sheet saved to', outfile, '\n')
    return(list(shipping_sheet = ss, file = outfile))
}

final_hsw_breeders_to_raw_ped <- function(
    pairs,   # path or dataframe - final breeders file
    colony_df, # path or dataframe - the current colony dataframe
    stem,
    outdir # where to save the pedigree file 
) {
    if (class(pairs)=='character') {
        pairs <- read.csv(pairs)
    }
    if (class(colony_df)=='character') {
        colony_df <- read.csv(colony_df)
    }

    paired_rfids <- c(pairs$dam_rfid, pairs$sire_rfid)
    ped_cols <- c('generation','rfid','animalid','accessid','sex','coatcolor','earpunch','dam','sire','comments')
    ped <- colony_df[colony_df$rfid %in% paired_rfids, ped_cols]
    ped <- ped[order(ped$animalid),]
    gen <- ped$generation[1]

    outfile <- file.path(outdir, paste0(stem, gen, '.csv'))
    write.csv(ped, outfile, row.names=F, quote=F, na='')
    cat('Raw HS West pedigree saved to', outfile, '\n')
    
    return(ped)
}