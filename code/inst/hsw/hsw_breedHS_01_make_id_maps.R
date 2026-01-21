library(argparse)

# get arguments
all_args <- commandArgs(trailingOnly = T)

# find the arg_parser_argument
argparse_idx <- which(all_args == '--arg_parser') + 1
arg_parser <- all_args[argparse_idx]

# read in the argument parser
source(arg_parser)

# parse the command-line arguments
args <- breedhs_parser(commandArgs())
list2env(args, envir = .GlobalEnv)

# read in breedHS functions
code_dir <- dirname(utils)
utils <- basename(utils)
setwd(code_dir)
source(utils)

print(args)

# create output directories
dir.create(file.path(outdir, 'hsw'), recursive=T, showWarnings=F)
dir.create(file.path(outdir, 'wfu'), recursive=T, showWarnings=F)

#### initialize ID maps for pedigree samples ####

# read in the raw WFU pedigree
printout('Reading the raw WFU pedigree')
wfu_ped <- read_wfu_raw_ped(wfu_raw_ped)
colnames(wfu_ped) <- gsub('.', '', colnames(wfu_ped), fixed=T)
keep_cols <- c('Transpondernumber', 'SWID', 'IDF51', 'Sex', 'DamID', 'SireID', 'DamSWID', 'SireSWID', 'Generation')
new_colnames <- c('rfid', 'swid', 'accessid', 'sex', 'dam_accessid', 'sire_accessid', 'dam_swid', 'sire_swid', 'generation')
wfu_map <- wfu_ped[,keep_cols]
names(wfu_map) <- new_colnames
wfu_map$generation <- gsub('00$','',wfu_map$generation)
wfu_map$dam_rfid <- sapply(wfu_map[['dam_accessid']], function(x) accessid_to_rfid(x, wfu_map))
wfu_map$sire_rfid <- sapply(wfu_map[['sire_accessid']], function(x) accessid_to_rfid(x, wfu_map))
wfu_col_order <- c('generation', 'rfid', 'swid', 'accessid', 'sex', 
                   'dam_rfid', 'sire_rfid', 'dam_swid', 'sire_swid', 'dam_accessid', 'sire_accessid')
wfu_map <- wfu_map[, wfu_col_order]

# read in the raw HSW pedigree
printout('Reading the raw HSW pedigree')
hsw_map <- read.csv(hsw_raw_ped)
hsw_col_order <- c('generation', 'rfid', 'animalid', 'accessid', 'sex', 
                   'dam_rfid', 'sire_rfid', 'dam_animalid', 'sire_animalid', 'dam_accessid', 'sire_accessid')
hsw_map <- hsw_map[, hsw_col_order]

# read in current HSW assignments, if being used
if (exists('hsw_assignments')) {
    printout('Reading HSW assignments')
    hsw_assignments <- read.csv(hsw_assignments)
    
    if (!exists('prev_colony_df')) {
        stop('Please provide the previous HSW colony dataframe.')
    }
    printout('Reading the previous HSW colony dataframe')
    prev_df <- read.csv(prev_colony_df)

    # add RFIDs, access IDs for dams & sires of assigned HSW breeders
    printout('Adding assigned HSW breeders to HSW ID map')
    hsw_assigned <- hsw_assignments[hsw_assignments$assignment=='hsw_breeders',]
    names(hsw_assigned)[which(names(hsw_assigned)=='dam')] <- 'dam_animalid'
    names(hsw_assigned)[which(names(hsw_assigned)=='sire')] <- 'sire_animalid'

    hsw_assigned$dam_rfid <- sapply(hsw_assigned$dam_animalid, function(x) {
        prev_df$rfid[which(prev_df$animalid==x)]})
    hsw_assigned$sire_rfid <- sapply(hsw_assigned$sire_animalid, function(x) {
        prev_df$rfid[which(prev_df$animalid==x)]})
    hsw_assigned$dam_accessid <- sapply(hsw_assigned$dam_animalid, function(x) {
        prev_df$accessid[which(prev_df$animalid==x)]})
    hsw_assigned$sire_accessid <- sapply(hsw_assigned$sire_animalid, function(x) {
        prev_df$accessid[which(prev_df$animalid==x)]})

    # hsw_assigned$dam_rfid <- as.character(sapply(hsw_assigned$dam_animalid, function(x) {
    #     matches <- prev_df$rfid[which(prev_df$animalid==x)]
    #     if (length(matches) == 0) return(NA)
    #         return(matches[1])}))
    # hsw_assigned$sire_rfid <- as.character(sapply(hsw_assigned$sire_animalid, function(x) {
    #     matches <- prev_df$rfid[which(prev_df$animalid==x)]
    #     if (length(matches) == 0) return(NA)
    #         return(matches[1])}))
    # hsw_assigned$dam_accessid <- as.character(sapply(hsw_assigned$dam_animalid, function(x) {
    #     matches <- prev_df$accessid[which(prev_df$animalid==x)]
    #     if (length(matches) == 0) return(NA)
    #         return(matches[1])}))
    # hsw_assigned$sire_accessid <- as.character(sapply(hsw_assigned$sire_animalid, function(x) {
    #     matches <- prev_df$accessid[which(prev_df$animalid==x)]
    #     if (length(matches) == 0) return(NA)
    #         return(matches[1])}))

    # add new breeders the the HSW id map
    hsw_assigned <- hsw_assigned[,hsw_col_order]
    hsw_map <- rbind(hsw_map, hsw_assigned)

}

# read in WFU shipping sheets
wfu_ss_list <- list()
printout('Reading WFU shipping sheets')
cat('NOTE:', length(all_wfu_ss), 'WFU shipping sheets provided. Ensure this is the correct number in your arguments! \n')
for (ss_file in all_wfu_ss) {

    filename <- basename(as.character(ss_file))
    wfu_gen <- gsub('^0', '', sub('.*wfu_gen(\\d+).*', '\\1', filename))
    hsw_gen <- gsub('^0', '', sub('.*hsw_gen(\\d+).*', '\\1', filename))
    ss <- as.data.frame(read_excel(ss_file, sheet = 'for HSW'))
    # ss <- read_wfu_shipping_sheet(ss_file)
    keep_cols <- c('Transponder ID', 'Animal ID', 'Access ID', 'Sex', 'Dam', 'Sire')
    new_colnames <- c('rfid', 'swid', 'accessid', 'sex', 'dam_accessid', 'sire_accessid')
    ss <- ss[,keep_cols]
    colnames(ss) <- new_colnames
    ss[['generation']] <- hsw_gen
    wfu_ss_list[[filename]] <- ss
    cat('\t', 'Added WFU shipping sheet:', filename, '\n')
}


# read in HSW shipping sheets
printout('Reading HSW shipping sheets')
cat('NOTE:', length(all_hsw_ss), 'HSW shipping sheets provided. Ensure this is the correct number in your arguments! \n')
hsw_ss_list <- list()

for (ss_file in all_hsw_ss) {
    filename <- basename(as.character(ss_file))
    wfu_gen <- gsub('^0', '', sub('.*wfu_gen(\\d+).*', '\\1', filename))
    hsw_gen <- gsub('^0', '', sub('.*hsw_gen(\\d+).*', '\\1', filename))
    ss <- as.data.frame(read_excel(ss_file, sheet = 'for WFU'))
    keep_cols <- c('rfid', 'animalid', 'sex', 'dam', 'sire')
    new_colnames <- c('rfid', 'animalid', 'sex', 'dam_animalid', 'sire_animalid')
    ss <- ss[,keep_cols]
    colnames(ss) <- new_colnames
    ss[['generation']] <- wfu_gen
    hsw_ss_list[[filename]] <- ss
    cat('\t', 'Added HSW shipping sheet:', filename, '\n')

}

# concatenate all shipping sheets
printout('Concatenating shipping sheets')
wfu_shipped <- do.call(rbind, wfu_ss_list)
hsw_shipped <- do.call(rbind, hsw_ss_list)

# format WFU access IDs to remove underscores
printout('Formatting WFU access IDs')
wfu_shipped[['accessid']] <- gsub('_', '', wfu_shipped[['accessid']])
wfu_shipped[['dam_accessid']] <- gsub('_', '', wfu_shipped[['dam_accessid']])
wfu_shipped[['sire_accessid']] <- gsub('_', '', wfu_shipped[['sire_accessid']])

# get alternate IDs for dams/sires of shipped rats
printout('Getting alternate IDs for shipped WFU rats')
wfu_shipped[['dam_swid']] <- sapply(wfu_shipped[['dam_accessid']], function(x) accessid_to_swid(x, wfu_map))
wfu_shipped[['sire_swid']] <- sapply(wfu_shipped[['sire_accessid']], function(x) accessid_to_swid(x, wfu_map))
wfu_shipped[['dam_rfid']] <- sapply(wfu_shipped[['dam_accessid']], function(x) accessid_to_rfid(x, wfu_map))
wfu_shipped[['sire_rfid']] <- sapply(wfu_shipped[['sire_accessid']], function(x) accessid_to_rfid(x, wfu_map))

printout('Getting alternate IDs for shipped HSW rats')
hsw_shipped[['accessid']] <- sapply(hsw_shipped[['animalid']], animalid_to_accessid)
hsw_shipped[['dam_accessid']] <- sapply(hsw_shipped[['dam_animalid']], animalid_to_accessid)
hsw_shipped[['sire_accessid']] <- sapply(hsw_shipped[['sire_animalid']], animalid_to_accessid)
hsw_shipped[['dam_rfid']] <- sapply(hsw_shipped[['dam_accessid']], function(x) accessid_to_rfid(x, hsw_map))
hsw_shipped[['sire_rfid']] <- sapply(hsw_shipped[['sire_accessid']], function(x) accessid_to_rfid(x, hsw_map))

# reorder columns
wfu_shipped <- wfu_shipped[, wfu_col_order]
hsw_shipped <- hsw_shipped[, hsw_col_order]

# remove samples from shipping sheets if they are already incorporated into pedigrees
wfu_shipped <- wfu_shipped[!wfu_shipped$swid %in% hsw_map$animalid,]
hsw_shipped <- hsw_shipped[!hsw_shipped$animalid %in% wfu_map$swid,]


####################
#### WFU ID MAP ####
####################

# concatenate all ID maps, order by generation + SW ID
printout('Concatenating all WFU data')
names(hsw_shipped) <- wfu_col_order
wfu_map <- rbind(wfu_map, hsw_shipped)
wfu_map[['generation']] <- as.numeric(wfu_map[['generation']])
wfu_map <- wfu_map[order(wfu_map[['generation']], wfu_map[['swid']]),]
wfu_gen <- max(wfu_map[['generation']])

# remove rows without access IDs
wfu_map <- wfu_map[!is.na(wfu_map[['accessid']]),]
wfu_mapfile <- file.path(outdir, 'wfu', paste0('wfu_id_map_gen', wfu_gen, '.csv'))


####################
#### HSW ID MAP ####
####################

# concatenate all ID maps, order by generation + animal ID
printout('Concatenating all HSW data')
names(wfu_shipped) <- hsw_col_order
hsw_map <- rbind(hsw_map, wfu_shipped)
hsw_map[['generation']] <- as.numeric(hsw_map[['generation']])
hsw_map <- hsw_map[order(hsw_map[['generation']], hsw_map[['animalid']]),]
hsw_gen <- max(hsw_map[['generation']])

# remove rows without access IDs
hsw_map <- hsw_map[!is.na(hsw_map[['accessid']]),]
hsw_mapfile <- file.path(outdir, 'hsw', paste0('hsw_id_map_gen', hsw_gen, '.csv'))


## save to file ##
write.csv(wfu_map, wfu_mapfile, row.names=F, quote=F, na='')
cat('WFU ID map saved to', wfu_mapfile, '\n')

write.csv(hsw_map, hsw_mapfile, row.names=F, quote=F, na='')
cat('HSW ID map saved to', hsw_mapfile, '\n\n')

