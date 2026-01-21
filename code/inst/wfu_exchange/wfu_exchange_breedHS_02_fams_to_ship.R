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

datestamp <- format(Sys.time(),'%Y%m%d')

# create output directories
breedhs_outdir <- file.path(proj_dir, 'breedhs_out')
counts_outdir <- file.path(breedhs_outdir, 'bp_counts')
altpair_outdir <- file.path(proj_dir, 'alt_pairings')
kinship_outdir <- file.path(proj_dir, 'kinship')
colony_dir <- file.path(proj_dir, 'colony_edits')
dir.create(breedhs_outdir, showWarnings=F)
dir.create(counts_outdir, showWarnings=F)
dir.create(altpair_outdir, showWarnings=F)
dir.create(kinship_outdir, showWarnings=F)
dir.create(colony_dir, showWarnings=F)

# set pedigree directories
hsw_dir <- file.path(peds_dir, 'breedail', 'hsw')
wfu_dir <- file.path(peds_dir, 'breedail', 'wfu')
hsw_raw_dir <- file.path(peds_dir, 'raw', 'hsw')
wfu_raw_dir <- file.path(peds_dir, 'raw', 'wfu')
hsw_merged_dir <- file.path(peds_dir, 'merged', 'hsw')
wfu_merged_dir <- file.path(peds_dir, 'merged', 'wfu')

# set temporary pedigree directories to house peds used for hypothetical pairings
breedhs_tmp_outdir <- file.path(breedhs_outdir, 'hypothetical_exch_pairs')
dir.create(breedhs_tmp_outdir, showWarnings=F)
tmp_ped_dir <- file.path(breedhs_tmp_outdir, 'hypothetical_pedigrees')
hsw_tmpdir <- file.path(tmp_ped_dir, 'breedail', 'hsw')
hsw_raw_tmpdir <- file.path(tmp_ped_dir, 'raw', 'hsw')
hsw_merged_tmpdir <- file.path(tmp_ped_dir, 'merged', 'hsw')
hsw_idmap_tmpdir <- file.path(tmp_ped_dir, 'id_maps','hsw')
wfu_tmpdir <- file.path(tmp_ped_dir, 'breedail', 'wfu')
wfu_raw_tmpdir <- file.path(tmp_ped_dir, 'raw', 'wfu')
wfu_merged_tmpdir <- file.path(tmp_ped_dir, 'merged', 'wfu')
wfu_idmap_tmpdir <- file.path(tmp_ped_dir, 'id_maps','wfu')

dir.create(hsw_tmpdir, recursive=T, showWarnings=F)
dir.create(hsw_raw_tmpdir, recursive=T, showWarnings=F)
dir.create(hsw_merged_tmpdir, recursive=T, showWarnings=F)
dir.create(hsw_idmap_tmpdir, recursive=T, showWarnings=F)
dir.create(wfu_tmpdir, recursive=T, showWarnings=F)
dir.create(wfu_raw_tmpdir, recursive=T, showWarnings=F)
dir.create(wfu_merged_tmpdir, recursive=T, showWarnings=F)
dir.create(wfu_idmap_tmpdir, recursive=T, showWarnings=F)


# set up a temporary HSW ID map using hypothetical rats to send
# (these rats are used as representatives for their families, actual indivs will be assigned later depending on RATTACA predictions)
if (hsw_id_fams_from=='ped') { # if there is no HSW shipping sheet, use the current colony dataframe

    hsw_prev_gen <- hsw_last_gen - 1
    prev_colony_df <- gsub(hsw_last_gen, hsw_prev_gen, colony_df)
    printout('Reading the HSW colony dataframe')
    colony_df <- read.csv(colony_df)
    
    # read in the HSW ID map up through the most recent generation
    prev_hsw_id_map <- gsub(hsw_last_gen, hsw_prev_gen, hsw_id_map)
    printout('Reading the HSW ID map')
    prev_hsw_idmap <- read.csv(prev_hsw_id_map)

    # read in the raw WFU pedigree
    printout('Reading the WFU raw pedigree')
    wfu_ped <- read_wfu_raw_ped(wfu_raw_ped)
    colnames(wfu_ped) <- gsub('.', '', colnames(wfu_ped), fixed=T)
    
    # # initialize the WFU id map
    # keep_cols <- c('Transpondernumber', 'SWID', 'IDF51', 'Sex', 'DamID', 'SireID', 'DamSWID', 'SireSWID', 'Generation')
    # new_colnames <- c('rfid', 'swid', 'accessid', 'sex', 'dam_accessid', 'sire_accessid', 'dam_swid', 'sire_swid', 'generation')
    # wfu_map <- wfu_ped[,keep_cols]
    # names(wfu_map) <- new_colnames
    # wfu_map$generation <- gsub('00$','',wfu_map$generation)
    # wfu_map$dam_rfid <- sapply(wfu_map[['dam_accessid']], function(x) accessid_to_rfid(x, wfu_map))
    # wfu_map$sire_rfid <- sapply(wfu_map[['sire_accessid']], function(x) accessid_to_rfid(x, wfu_map))
    # wfu_col_order <- c('generation', 'rfid', 'swid', 'accessid', 'sex', 
    #                 'dam_rfid', 'sire_rfid', 'dam_swid', 'sire_swid', 'dam_accessid', 'sire_accessid')
    # wfu_map <- wfu_map[, wfu_col_order]

    # sample one male per breederpair from the HSW colony, format like a shipping sheet
    printout('Sampling rats from the HSW colony dataframe')
    tmp_hsw_rats <- one_rat_per_hsw_fam(
        colony_df = colony_df,
        sex = 'M')

    # add columns to the HSW "shipping sheet"
    keep_cols <- c('rfid', 'animalid', 'sex', 'dam', 'sire')
    tmp_hsw_rats <- tmp_hsw_rats[,keep_cols]
    colnames(tmp_hsw_rats) <-c('rfid', 'animalid', 'sex', 'dam_animalid', 'sire_animalid')
    tmp_hsw_rats[['generation']] <- wfu_last_gen
    tmp_hsw_rats[['accessid']] <- sapply(tmp_hsw_rats[['animalid']], animalid_to_accessid)
    tmp_hsw_rats[['dam_accessid']] <- sapply(tmp_hsw_rats[['dam_animalid']], animalid_to_accessid)
    tmp_hsw_rats[['sire_accessid']] <- sapply(tmp_hsw_rats[['sire_animalid']], animalid_to_accessid)
    tmp_hsw_rats[['dam_rfid']] <- sapply(tmp_hsw_rats[['dam_accessid']], function(x) accessid_to_rfid(x, colony_df))
    tmp_hsw_rats[['sire_rfid']] <- sapply(tmp_hsw_rats[['sire_accessid']], function(x) accessid_to_rfid(x, colony_df))

    # add sampled HSW rats into the WFU pedigree
    printout('Adding sampled HSW rats to the temp WFU raw pedigree')
    add_hsw_rats_to_wfu_raw_ped(
        ped = wfu_raw_ped,
        hsw_shipping_sheet = colony_df[colony_df$rfid %in% tmp_hsw_rats$rfid,],
        add_to_gen = wfu_last_gen,
        outdir = wfu_raw_tmpdir)
    wfu_tmp_raw_ped <- list.files(wfu_raw_tmpdir, pattern='complete', full.names=T)

    # add HSW rats to the WFU ID map
    printout('Adding HSW rats to the WFU ID map')
    tmp_hsw_idmap <- tmp_hsw_rats
    names(tmp_hsw_idmap)[which(names(tmp_hsw_idmap)=='dam')] <- 'dam_animalid'
    names(tmp_hsw_idmap)[which(names(tmp_hsw_idmap)=='sire')] <- 'sire_animalid'
    tmp_hsw_idmap$dam_accessid <- sapply(tmp_hsw_idmap$dam_animalid, animalid_to_accessid)
    tmp_hsw_idmap$sire_accessid <- sapply(tmp_hsw_idmap$sire_animalid, animalid_to_accessid)
    tmp_hsw_idmap$dam_rfid <- sapply(tmp_hsw_idmap$dam_accessid, function(x) accessid_to_rfid(x, colony_df))
    tmp_hsw_idmap$sire_rfid <- sapply(tmp_hsw_idmap$sire_accessid, function(x) accessid_to_rfid(x, colony_df))

    hsw_col_order <- c('generation', 'rfid', 'animalid', 'accessid', 'sex', 
                    'dam_rfid', 'sire_rfid', 'dam_animalid', 'sire_animalid', 'dam_accessid', 'sire_accessid')
    tmp_hsw_idmap <- tmp_hsw_idmap[,hsw_col_order]
    hsw_idmap <- rbind(prev_hsw_idmap, tmp_hsw_idmap)

    tmp_hsw_id_map <- gsub('.csv', '_tmp_hypothetical_fams.csv', basename(hsw_id_map))
    tmp_hsw_id_map <- file.path(hsw_idmap_tmpdir, tmp_hsw_id_map)
    write.csv(hsw_idmap, tmp_hsw_id_map, row.names=F, quote=F, na='')
    cat('Updated (temporary) HSW ID map saved to', tmp_hsw_id_map, '\n')

} else if (hsw_id_fams_from=='ss') {
# use a shipping sheet if available

    stop('NEED TO SET UP CODE FOR IDENTIFYING PAIRS FROM AN HSW SHIPPING SHEET')
}


# convert the WFU raw pedigree into breedail format as part of the HSW pedigree
printout('Converting the raw WFU pedigree')
format_wfu_raw_ped(
    ped = wfu_tmp_raw_ped,
    wfu_map = wfu_id_map,
    outdir = wfu_tmpdir
)

# convert the HSW raw pedigree into breedail format
printout('Converting the raw HSW pedigree')
format_hsw_raw_ped(
    ped = hsw_raw_ped,
    wfu_map = wfu_id_map, 
    hsw_map = prev_hsw_id_map, 
    return_df = F, 
    outdir = hsw_tmpdir)

# incorporate WFU generations into the HSW pedigree (prior to founding HSW)
printout('Incorporating WFU generations into the HSW pedigree')
wfu_raw_to_hsw(
    ped = wfu_raw_ped, 
    wfu_map = wfu_id_map,
    outdir = hsw_tmpdir)

# check the WFU pedigree for errors
printout('Checking the WFU pedigree for errors')
wfu_errors <- find.ped.errors(
    first_gen = wfu_first_gen, 
    last_gen = wfu_last_gen,   
    data_dir = wfu_tmpdir,   
    file_stem = wfu_stem, 
    write_file = TRUE,
    return_ids = FALSE,
    print_ids = TRUE)

printout('Checking the HSW pedigree for errors')
hsw_errors <- find.ped.errors(
    first_gen = hsw_first_gen,  
    last_gen = hsw_prev_gen,   
    data_dir = hsw_tmpdir,   
    file_stem = hsw_stem, 
    write_file = TRUE,
    return_ids = FALSE,
    print_ids = TRUE)

# merge the HSW pedigree into the WFU pedigree
printout('Merging the HSW pedigree into the WFU pedigree')
# note: last_gen_2 is hsw_PREV_gen here, since hsw_last_gen was directly appended to the WFU pedigree
# here, all HSW gens before hsw_last_gen will be merged into the previous wfu generation
wfu_merged_ped <- merge.pedigrees(
    merge_into = 1,
    ped_map = ped_map, 
    ex_1_2 = wfu_to_hsw,
    ex_2_1 = hsw_to_wfu,
    dir_1 = wfu_tmpdir,
    stem_1 = wfu_stem,
    first_gen_1 = wfu_first_gen,
    last_gen_1 = wfu_last_gen,
    dir_2 = hsw_tmpdir,
    stem_2 = hsw_stem,
    first_gen_2 = hsw_first_gen,
    last_gen_2 = hsw_prev_gen, 
    as_df = T,
    out_dir = wfu_merged_tmpdir,
    out_stem = 'wfu_merged')
current_ped <- file.path(wfu_merged_tmpdir, paste0('wfu_merged', wfu_last_gen, '.csv'))

printout('Translating the merged pedigree')
translated_ped <- translate.merged.ped(
    ped = wfu_merged_ped,
    id_map = wfu_merged_ped$id.map,
    wfu_map = wfu_id_map,
    hsw_map = tmp_hsw_id_map)
translated_ped <- list.files(wfu_merged_tmpdir, pattern='all_ids', full.names=T)

# estimate (or read in) kinship, copy files to the final output directory
skip_kinship <- Sys.getenv('skip_k') == 'true'

if (!skip_kinship) {

    printout('Estimating kinship')
    k <- current.kinship(
        df = translated_ped, 
        first_gen = wfu_first_gen,
        last_gen = wfu_last_gen,
        data_dir = wfu_merged_tmpdir,
        file_stem = 'wfu_merged',
        hsw_map = tmp_hsw_id_map,
        wfu_map = wfu_id_map,
        id_map = wfu_merged_ped$id.map)

    kfiles <- list.files(wfu_merged_tmpdir, pattern='kinship', full.names=T)
    kfile <- kfiles[which(grepl('all_pairings',kfiles))]
    k <- read.csv(kfile)
    file.copy(from = kfile, to = breedhs_tmp_outdir)

    
} else {

    printout(paste('Reading kinship from file:', kinship_file))
    k <- read.csv(kinship_file)
}


# identify best breeder pairs
printout('Pairing breeders')
breedpairs <- select.breeders(
    first_gen = wfu_first_gen,
    last_gen = wfu_last_gen,
    data_dir = wfu_merged_tmpdir,
    out_dir = breedhs_tmp_outdir,
    file_stem = 'wfu_merged',
    one_per_sibship=T)

# translate breederpair file
printout('Translating the pairing file')

breedpairs_animalid <- translate.merged.ids(
    input = breedpairs,    
    id_map = wfu_merged_ped$id.map, 
    cols = c('dam','sire'),   
    from = rep('merged_id', 2), 
    to = rep('animalid', 2)) 

# save a list of breederpairs to pair with WFU
fams_to_ship <- sapply(breedpairs_animalid$pairs$sire, get_hsw_fam)
fams_file <- file.path(breedhs_outdir, paste0('hsw_gen',hsw_last_gen,'_fams_to_pair_w_wfu_gen',wfu_last_gen,'_',datestamp))
writeLines(sort(fams_to_ship), fams_file)
cat('List of priority HSW families to pair with WFU saved to', fams_file, '\n')


printout('Hypothetical pairings complete. Use the proposed pairs file to assign HSW families to ship to WFU')