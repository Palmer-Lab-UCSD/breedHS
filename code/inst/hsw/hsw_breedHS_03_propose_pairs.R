library(argparse)

# get arguments
all_args <- commandArgs(trailingOnly = T)

print('all_args:')
print(all_args)

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

# set up variables, paths, directories
if (!exists('avail_ids')) avail_ids <- NULL
if (dir.exists(avail_ids)) avail_ids <- NULL  # set NULL if no file is provided
if (is.character(avail_ids) && file.exists(avail_ids)) avail_ids <- readLines(avail_ids)
if (!exists('proposed_breeders_file')) proposed_breeders_file <- NULL
hsw_merged_stem <- 'hsw_merged'

breedhs_outdir <- file.path(proj_dir, 'breedhs_out')
counts_outdir <- file.path(breedhs_outdir, 'bp_counts')
altpair_outdir <- file.path(proj_dir, 'alt_pairings')
kinship_outdir <- file.path(proj_dir, 'kinship')
colony_dir <- file.path(proj_dir, 'colony_edits')
dir.create(breedhs_outdir, showWarnings=F)
dir.create(counts_outdir, showWarnings=F)
dir.create(altpair_outdir, showWarnings=F)
dir.create(kinship_outdir, showWarnings=F)

hsw_dir <- file.path(peds_dir, 'breedail', 'hsw')
hsw_raw_dir <- file.path(peds_dir, 'raw', 'hsw')
wfu_dir <- file.path(peds_dir, 'breedail', 'wfu')
wfu_raw_dir <- file.path(peds_dir, 'raw', 'wfu')
hsw_merged_dir <- file.path(peds_dir, 'merged', 'hsw', paste0('gen', hsw_last_gen))
dir.create(hsw_merged_dir, recursive=T, showWarnings=F)
current_ped <- file.path(hsw_dir, paste0(hsw_stem, hsw_last_gen, '.csv'))


# merge the WFU pedigree into the HSW pedigree
printout('Merging pedigrees')

merged_ped <- merge.pedigrees(
    merge_into = 'hsw',
    ped_map = ped_map, 
    ex_wfu_hsw = wfu_to_hsw,
    ex_hsw_wfu = hsw_to_wfu,
    dir_wfu = wfu_dir,
    stem_wfu = wfu_stem,
    first_gen_wfu = wfu_first_gen,
    last_gen_wfu = wfu_last_gen,
    dir_hsw = hsw_dir,
    stem_hsw = hsw_stem,
    first_gen_hsw = hsw_first_gen,
    last_gen_hsw = hsw_last_gen,
    as_df = T,
    out_dir = hsw_merged_dir,
    out_stem = hsw_merged_stem,
    verbose = FALSE) # set to TRUE for troubleshooting

printout('Checking the merged pedigree for errors')
find.ped.errors(
    hsw_first_gen,
    hsw_last_gen,
    data_dir = hsw_merged_dir,
    file_stem = hsw_merged_stem,
    write_file = FALSE,
    return_ids = FALSE,
    print_ids = FALSE)

printout('Translating the merged pedigree')
translated_ped <- translate.merged.ids(
    input = merged_ped,    
    id_map = merged_ped$id.map, 
    cols = c('id','sire','dam'),   
    from = rep('merged_id', 3), 
    to = rep('animalid', 3)) 

translated_ped <- translate.merged.ids(
    input = merged_ped,    
    id_map = merged_ped$id.map, 
    cols = c('id','sire','dam'),   
    from = rep('merged_id', 3), 
    to = rep('accessid', 3)) 

# estimate (or read in) kinship, copy files to the final output directory
skip_kinship <- Sys.getenv('skip_k') == 'true'

if (!skip_kinship) {

    printout('Estimating kinship')
    k <- current.kinship(
        df = current_ped, 
        first_gen = hsw_first_gen,
        last_gen = hsw_last_gen,
        data_dir = hsw_merged_dir,
        file_stem = hsw_merged_stem,
        hsw_map = hsw_id_map,
        wfu_map = wfu_id_map,
        id_map = merged_ped$id.map,
        verbose = FALSE) # set TRUE for troubleshooting

    kfiles <- list.files(hsw_merged_dir, pattern='kinship', full.names=T)
    kfile <- read.csv(kfiles[which(grepl('all_pairings',kfiles))])
    
    copy_results <- sapply(kfiles, function(file) {
        result <- file.copy(from = file, to = kinship_outdir, overwrite = TRUE)
        cat('Copying', basename(file), 'into', kinship_outdir, '\n')
    })

} else {

    printout(paste('Reading kinship from file:', kinship_file))
    kfile <- read.csv(kinship_file)
}

# plot pairwise kinship for the entire generation
plot_k_hist(
    kinship = kfile,
    pop = pop,
    gen = hsw_last_gen,
    sample = 'all',
    out_dir = kinship_outdir)

# identify (or read in) proposed breeder pairs
if (!is.null(proposed_pairs_file) && file.exists(proposed_pairs_file)) {
    printout(paste('Reading pairings from file:', proposed_pairs_file))
    breedpairs <- read.csv(proposed_pairs_file)
} else {

    # identify best breeder pairs
    printout('Pairing breeders')

    breedpairs <- select.breeders(
        first_gen = hsw_first_gen,
        last_gen = hsw_last_gen,
        data_dir = hsw_merged_dir,
        out_dir = breedhs_outdir,
        file_stem = hsw_merged_stem,
        one_per_sibship=TRUE,
        verbose = FALSE) # set TRUE for troubleshooting

    # translate breederpair file
    printout('Translating the pairing file')

    breedpairs_accessid <- translate.merged.ids(
        input = breedpairs,    
        id_map = merged_ped$id.map, 
        cols = c('dam','sire'),   
        from = rep('merged_id', 2), 
        to = rep('accessid', 2)) 

    breedpairs_animalid <- translate.merged.ids(
        input = breedpairs,    
        id_map = merged_ped$id.map, 
        cols = c('dam','sire'),   
        from = rep('merged_id', 2), 
        to = rep('animalid', 2)) 

    # create the final breederpair file
    printout('Producing final breederpair file')

    breedpairs <- make_hsw_breeder_file(
        pairs = breedpairs_animalid, 
        df = colony_df, 
        wfu_ss = NULL,    
        verbose = FALSE, # set TRUE for troubleshooting
        outdir = breedhs_outdir)$pairs

}


# identify unpaired breeders
bp_counts <- count_breeder_pairs(
    assignments = hsw_assignments,
    breederpairs = breedpairs,
    outdir = counts_outdir)

unpaired_ids <- bp_counts$assigned.but.not.paired
n_paired <- nrow(breedpairs)
all_paired_ids <- c(breedpairs$dam_animalid, breedpairs$sire_animalid)

print('avail_ids1:')
print(avail_ids)
if (nrow(breedpairs) < n_pairs) {

    printout('Proposing additional pairings')

    if (length(unpaired_ids) > 0) cat('unpaired_ids:', unpaired_ids, '\n')
    all_pairs <- add_extra_pairs_to_breeder_file (
        breedpairs = breedpairs, # pairings so far
        kinship = kfile,  # kinship for the whole colony
        colony_df = colony_df,
        avail_ids = avail_ids,   # currently available ids
        unpaired_ids = if(length(unpaired_ids) > 0) unpaired_ids else NULL, # any breeders that still need to be paired
        n_extra = n_pairs - nrow(breedpairs), # make up the difference between proposed pairs and desired pairs
        outdir = breedhs_outdir,
        verbose = TRUE)
    
    breedpairs <- all_pairs$pairs

    # identify unpaired breeders
    printout('Re-checking for unpaired breeders')
    bp_counts <- count_breeder_pairs(
        assignments = hsw_assignments,
        breederpairs = all_pairs$pairs,
        outdir = counts_outdir)

    all_paired_ids <- c(breedpairs$dam_animalid, breedpairs$sire_animalid)
    avail_ids <- setdiff(avail_ids, c(all_paired_ids))
    print('avail_ids2:')
    print(avail_ids)
}

# get best aternative pairings

# printout('Identifying best alternative pairings for all breeders')
# best_alt_pairs(
#     pairs = breedpairs, 
#     kinship = kfile, 
#     avail_ids = avail_ids, 
#     ids_to_replace = NULL,
#     ids_to_pair = all_paired_ids,
#     n_best = 4,
#     verbose = TRUE,
#     outdir = altpair_outdir)

printout('Identifying best replacement pairings among available rats')
best_alt_pairs(
    pairs = breedpairs, 
    kinship = kfile, 
    avail_ids = avail_ids, 
    ids_to_replace = all_paired_ids,
    ids_to_pair = NULL,
    n_best = 4,
    verbose = FALSE,
    outdir = altpair_outdir)

printout('Identifying best new pairings among available rats')
best_alt_pairs(
    pairs = breedpairs, 
    kinship = kfile, 
    avail_ids = avail_ids, 
    ids_to_replace = NULL,
    ids_to_pair = NULL,
    n_best = 6,
    verbose = FALSE,
    outdir = altpair_outdir)


printout('Breeder pairing complete! Pairing file is ready for use in the colony')

