library(argparse)

# get arguments
all_args <- commandArgs(trailingOnly = T)

print('all_args:')
print(all_args)

# find the arg_parser_argument
argparse_idx <- which(all_args == '--arg_parser') + 1
print('argparse_idx:')
print(argparse_idx)
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
if (!exists('replace_ids')) replace_ids <- NULL
hsw_merged_stem <- 'hsw_merged'

breedhs_outdir <- file.path(proj_dir, 'breedhs_out')
counts_outdir <- file.path(breedhs_outdir, 'bp_counts')
altpair_outdir <- file.path(proj_dir, 'alt_pairings')
kinship_outdir <- file.path(proj_dir, 'kinship')
colony_dir <- file.path(proj_dir, 'colony_edits')

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
    as_df = TRUE,
    out_dir = hsw_merged_dir,
    out_stem = hsw_merged_stem,
    verbose = TRUE) # set to TRUE for troubleshooting


printout('Checking the merged pedigree for errors')
find.ped.errors(
    hsw_first_gen,
    hsw_last_gen,
    data_dir = hsw_merged_dir,
    file_stem = hsw_merged_stem,
    write_file = T,
    return_ids = F,
    print_ids = F)

printout('Translating the merged pedigree')
translated_ped <- translate.merged.ped(
    ped = merged_ped,
    id_map = merged_ped$id.map,
    wfu_map = wfu_id_map,
    hsw_map = hsw_id_map)

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
        verbose = TRUE) # set TRUE for troubleshooting

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
    pop = 'hsw',
    gen = hsw_last_gen,
    sample = 'all',
    out_dir = kinship_outdir)

# identify (or read in) best breeder pairs
skip_pairing <- Sys.getenv('skip_pairing') == 'true'

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

# create the proposal breederpair file
printout('Producing proposed breederpair file')

breedpairs <- make_hsw_breeder_file(
    pairs = breedpairs_animalid, 
    df = colony_df, 
    wfu_ss = wfu_shipping_sheet,    
    outdir = breedhs_outdir,
    verbose = FALSE)$pairs

### OPTIONAL: alternative pairs ###

# identify unpaired breeders
printout('Checking for unpaired breeders')
bp_counts <- count_breeder_pairs(
    assignments = hsw_assignments,
    breederpairs = breedpairs,
    outdir = counts_outdir)

all_paired_ids <- c(breedpairs$dam_animalid, breedpairs$sire_animalid)
unpaired_ids <- bp_counts$assigned_but_not_paired
if (!is.null(avail_ids)) avail_ids <- readLines(avail_ids)
# for exchange generations only: treat unpaired IDs as "available"
# these are leftovers after pairing w/ WFU rats, they do not need priority over other available rats
avail_ids <- unique(c(unpaired_ids, avail_ids))
unpaired_ids <- NULL

if (nrow(breedpairs) < n_pairs) {

    printout('Proposing additional pairings')

    if (length(unpaired_ids)>0) cat('unpaired_ids:', unpaired_ids, '\n')
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

}

# NOTE: alt pairs using WFU IDs simply gets too complicated/messy
# for WFU exchange generations, simply print out the kinship pairing file
cat('Use the kinship pairs file to identify alternative/replacement pairs if needed \n')

# printout('Identifying best replacement pairings among available rats')
# replacements <- best_alt_pairs(
#     pairs = breedpairs, 
#     kinship = kfile, 
#     avail_ids = avail_ids, 
#     ids_to_replace = all_paired_ids,
#     ids_to_pair = NULL,
#     n_best = 4,
#     return_pairs = TRUE,
#     outdir = altpair_outdir,
#     wfu_map = wfu_id_map,
#     hsw_map = hsw_id_map,
#     verbose = TRUE)

# printout('Identifying best new pairings among available rats')
# new_pairs <- best_alt_pairs(
#     pairs = breedpairs, 
#     kinship = kfile, 
#     avail_ids = avail_ids, 
#     ids_to_replace = NULL,
#     ids_to_pair = NULL,
#     n_best = 12,
#     return_pairs = TRUE,
#     wfu_map = wfu_id_map,
#     hsw_map = hsw_id_map,
#     outdir = altpair_outdir,
#     verbose = TRUE)


printout('Breeder pairing complete! Pairing file is ready for use in the colony')
