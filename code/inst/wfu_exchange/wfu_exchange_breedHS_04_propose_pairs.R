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
if (!exists('replace_ids')) replace_ids <- NULL

breedhs_outdir <- file.path(proj_dir, 'breedhs_out')
counts_outdir <- file.path(breedhs_outdir, 'bp_counts')
altpair_outdir <- file.path(proj_dir, 'alt_pairings')
kinship_outdir <- file.path(proj_dir, 'kinship')
colony_dir <- file.path(proj_dir, 'colony_edits')

hsw_dir <- file.path(peds_dir, 'breedail', 'hsw')
hsw_raw_dir <- file.path(peds_dir, 'raw', 'hsw')
wfu_dir <- file.path(peds_dir, 'breedail', 'wfu')
wfu_raw_dir <- file.path(peds_dir, 'raw', 'wfu')

# set up output directories depending on the direction of merge
if (pop == 'hsw') {
    first_gen <- hsw_first_gen
    gen <- hsw_last_gen
    merged_stem <- 'hsw_merged'
    merged_dir <- file.path(peds_dir, 'merged', 'hsw', paste0('gen', hsw_last_gen))
    current_raw_ped <- file.path(hsw_raw_dir, paste0(hsw_raw_stem, hsw_last_gen, '.csv'))
    current_ped <- file.path(hsw_dir, paste0(hsw_stem, hsw_last_gen, '.csv'))
    out_dir <- hsw_dir
} else if (pop == 'wfu') {
    first_gen <- wfu_first_gen
    gen <- wfu_last_gen
    merged_stem <- 'wfu_merged'
    merged_dir <- file.path(peds_dir, 'merged', 'wfu', paste0('gen', wfu_last_gen))
    current_raw_ped <- file.path(wfu_raw_dir, paste0(wfu_raw_stem, wfu_last_gen, '.csv'))
    current_ped <- file.path(wfu_dir, paste0(wfu_stem, wfu_last_gen, '.csv'))
    out_dir <- wfu_dir
}
current_dir <- file.path(out_dir, paste0('gen', gen))
dir.create(merged_dir, showWarnings=F)
dir.create(current_dir, showWarnings=F)


# merge the WFU pedigree into the HSW pedigree
printout('Merging pedigrees')
merged_ped <- merge.pedigrees(
    merge_into = ifelse(pop == 'hsw', 2, 1),
    ped_map = ped_map, 
    ex_1_2 = wfu_to_hsw,
    ex_2_1 = hsw_to_wfu,
    dir_1 = wfu_dir,
    stem_1 = wfu_stem,
    first_gen_1 = wfu_first_gen,
    last_gen_1 = wfu_last_gen,
    dir_2 = hsw_dir,
    stem_2 = hsw_stem,
    first_gen_2 = hsw_first_gen,
    last_gen_2 = hsw_last_gen,
    as_df = T,
    out_dir = merged_dir,
    out_stem = merged_stem)

printout('Checking the merged pedigree for errors')
find.ped.errors(
    first_gen,
    gen,
    data_dir = merged_dir,
    file_stem = merged_stem,
    write_file = F,
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
        first_gen = first_gen,
        last_gen = gen,
        data_dir = merged_dir,
        file_stem = merged_stem,
        hsw_map = hsw_id_map,
        wfu_map = wfu_id_map,
        id_map = merged_ped$id.map)

    kfiles <- list.files(merged_dir, pattern='kinship', full.names=T)
    kfile <- read.csv(kfiles[which(grepl('all_pairings',kfiles))])
    
    copy_results <- sapply(kfiles, function(file) {
        result <- file.copy(from = file, to = kinship_outdir, overwrite = TRUE)
        print(paste("Copying", basename(file), "into", kinship_outdir))
    })

} else {

    printout(paste('Reading kinship from file:', kinship_file))
    kfile <- read.csv(kinship_file)
}

# plot pairwise kinship for the entire generation
plot_k_hist(
    kinship = kfile,
    pop = pop,
    gen = gen,
    sample = 'all',
    out_dir = kinship_outdir)

# identify best breeder pairs
printout('Pairing breeders')

breedpairs <- select.breeders(
    first_gen = wfu_first_gen,
    last_gen = wfu_last_gen,
    data_dir = merged_dir,
    out_dir = breedhs_outdir,
    file_stem = merged_stem,
    one_per_sibship=T)

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
printout('Producing proposal breederpair file')

breedpairs <- make_wfu_breeder_file(
    pairs_animalid = breedpairs_animalid, 
    pairs_accessid = breedpairs_accessid, 
    gen = gen,
    wfu_map = wfu_id_map,
    hsw_map = hsw_id_map,
    outdir = breedhs_outdir)



### OPTIONAL: alternative pairs ###

# # identify unpaired breeders
# bp_counts <- count_breeder_pairs(
#     assignments = assignment_file,
#     breederpairs = breedpairs,
#     outdir = counts_outdir)

# all_paired_ids <- c(breedpairs$dam_animalid, breedpairs$sire_animalid)
# unpaired_ids <- bp_counts$assigned_but_not_paired
# replacement_ids <- c(unpaired_ids, readLines(avail_ids))

# printout('Identifying best replacement pairings among available rats')

# replacements <- best_alt_pairs(
#     pairs = breedpairs, 
#     kinship = kfile, 
#     avail_ids = replacement_ids, 
#     ids_to_replace = all_paired_ids,
#     ids_to_pair = NULL,
#     n_best = 4,
#     return_pairs = TRUE,
#     outdir = altpair_outdir)

# printout('Identifying best new pairings among available rats')
# new_pairs <- best_alt_pairs(
#     pairs = breedpairs, 
#     kinship = kfile, 
#     avail_ids = replacement_ids, 
#     ids_to_replace = NULL,
#     ids_to_pair = NULL,
#     n_best = 6,
#     return_pairs = TRUE,
#     outdir = altpair_outdir)


printout('Breeder pairing complete! Pairing file is ready for use in the colony')

