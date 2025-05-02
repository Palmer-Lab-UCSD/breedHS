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
altpair_outdir <- file.path(proj_dir, 'alt_pairings')
kinship_outdir <- file.path(proj_dir, 'kinship')
colony_dir <- file.path(proj_dir, 'colony')

hsw_dir <- file.path(peds_dir, 'breedail', 'hsw')
hsw_raw_dir <- file.path(peds_dir, 'raw', 'hsw')
hsw_current_raw_ped <- file.path(hsw_raw_dir, paste0(hsw_raw_stem, hsw_last_gen, '.csv'))
wfu_dir <- file.path(peds_dir, 'breedail', 'wfu')
wfu_raw_dir <- file.path(peds_dir, 'raw', 'wfu')

# set up output directories depending on the direction of merge
if (pop == 'hsw') {
    gen <- hsw_last_gen
    merged_stem <- 'hsw_merged'
    merged_dir <- file.path(peds_dir, 'merged', 'hsw', paste0('gen', hsw_last_gen))
    current_raw_ped <- file.path(hsw_raw_dir, paste0(hsw_raw_stem, hsw_last_gen, '.csv'))
    out_dir <- hsw_dir
} else if (pop == 'wfu') {
    gen <- wfu_last_gen
    merged_stem <- 'wfu_merged'
    merged_dir <- file.path(peds_dir, 'merged', 'wfu', paste0('gen', wfu_last_gen))
    current_raw_ped <- file.path(wfu_raw_dir, paste0(wfu_raw_stem, wfu_last_gen, '.csv'))
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
    hsw_first_gen,
    hsw_last_gen,
    data_dir = merged_dir,
    file_stem = merged_stem,
    write_file = F,
    return_ids = F,
    print_ids = F)

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

# estimate kinship, copy files to the final output directory
printout('Estimating kinship')

k <- current.kinship(
    df = hsw_current_raw_ped, 
    first_gen = hsw_first_gen,
    last_gen = hsw_last_gen,
    data_dir = merged_dir,
    file_stem = merged_stem,
    id_map = merged_ped$id.map)

kfiles <- list.files(merged_dir, pattern='kinship', full.names=T)
sapply(kfiles, function(file) file.copy(from = file, to = kinship_outdir))
kfile <- read.csv(kfiles[which(grepl('all_pairings',kfiles))])

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
    first_gen = hsw_first_gen,
    last_gen = hsw_last_gen,
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

# create the final breederpair file
printout('Producing final breederpair file')

# WFU sample sheet parameter to automate pairings file creation 
wfu_ss_input <- NULL
if (exists('wfu_ss') && !is.null(wfu_ss) && !is.na(wfu_ss)) {
  wfu_ss_input <- wfu_ss
}

breedpairs <- create_hsw_breeder_file(
    pairs = breedpairs_animalid, 
    df = colony_df, 
    wfu_ss = wfu_ss_input,    
    outdir = breedhs_outdir)$pairs

# identify unpaired breeders
bp_counts <- count_breeder_pairs(
    assignments = assignment_file,
    breederpairs = breedpairs)

unpaired_ids <- bp_counts$assigned.but.not.paired
n_paired <- nrow(breedpairs)
all_paired_ids <- c(breedpairs$dam_animalid, breedpairs$sire_animalid)

# get best aternative pairings

printout('Identifying best alternative pairings for all breeders')
best_alt_pairs(
    pairs = breedpairs, 
    kinship = kfile, 
    avail_ids = avail_ids, 
    ids_to_replace = NULL,
    ids_to_pair = all_paired_ids,
    n_best = 4,
    outdir = altpair_outdir)

printout('Identifying best replacement pairings among available rats')
best_alt_pairs(
    pairs = breedpairs, 
    kinship = kfile, 
    avail_ids = avail_ids, 
    ids_to_replace = all_paired_ids,
    ids_to_pair = NULL,
    n_best = 4,
    outdir = altpair_outdir)

printout('Identifying best new pairings among available rats')
best_alt_pairs(
    pairs = breedpairs, 
    kinship = kfile, 
    avail_ids = avail_ids, 
    ids_to_replace = NULL,
    ids_to_pair = NULL,
    n_best = 6,
    outdir = altpair_outdir)


printout('Breeder pairing complete! Pairing file is ready for use in the colony')

