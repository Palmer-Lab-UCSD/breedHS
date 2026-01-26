library(argparse)

# get arguments
all_args <- commandArgs(trailingOnly = T)

# find the arg_parser argument
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

cat(paste0('Producing final breeders file for HS West gen', generation), '\n')
datestamp <- format(Sys.time(),'%Y%m%d-%H:%M:%S')
cat(datestamp, '\n\n')

# set up variables, paths, directories
breedhs_outdir <- file.path(proj_dir, 'breedhs_out')
altpair_outdir <- file.path(proj_dir, 'alt_pairings')
kinship_outdir <- file.path(proj_dir, 'kinship')
colony_dir <- file.path(proj_dir, 'colony_edits')
hsw_raw_dir <- file.path(peds_dir, 'raw', 'hsw')
pairs_proposed <- read.csv(pairs_proposed)
kfile <- list.files(kinship_outdir, pattern='kinship_all_pairings', full.names=T)


# read in any alternative pairings files filled out in the colony
colony_files <- list.files(colony_dir, full.names=T)
for (file in colony_files) {
    filestem <- basename(file)
    if (grepl('breeders_proposed', filestem)) {
        cat('Using pairings file:', file, '\n\n')
        colony_pairs <- read.csv(file) } 
    if (grepl('alternative_pairs', filestem)) {
        cat('Using alternative pairs file:', file, '\n')
        colony_alt_pairs <- read.csv(file) }
    if (grepl('replacement_pairs', filestem)) {
        cat('Using replacement pairs file:', file, '\n')
       colony_rep_pairs <- read.csv(file) }
    if (grepl('new_pairs', filestem)) {
        cat('Using novel pairs file:', file, '\n')
        colony_new_pairs <- read.csv(file) }
}

if (!exists('colony_pairs')) {
    cat('ERROR: No pairings file. Please provide the pairings file used in the colony \n')
    quit(status=1)
}
if (!exists('colony_alt_pairs')) { 
    colony_alt_pairs <- NULL
    cat('No alternative pairs file used \n')
    }
if (!exists('colony_rep_pairs')) { 
    colony_rep_pairs <- NULL
    cat('No replacement pairs file used \n')
    }
if (!exists('colony_new_pairs')) { 
    cat('No novel pairs file used \n')
    colony_new_pairs <- NULL
    }


# update the breeder file with colony edits
printout('Creating final breeder file')
breedpairs <- final_breeder_file(
    proposed_pairs = pairs_proposed, 
    colony_pairs = colony_pairs,
    colony_df = colony_df,
    alt_pairs = colony_alt_pairs,
    rep_pairs = colony_rep_pairs,
    new_pairs = colony_new_pairs,
    pop = pop,
    gen = generation,
    outdir = proj_dir)


# plot kinship 
printout('Plotting pairwise kinship for all breeder pairings')

plot_k_hist(
    kinship = kfile, 
    pop = pop,
    gen = generation,
    sample = 'all',
    out_dir = kinship_outdir)

plot_k_hist(
    kinship = breedpairs, 
    pop = pop,
    gen = generation,
    sample = 'breederpairs',
    out_dir = kinship_outdir)

plot_k_network(
    kinship = breedpairs, 
    pop = pop,
    gen = generation,
    sample = 'breederpairs',
    out_dir = kinship_outdir)


# finalize the current generation's pedigree
printout(paste0('Producing the finalized raw pedigree for gen', generation))

hsw_raw_ped <- final_hsw_breeders_to_raw_ped(
    pairs = breedpairs,
    colony_df = colony_df,
    prev_df = prev_colony_df,
    stem = hsw_raw_stem,
    outdir = hsw_raw_dir)

printout(paste0('Producing the finalized raw pedigree for gen', generation))
concat_peds(
    directory = hsw_raw_dir,
    stem = hsw_raw_stem,
    outdir = hsw_raw_dir,
    outstem = 'hsw_raw_ped_complete',
    return_df = FALSE)

printout('Breederpair updates complete!')
