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

# create output directories
breedhs_outdir <- file.path(proj_dir, 'breedhs_out')
altpair_outdir <- file.path(proj_dir, 'alt_pairings')
kinship_outdir <- file.path(proj_dir, 'kinship')
colony_dir <- file.path(proj_dir, 'colony')
dir.create(breedhs_outdir, showWarnings=F)
dir.create(altpair_outdir, showWarnings=F)
dir.create(kinship_outdir, showWarnings=F)
dir.create(colony_dir, showWarnings=F)

# intermediate directories
hsw_dir <- file.path(peds_dir, 'breedail', 'hsw')
hsw_raw_dir <- file.path(peds_dir, 'raw', 'hsw')
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


# convert the WFU raw pedigree into breedail format as part of the HSW pedigree
# (this formats WFU generations 00-47)
printout('Converting the raw WFU pedigree')
format_wfu_raw_ped(
    ped = wfu_raw_ped,
    outdir = wfu_dir
)

# convert the most recent HSW assignments into a raw pedigree
printout('Converting HSW assignments to pedigree')
hsw_current_raw_ped <- assignment_to_raw_ped(
    assignments = assignment_file, 
    outdir = hsw_raw_dir)

# convert the most recent HSW generation to breedail format
printout('Converting the raw HSW pedigree')
hsw_current_ped <- format_hsw_raw_ped(
    df = hsw_current_raw_ped,
    wfu_map = NULL, 
    return_df = T, 
    outdir = hsw_dir)

# convert the previous pairings file into a raw pedigree
# NOTE: this should usually be done at the end of pairing, but was not done in gen104 so doing it now
printout('Converting gen104 pairings to raw HSW pedigree')
hsw_prev_raw_ped <- final_hsw_breeders_to_raw_ped(
    pairs = prev_pairs,
    colony_df = prev_colony_df,
    stem = hsw_raw_stem,
    outdir = hsw_raw_dir)

# convert the previous generation's raw pedigree to breedail format
printout('Converting gen104 raw pedigree')
hsw_prev_ped <- format_hsw_raw_ped(
    df = hsw_prev_raw_ped,
    wfu_map = NULL,
    return_df = F,
    outdir = hsw_dir)

# check the HSW pedigree for errors
find.ped.errors(
    first_gen = 95,  # the HSW founder generation
    last_gen = hsw_last_gen,   # the current HSW generation
    data_dir = hsw_dir,   
    file_stem = hsw_stem, 
    write_file = FALSE,
    return_ids = FALSE,
    print_ids = TRUE)

printout('Pedigree setup complete. Ready for pair selection')