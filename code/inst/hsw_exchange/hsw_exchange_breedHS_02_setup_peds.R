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
counts_outdir <- file.path(breedhs_outdir, 'bp_counts')
altpair_outdir <- file.path(proj_dir, 'alt_pairings')
kinship_outdir <- file.path(proj_dir, 'kinship')
colony_dir <- file.path(proj_dir, 'colony_edits')
dir.create(breedhs_outdir, showWarnings=F)
dir.create(counts_outdir, showWarnings=F)
dir.create(altpair_outdir, showWarnings=F)
dir.create(kinship_outdir, showWarnings=F)
dir.create(colony_dir, showWarnings=F)

# create pedigree directories
breedail_dir <- file.path(peds_dir, 'breedail')
hsw_dir <- file.path(breedail_dir, 'hsw')
hsw_raw_dir <- file.path(peds_dir, 'raw', 'hsw')
wfu_dir <- file.path(breedail_dir, 'wfu')
wfu_raw_dir <- file.path(peds_dir, 'raw', 'wfu')
dir.create(breedail_dir, recursive=F, showWarnings=F)
dir.create(hsw_dir, recursive=F, showWarnings=F)
dir.create(wfu_dir, recursive=F, showWarnings=F)

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
dir.create(merged_dir, recursive=T, showWarnings=F)
dir.create(current_dir, showWarnings=F)

# # incorporate HSW rats into the WFU raw pedigree - not executed gen107
# printout('Adding HSW rats to the raw WFU pedigree')
# add_hsw_rats_to_wfu_raw_ped(
#     ped = wfu_raw_ped,
#     hsw_shipping_sheet = hsw_shipping_sheet,
#     add_to_gen = wfu_last_gen,
#     outdir = wfu_raw_dir)

# # new raw pedigree file
# datestamp <- format(Sys.time(),'%Y%m%d')
# wfu_raw_ped <- file.path(wfu_raw_dir, paste0('wfu_raw_ped_complete_', datestamp, '.csv'))


# convert the WFU raw pedigree into breedail format as part of the HSW pedigree
printout('Converting the raw WFU pedigree')
format_wfu_raw_ped(
    ped = wfu_raw_ped,
    wfu_map = wfu_id_map,
    outdir = wfu_dir
)

# save the complete formatted WFU pedigree
concat_peds(
    directory = wfu_dir,
    stem = wfu_stem,
    outdir = file.path(wfu_dir, paste0('gen', wfu_last_gen)),
    outstem = 'wfu_ped_complete',
    return_df=FALSE
)

# convert the most recent HSW assignments into a raw pedigree
printout('Converting HSW assignments to pedigree')
hsw_current_raw_ped <- assignment_to_raw_ped(
    assignments = hsw_assignments,
    hsw_map = hsw_id_map, 
    wfu_ss = wfu_shipping_sheet,
    wfu_map = wfu_id_map,
    outdir = hsw_raw_dir)

# # incorporate WFU rats into the HSW raw pedigree
# printout('Adding WFU rats to the raw HSW pedigree')
# wfu_into_hsw_gen(
#     hsw_raw = hsw_assignments,
#     wfu_ss = wfu_shipping_sheet,
#     wfu_sheet = 'for ',
#     ss_gen = wfu_last_gen,
#     wfu_prev = wfu_raw_ped,
#     outdir = hsw_raw_dir,
#     return_df = FALSE)

# convert the HSW raw pedigree into breedail format
printout('Converting the raw HSW pedigree')
for (gen in 95:hsw_last_gen) {

    gen <- ifelse(nchar(gen)==2, paste0('0',gen), gen)
    hsw_raw_ped_use <- paste0('hsw_raw_gen', gen, '.csv')
    hsw_raw_ped_use <- file.path(hsw_raw_dir, hsw_raw_ped_use)

    format_hsw_raw_ped(
        ped = hsw_raw_ped_use,
        wfu_map = wfu_id_map, 
        hsw_map = hsw_id_map, 
        return_df = F,  
        outdir = hsw_dir)

}

# incorporate WFU generations into the HSW pedigree (prior to founding HSW)
printout('Incorporating WFU generations into the HSW pedigree')
wfu_raw_to_hsw(
    ped = wfu_raw_ped, 
    wfu_map = wfu_id_map,
    outdir = hsw_dir)

# save the complete formatted HSW pedigree
concat_peds(
    directory = hsw_dir,
    stem = hsw_stem,
    outdir = file.path(hsw_dir, paste0('gen', hsw_last_gen)),
    outstem = 'hsw_ped_complete',
    return_df=FALSE
)

# check the WFU pedigree for errors
printout('Checking the WFU pedigree for errors')
find.ped.errors(
    first_gen = wfu_first_gen,  # the WFU founder generation
    last_gen = wfu_last_gen,   # the current WFU generation
    data_dir = wfu_dir,   
    file_stem = wfu_stem, 
    write_file = TRUE,
    return_ids = FALSE,
    print_ids = TRUE)

# check the HSW pedigree for errors
printout('Checking the HSW pedigree for errors')
find.ped.errors(
    first_gen = hsw_first_gen,  # the HSW founder generation
    last_gen = hsw_last_gen,   # the current HSW generation
    data_dir = hsw_dir,   
    file_stem = hsw_stem, 
    write_file = TRUE,
    return_ids = FALSE,
    print_ids = TRUE)

# save copies of complete pedigrees alongside the other colony's pedigree
hsw_ped <- list.files(file.path(hsw_dir, paste0('gen', hsw_last_gen)),full.names=T)
wfu_ped <- list.files(file.path(wfu_dir, paste0('gen', wfu_last_gen)),full.names=T)
file.copy(from = hsw_ped, to = file.path(wfu_dir, paste0('gen', wfu_last_gen)))
file.copy(from = wfu_ped, to = file.path(hsw_dir, paste0('gen', hsw_last_gen)))

printout('Pedigree setup complete. Ready for pair selection')