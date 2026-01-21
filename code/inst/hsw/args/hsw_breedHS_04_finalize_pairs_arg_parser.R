library(argparse)

# create a parser closure for input arguments
breedhs_parser <- argument_parser(
    argument_def('breedhs', type = 'character', default_val = 'breedhs',
        help = 'arbitrary positional argument included to enable argument parsing'),
    argument_def('--pop', type = 'character', required = TRUE,
        help = 'the population being analyzed: the population to merge pedigree info into'),
    argument_def('--peds_dir', type = 'character', required = TRUE, 
        help = 'path to the upper directory for pedigrees'),
    argument_def('--hsw_raw_stem', type = 'character', required = TRUE,
        help = 'the file stem name for raw HSW pedigrees'),
    argument_def('--generation', type = 'integer', required = TRUE,
        help = 'generic generation variable to accommodate a specific desired generation'),
    argument_def('--utils', type = 'character', required = TRUE,
        help = 'path to the main breedHS.R script'),
    argument_def('--replace_ids', type = 'character', required = FALSE,
        help = 'the path to the file listing animal IDs that need replacement for pairing'),
    argument_def('--pairs_proposed', type = 'character', required = TRUE,
        help = 'path to the original proposed pairings file produced before colony pairing'),
    argument_def('--proj_dir', type = 'character', required = TRUE,
        help = 'the path to the top directory for the current breeders generation'),
    argument_def('--colony_df', type = 'character', required = TRUE,
        help = 'path to the most recent HSW colony dataframe of ALL rats (not only breeders)'),
    argument_def('--prev_colony_df', type = 'character', required = TRUE,
        help = 'path to the previous generations HSW colony dataframe of ALL rats (not only breeders)'),
    argument_def('--arg_parser', type = 'character', required = TRUE,
        help = 'path to the argument parser R script')

)
