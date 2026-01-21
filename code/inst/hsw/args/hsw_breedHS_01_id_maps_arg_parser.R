library(argparse)

# create a parser closure for input arguments
breedhs_parser <- argument_parser(
    argument_def('breedhs', type = 'character', default_val = 'breedhs',
        help = 'arbitrary positional argument included to enable argument parsing'),
    argument_def('--wfu_raw_ped', type = 'character', required = TRUE,
        help = 'the file path to the complete raw WFU pedigree'),
    argument_def('--all_wfu_ss', type = 'character', required = TRUE, nargs=3,
        help = 'bash array of file paths to WFU shipping sheets'),
    argument_def('--hsw_raw_ped', type = 'character', required = TRUE,
        help = 'the file path to the complete raw HSW pedigree'),
    argument_def('--all_hsw_ss', type = 'character', required = FALSE, nargs=1,
        help = 'bash array of file paths to WFU shipping sheets'),
    argument_def('--hsw_assignments', type = 'character', required = FALSE,
        help = 'path to HSW assignments file'),
    argument_def('--prev_colony_df', type = 'character', required = FALSE,
        help = 'path to the previous generations HSW colony dataframe'),
    argument_def('--utils', type = 'character', required = TRUE,
        help = 'path to the main breedHS.R script'),
    argument_def('--outdir', type = 'character', required = TRUE,
        help = 'output directory path'),
    argument_def('--arg_parser', type = 'character', required = TRUE,
        help = 'path to the argument parser R script')
)
