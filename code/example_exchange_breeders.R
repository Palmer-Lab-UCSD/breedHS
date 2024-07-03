# EXAMPlE CODE TO EXCHANGE BREEDERS BETWEEN TWO PEDIGREES

# In this example Population 1 is the original, founder population.
# This population was bred for 31 generations. Some animals from generation 32
# were used to found Population 2. These two populations were maintained separately
# with different generation times, such that at the time of breeder exchange,
# Population 1 is at generation 34 and Population 2 is at generation 35.

# to select breeders for an exchange requires setting multiple options:

# dir_1 
#   The path directory housing all pedigree files for population 1.
#   Pop 1 should be the 'original' population used to found the second
#   population used in for exchange. SFiles must be separated one per generation, 
#   and all files must have the same stem name followed by a generation number, 
#   e.g. 'ped0', 'ped1', 'ped2', etc

# stem_1
#   The file name stem used for pedigree files for population 1. All files must
#   have the same stem name followed by a generation number

# dir_2
#   The directory housing all pedigree files for pop 2. Same formatting requirements
#   as for pop 1.

# stem_2
#   The file name stem used for pedigree files in pop 2. Same formatting requirements
#   as for pop 1

# first_gen_1
#   The first generation to include from pop 1. All generations starting from this
#   first one to the last will be used to merge pedigrees and identify breeders.
#   Note, it is okay if both populations share data from previous generations
#   (for example, Pop 2 started at generation 32 and keeps Pop 1 pedigree info for
#   generations 0-31). These overlapping generations can be included, and duplicate
#   entries will be automatically removed from analysis

# last_gen_1
#   The final generation to include from Pop 1. This will be the most recent 
#   generation, from which breeders will be selected to send to Pop 2.

# first_gen_2
#   The first generation to include from Pop 2.

# last_gen_2
#   The final generation to include from Pop 2, from which breeders will be 
#   selected to send to Pop 1.

# out_dir
#   The output directory in which to save files listing breeder pairs

# out_stem
#   The stem file name used for output files

# one_per_sibship
#   A boolean option denoting whether to allow only one or more than one individual 
#   per sibship to be used for breeding. The default is TRUE, meaning only one 
#   offspring per sibship will be selected. If FALSE, multiple siblings can be 
#   chosen as breeders

## note: The function exchange.breeders() is coded to run four rounds of breeder
#   selection. Outputs will display the round from which each breeder pair was
#   identified



# first, set your working directory to read in required breedail and utility functions
workdir <- '/path/to/breedail/code/'
setwd(workdir)
source('utils.R')

# next, set your data directories and options:
ped_dir <- '/path/to/breedail/pedigrees'
pop1_dir <- file.path(ped_dir, 'pop1')
pop2_dir <- file.path(ped_dir, 'pop2')
outdir <- '/path/to/breedail/breedpairs'
pop1_stem <- 'pop1_gen'
pop2_stem <- 'pop2_gen'
pop1_first_gen <- 0
pop1_last_gen <- 34
pop2_first_gen <- 0
pop2_last_gen <- 35
out_stem <- 'pop1_pop2'
one_per_sibship <- FALSE


# identify breeder pairs using the previously-set options
exchange.breeders(
    dir_1 = pop1_dir, 
    stem_1 = pop1_stem, 
    dir_2 = pop2_dir, 
    stem_2 = pop2_stem,
    first_gen_1 = pop1_first_gen, 
    last_gen_1 = pop1_last_gen, 
    first_gen_2 = pop2_first_gen, 
    last_gen_2 = pop2_last_gen, 
    out_dir = outdir, 
    out_stem = out_stem, 
    one_per_sibship = one_per_sibship)


# descriptions of output files:

# 'M1F2' file
#   Files including the flag 'M1F2' in their name show breeding pairs with males
#   from Pop 1 and females from Pop 2. Each pair's kinship and the breedail round 
#   during which the pair was identified are also included. The remainder of the 
#   file name will include the output stem set by the user and the number of 
#   offspring per sibship.

# 'M2F1' file
#   Files with 'M2F1' in their name show breeding pairs with males from Pop 2
#   and females from Pop 1, formatted identically as the M1F2 file
