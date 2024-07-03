## BREEDER SELECTION CODE
## breedail written by Dr. Andrew Skol (askol at uchicago.edu) 
## edited for use at HS West by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

#######################################################################

## This script selects breeders from Heterogeneous Stock rats using kinship
## data based on a pedigree. This is the script you will use to obtain a list 
## of breeder pairs for each generation. 

## It calls on three other scripts: (1) kinship.R, which calculates pairwise
## kinship coefficients for members of the current generation, (2) find_mates.R, 
## which pairs breeders with the smallest average inbreeding coefficients 
## (and can also try to pair animals that will produce the minimum average 
## inbreeding coefficient in the offspring), and (3) utils.R, which provides
## utility functions that wrap around the original breedail functions found in 
## the first two files. 
## Note: files (1) and (2) are unedited code as distributed with breedail. File
## (3) was written by Ben Johnson to simplify execution of the program 


#########################################################################
## ARGUMENTS ##
#########################################################################

## First set the working directory (see function "setwd") to the
## location where the code and scripts are.

setwd('/path/to/breedail/code/')
source('utils.R')

## Next set the directory that holds all pedigree files and the stem name 
## for all pedigree files. Each file must have the same stem name, followed 
## by consecutive generation numbers, e.g. 'ped2', 'ped3', etc. Each file must
## have IDs for the current generation stored in the first column, and must 
## contain columns labeled 'dam', 'sire', and 'sex' anywhere in the file.
## Missing values must be encoded as '?'
## This example uses the example dataset originally provided with breedail, 
## which entails generations 2-10 of an example pedigree.

data.dir <- file.path('..', 'data')
file.stem <- 'ped'

## Set your desired output directory, where results will be saved

out.dir <- file.path('..', 'breedpairs')

## Set your desired first and last generations to analyze. Kinship will be 
## calculated across all generations. Breeders will be selected from the 
## final generation

first.gen <- 2
last.gen <- 10

## Set the desired number of selection rounds. Breedail will iterate this 
## number of times to identify breeder pairs. Usually a handful of iterations 
## are needed

num.rounds <- 4


#########################################################################
## ANALYSIS ##
#########################################################################

## Execute this function to conduct all breedail breeder selection.
## The function takes as input all of the variables declared above, plus the
## boolean option 'one_per_sibship'. If one_per_sibship=TRUE, the first round
## of breeder selection will be restricted to only one individual per sibship.
## If one_per_sibship=FALSE, the first round will allow multiple siblings to
## used as breeders. If this option is not explicitly set, the function will
## by default set to TRUE

breedpairs <- select.breeders(first.gen, 
                      last.gen, 
                      data.dir, 
                      out.dir,
                      file.stem, 
                      num.rounds, 
                      one_per_sibship=T)
