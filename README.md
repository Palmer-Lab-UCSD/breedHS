# Breeder selection in HS rats :rat: :cupid: :rat:

This repository contains R source code to select breeders in an
outbred line so that inbreeding in the offspring is
minimized. This repository is forked from [breedail](https://github.com/pcarbo/breedail),
which was originally developed for advanced intercross lines in mice, and is currently
used by the [HS West](https://ratgenes.org/) colony, managed by 
[Dr. Abraham Palmer's lab](https://palmerlab.org/) at UC San Diego, for breeding 
Heterogeneous Stock (HS) rats.

The core breedail algorithm is unchanged here, but this repository provides three main 
additions to the original breedail repo:

1. Utility functions to make breeder selection more user-friendly, without the need
for substantial code and manual re-runs (such as in the [original example](code/orig_breedail_example.R) 
code provided with breedail).

2. A novel algorithm for merging pedigrees following breeder exchanges between separate populations.

3. Functions to simulate forward breeding and breeder selection in theoretical future generations.

## Background: HS rats :rat:

The HS West colony (HSW) was founded using HS rats provided to UCSD from Wake Forest University (WFU),from 
a subset of WFU generation 42. That is, some WFU generation 42 rats were used to produce 
WFU generation 43, and other WFU generation 42 rats were sent to HSW (event 'A' in the figure below). 
The HSW colony follows a naming convention that counts HSW generations based on the number of generations since 
the *original* founding of the HS rats strain, not since the founding of the WFU colony. So, WFU gen042 rats sent to
HSW comprise HSW gen096. Both colonies have been maintained separately (with occasional exchanges of rats betweeen 
the two, see below), following separate generation times. As such, there is no direct equivalence between any WFU 
generations after WFU gen042 and any HSW generations after HSW gen096. While an HSW colony didn't exist prior to 
gen096, for the sake of a consistent numbering convention required for breedail, WFU generations before gen042 can 
also be considered HSW generations 54-95. 

<img src="readme_resources/hsw_colony_history_breedHS.png" width=350 height=400>

The two colonies regularly exchange animals to breed between them in order to minimize drift between colonies and 
maintain outbreeding across the greater HS rat population. In March 2024, WFU sent 43 male rats to HSW (event 'B' above) 
as the first exchange. Different generation times and ID conventions between the colonies require a complicated 
merging algorithm in order to produce a joint pedigree usable for breedail. Here, breedHS automates these merges to 
enable kinship estimation and breeder selection using breedail. 

## Files

* [kinship.R](code/kinship.R) defines functions to calculate and update pairwise kinship coefficients. These are 
core breedail functions unmodified from the original repository.

* [find_mates.R](code/find_mates.R) defines functions to pairs breeders in a given generation of the pedigree in such 
a way that they will produce the smallest average inbreeding coefficients in the offspring. These are 
core breedail functions unmodified from the original repository.

* [utils.R](code/utils.R) defines pre- and post-processing functions to format data input and output from breedail, and
wrapper functions to facilitate the execution of breedail's core algorithm. These are all generalizable to any 
population using breedail for kinship estimation and breeder selection, and/or any pair of populations for whom a 
common, joint pedigree is needed. 
Some important functions:
    -   `find.ped.errors`:
    -   `format.pedigree`:
    -   `select.breeders`:
    -   `mate.breeders`:
    -   `format.pedigree.for.merge`:
    -   `merge.pedigrees`:

* [HSW_utils.R](code/HSW_utils.R) defines pre-processing functions specifically tailored to data formatting conventions at HSW and WFU.
While these are not generalizable to other populations as-is, they may provide a useful template for other 
investigators interested in using breedail. 

## Example data
:construction: :construction: 
A comprehensive example workflow is currently in progress. Come back soon!
:construction: :construction: 

## Usage
:construction: :construction: 
A comprehensive example workflow is currently in progress. Come back soon!
:construction: :construction: 

See [orig_breedail_example.R](code/orig_breedail_example.R) for an example of how to use these functions to select 
mouse breeders in generation F10 of an advanced intercross line. The pedigree data for this mouse population are found
in [data](data). This script should generate matings similar to those in [F10pairs.txt](results/F10pairs.txt).
