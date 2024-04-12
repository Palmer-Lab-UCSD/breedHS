# Breeder selection in HS rats :rat: :cupid: :rat:
<br>
This repository contains R source code to select breeders in an
outbred line so that inbreeding in the offspring is
minimized. This repository is forked from [breedail](https://github.com/pcarbo/breedail),
which was originally developed for advanced intercross lines in mice, and is currently
used by the [HS West](https://ratgenes.org/) colony, managed by 
[Dr. Abraham Palmer's lab](https://palmerlab.org/) at UC San Diego, for breeding HS rats.
<br>
The core breedail algorithm is unchanged here, but this repository provides two main 
additions to the original breedail repo:
1. Utility functions to make breeder selection more user-friendly, without the need
for substantial code and manual re-runs (such as in the [original example](code/orig_breedail_example.R) 
code provided with breedail).
2. A novel algorithm for selecting breeders to exchange between two separate populations.
<br>

* [kinship.R](code/kinship.R) defines functions to calculate and
update pairwise kinship coefficients.

* [find_mates.R](code/find_mates.R) pairs breeders in a given
generation of the advanced intercross line in such a way that they
will produce the smallest average inbreeding coefficients in the
offspring.

See [orig_breedail_example.R](code/orig_breedail_example.R) for an example of how to use these
functions to select mouse breeders in generation F10 of an advanced
intercross line. The pedigree data for this mouse population are found
in [data](data). This script should generate matings similar to those
in [F10pairs.txt](results/F10pairs.txt).


