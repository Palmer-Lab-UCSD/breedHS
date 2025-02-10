## PEDIGREE-RELATED FUNCTIONS SPECIFIC TO WAKE FOREST
## written by Dr. Ben Johnson (bbjohnson@health.ucsd.edu)

## These are used in conjunction with utils.R, which is more general
## These functions are designed to adhere specifically to data and organizational
## conventions used in-house at the WFU colony. They are likely not useful
## to outside users, though utils.R is of general utility once pedigrees have already
## been properly formatted for use in breedail

source('utils.R')
library(readxl)



# convert WFU SW.ID to access ID
swid_to_accessid <- function(id, wfu_map) {

    accessid <- wfu_map[wfu_map[['swid']] == id,][['accessid']]
    return(accessid)
}

# convert access ID to WFU SW.ID
accessid_to_swid <- function(id, wfu_map) {

    swid <- wfu_map[wfu_map[['accessid']] == id,][['swid']]
    return(swid)
}
