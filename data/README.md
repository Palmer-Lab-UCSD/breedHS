Readme: column keys for all files in `breedHS/code/data`

`ped_map.csv`
1. `wfu_gen`: All WFU generations with a corresponding HSW generation, i.e. the
    pedigree history shared between the two populations, denoted in WFU generations
2. `hsw_gen`: all HSW generations with a corresponding WFU generation, i.e. the
    pedigree history shared between the two populations, denoted in HSW generations  

`wfu_hsw_exchange.csv`
1. `wfu_from`: the complete pedigree history at WFU (all WFU generations). Any rows 
    in this column with a filled `hsw_to` value reflect WFU generations from which rats
    were sent to HSW.
2. `hsw_to`: all HSW generations that received breeders from WFU, on the row corresponding 
    to the generation (`wfu_from`) from which the rats were received  
    
`hsw_wfu_exchange.csv`
1. `hsw_from`: the complete pedigree history at HSW (all HSW generations). Any rows 
    in this column with a filled `wfu_to` value reflect HSW generations from which rats
    were sent to WFU.
2. `wfu_to`: all WFU generations that received rats from HSW, on the row corresponding 
    to the generation (`hsw_from`) from which the rats were received