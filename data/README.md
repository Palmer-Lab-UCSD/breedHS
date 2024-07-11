`ped_map.csv`
1. `pop1_gen`: All pop1 generations with a corresponding pop2 generation, i.e. the
    pedigree history shared between the two populations, denoted in pop1 generations
2. `pop2_gen`: all pop2 generations with a corresponding pop1 generation, i.e. the
    pedigree history shared between the two populations, denoted in HSW generations  

`pop1_pop2_exchange.csv`
1. `pop1_from`: the complete pedigree history at pop1 (all pop1 generations). Any rows 
    in this column with a filled `pop2_to` value reflect pop1 generations from which rats
    were sent to pop2.
2. `pop2_to`: all pop2 generations that received breeders from pop1, on the row corresponding 
    to the generation (`pop1_from`) from which the rats were received  
    
`pop2_pop1_exchange.csv`
1. `pop2_from`: the complete pedigree history at pop2 (all pop2 generations). Any rows 
    in this column with a filled `pop1_to` value reflect pop2 generations from which rats
    were sent to pop1.
2. `pop1_to`: all pop1 generations that received rats from pop2, on the row corresponding 
    to the generation (`pop1_from`) from which the rats were received