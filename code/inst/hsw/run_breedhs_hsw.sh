#!/bin/bash
proj_dir=/tscc/projects/ps-palmer/hs_rats/rattaca/breeders/gen105_106/test
breedhs_args=${proj_dir}/code/breedHS_args_gen105_106.cfg
source ${breedhs_args}
workdir=/tscc/projects/ps-palmer/hs_rats/rattaca/breeders/gen105_106/test/code
source activate ${conda_env}

# if no output directory exists, set up breedHS directories and pedigrees
if [ ! -d "${breedhs_out}" ]; then

    Rscript ${workdir}/breedHS_01_setup_peds.R \
        --pop ${pop} \
        --utils ${utils} \
        --arg_parser ${arg_parser} \
        --peds_dir ${peds_dir} \
        --wfu_raw_ped ${wfu_raw_ped} \
        --assignment_file ${assign_file} \
        --hsw_first_gen ${hsw_first_gen} \
        --hsw_last_gen ${hsw_last_gen} \
        --hsw_stem ${hsw_stem} \
        --hsw_raw_stem ${hsw_raw_stem} \
        --prev_pairs ${prev_pairs} \
        --prev_colony_df ${prev_colony_df} \
        --proj_dir ${proj_dir} \
        breedhs

fi

# if no proposed alternative pairings have been made,
# initiate the first round of breeder pairing
if [ -z "${proposed_pairs+x}" ]; then

    Rscript ${workdir}/breedHS_02_propose_pairs.R \
        --pop ${pop} \
        --n_pairs ${n_pairs} \
        --arg_parser ${arg_parser} \
        --inputs ${inputs} \
        --peds_dir ${peds_dir} \
        --ped_map ${ped_map} \
        --hsw_to_wfu ${hsw_to_wfu} \
        --wfu_to_hsw ${wfu_to_hsw} \
        --hsw_stem ${hsw_stem} \
        --hsw_raw_stem ${hsw_raw_stem} \
        --hsw_first_gen ${hsw_first_gen} \
        --hsw_last_gen ${hsw_last_gen} \
        --wfu_raw_ped ${wfu_raw_ped} \
        --wfu_id_map ${wfu_id_map} \
        --wfu_stem ${wfu_stem} \
        --wfu_raw_stem ${wfu_raw_stem} \
        --wfu_first_gen ${wfu_first_gen} \
        --wfu_last_gen ${wfu_last_gen} \
        --utils ${utils} \
        --colony_df ${colony_df} \
        --assignment_file ${assign_file} \
        --avail_ids ${avail_ids} \
        --proj_dir ${proj_dir} \
        breedhs

# if proposed_pairs is set, then log all colony pairings into the final breeders file
elif [ ! -z ${proposed_pairs+x} ]; then
    
    echo "hsw_last_gen: ${hsw_last_gen}"
    Rscript ${workdir}/breedHS_03_finalize_pairs.R \
        --pop ${pop} \
        --generation ${hsw_last_gen} \
        --colony_dir ${colony_dir} \
        --pairs_proposed ${proposed_pairs} \
        --colony_pairs ${colony_pairs} \
        --colony_df ${colony_df} \
        --utils ${utils} \
        --proj_dir ${proj_dir} \
        --arg_parser ${arg_parser} \
        breedhs

fi

conda deactivate