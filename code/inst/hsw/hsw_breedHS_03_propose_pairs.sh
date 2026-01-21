#!/bin/bash

START=$(date +%s)
source ${breedhs_args}
source activate ${conda_env}

echo ""
echo "--------------------------------------------------------"
echo "---------------- HS WEST BREEDER PAIRING ---------------"
echo "--------------------------------------------------------"
echo ""
echo "Proposing breeder pairs for HSW generation: ${hsw_last_gen}"
echo ""
echo "$(date +"%Y%m%d-%H:%M:%S")"
echo "user: ${user}"
echo ""


if [[ "$skip_k" == "false" ]]; then

    Rscript ${code}/hsw_breedHS_03_propose_pairs.R \
        --pop ${pop} \
        --n_pairs ${n_pairs} \
        --arg_parser ${arg_parser_03} \
        --inputs ${inputs} \
        --peds_dir ${peds_dir} \
        --ped_map ${ped_map} \
        --hsw_to_wfu ${hsw_to_wfu} \
        --wfu_to_hsw ${wfu_to_hsw} \
        --hsw_id_map ${hsw_id_map} \
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
        --colony_df ${colony_df} \
        --hsw_assignments ${hsw_assignments} \
        --wfu_shipping_sheet ${wfu_shipping_sheet} \
        --avail_ids ${avail_ids} \
        --utils ${utils} \
        --proj_dir ${proj_dir} \
        breedhs

elif [[ "$skip_k" == "true" ]]; then

    kinship_file=${proj_dir}/kinship/hsw_merged_gen${hsw_last_gen}_kinship_all_pairings.csv
    echo "Using kinship file: ${kinship_file}"
    echo ""

    Rscript ${code}/hsw_breedHS_03_propose_pairs.R \
        --pop ${pop} \
        --n_pairs ${n_pairs} \
        --arg_parser ${arg_parser_03} \
        --inputs ${inputs} \
        --peds_dir ${peds_dir} \
        --ped_map ${ped_map} \
        --hsw_to_wfu ${hsw_to_wfu} \
        --wfu_to_hsw ${wfu_to_hsw} \
        --hsw_id_map ${hsw_id_map} \
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
        --colony_df ${colony_df} \
        --hsw_assignments ${hsw_assignments} \
        --wfu_shipping_sheet ${wfu_shipping_sheet} \
        --avail_ids ${avail_ids} \
        --utils ${utils} \
        --proj_dir ${proj_dir} \
        --kinship_file ${kinship_file} \
        breedhs

fi

conda deactivate


END=$(date +%s)
echo "Time elapsed for breeder pairing: $(( $END - $START )) seconds"
