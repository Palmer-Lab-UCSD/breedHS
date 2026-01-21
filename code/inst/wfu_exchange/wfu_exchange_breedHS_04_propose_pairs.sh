#!/bin/bash

START=$(date +%s)
source ${breedhs_args}
source activate ${conda_env}

echo ""
echo "----------------------------------------------------------"
echo "--------------- WAKE FOREST BREEDER PAIRING --------------"
echo "----------- (WHEN RECEIVING RATS FROM HS WEST) -----------"
echo "----------------------------------------------------------"
echo ""
echo "Proposing breeder pairs for WFU generation: ${wfu_last_gen}"
echo ""
echo "$(date +"%Y%m%d-%H:%M:%S")"
echo "user: ${user}"
echo ""

echo "skip_k: ${skip_k}"
echo "hsw_raw_ped: ${hsw_raw_ped}"

if [[ "$skip_k" == "false" ]]; then

    Rscript ${code}/wfu_exchange_breedHS_04_propose_pairs.R \
        --pop ${pop} \
        --n_pairs ${n_pairs} \
        --arg_parser ${arg_parser_04} \
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
        --utils ${utils} \
        --proj_dir ${proj_dir} \
        breedhs

elif [[ "$skip_k" == "true" ]]; then

    kinship_file=${proj_dir}/kinship/wfu_merged_gen${wfu_last_gen}_kinship_all_pairings.csv
    echo "Using kinship file: ${kinship_file}"
    echo ""

    Rscript ${code}/wfu_exchange_breedHS_04_propose_pairs.R \
        --pop ${pop} \
        --n_pairs ${n_pairs} \
        --arg_parser ${arg_parser_04} \
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
        --utils ${utils} \
        --proj_dir ${proj_dir} \
        --kinship_file ${kinship_file} \
        breedhs

fi

conda deactivate


END=$(date +%s)
echo "Time elapsed for breeder pairing: $(( $END - $START )) seconds"
