#!/bin/bash

START=$(date +%s)
source ${breedhs_args}
source activate ${conda_env}

echo ""
echo "--------------------------------------------------------"
echo "--------------- HSW WEST BREEDER PAIRING ---------------"
echo "----------- (AND EXCHANGE WITH WAKE FOREST) ------------"
echo "--------------------------------------------------------"
echo ""
echo "Identifying HSW gen${hsw_ship_gen} families to pair with WFU gen${wfu_exchange_gen}"
echo ""
echo "$(date +"%Y%m%d-%H:%M:%S")"
echo "user: ${user}"
echo ""


if [[ "$skip_k" == "false" ]]; then

    Rscript ${code}/hsw_exchange_breedHS_02_fams_to_ship.R \
        --pop ${pop} \
        --utils ${utils} \
        --peds_dir ${peds_dir} \
        --wfu_raw_ped ${wfu_raw_ped} \
        --wfu_first_gen ${wfu_first_gen} \
        --wfu_last_gen ${wfu_last_gen} \
        --wfu_stem ${wfu_stem} \
        --wfu_raw_stem ${hsw_raw_stem} \
        --hsw_raw_ped ${hsw_raw_ped} \
        --hsw_first_gen ${hsw_first_gen} \
        --hsw_last_gen ${hsw_last_gen} \
        --hsw_stem ${hsw_stem} \
        --hsw_raw_stem ${hsw_raw_stem} \
        --wfu_id_map ${wfu_id_map} \
        --hsw_id_map ${hsw_id_map} \
        --ped_map ${ped_map} \
        --hsw_to_wfu ${hsw_to_wfu} \
        --wfu_to_hsw ${wfu_to_hsw} \
        --hsw_id_fams_from ${hsw_id_fams_from} \
        --wfu_id_fams_from ${wfu_id_fams_from} \
        --proj_dir ${proj_dir} \
        --colony_df ${colony_df} \
        --arg_parser ${arg_parser_02} \
        breedhs

elif [[ "$skip_k" == "true" ]]; then

    kinship_file=${proj_dir}/kinship/wfu_merged_gen${wfu_last_gen}_kinship_all_pairings.csv
    echo "Using kinship file: ${kinship_file}"
    echo ""

    Rscript ${code}/hsw_exchange_breedHS_02_fams_to_ship.R \
        --pop ${pop} \
        --utils ${utils} \
        --peds_dir ${peds_dir} \
        --wfu_raw_ped ${wfu_raw_ped} \
        --wfu_first_gen ${wfu_first_gen} \
        --wfu_last_gen ${wfu_last_gen} \
        --wfu_stem ${wfu_stem} \
        --wfu_raw_stem ${hsw_raw_stem} \
        --hsw_raw_ped ${hsw_raw_ped} \
        --hsw_first_gen ${hsw_first_gen} \
        --hsw_last_gen ${hsw_last_gen} \
        --hsw_stem ${hsw_stem} \
        --hsw_raw_stem ${hsw_raw_stem} \
        --wfu_id_map ${wfu_id_map} \
        --hsw_id_map ${hsw_id_map} \
        --ped_map ${ped_map} \
        --hsw_to_wfu ${hsw_to_wfu} \
        --wfu_to_hsw ${wfu_to_hsw} \
        --hsw_id_fams_from ${hsw_id_fams_from} \
        --wfu_id_fams_from ${wfu_id_fams_from} \
        --proj_dir ${proj_dir} \
        --colony_df ${colony_df} \
        --arg_parser ${arg_parser_02} \
        --kinship_file ${kinship_file} \
        breedhs

fi

conda deactivate


END=$(date +%s)
echo "Time elapsed for breeder pairing: $(( $END - $START )) seconds"
