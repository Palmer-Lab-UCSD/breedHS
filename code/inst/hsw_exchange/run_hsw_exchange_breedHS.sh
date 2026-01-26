#!/bin/bash

# initialize variables
breedhs_args=""
step=""
skip_k="false"
proposed_pairs="NULL"

usage() {
    echo "Usage: $0 --step <step_name> [--skip_k]"
    echo ""
    echo "Arguments:"
    echo " --args <filepath>      Required. File path to the breedHS arguments file"
    echo "  --step <step_name>    Required. Analysis step to run: 'make_id_maps', 'fams_to_ship', 'setup_peds', 'propose_pairs', or 'finalize_pairs'"
    echo "  --skip_k              Optional. Skip kinship calculations in propose_pairs step, read kinship from file instead"
    echo ""
    echo "Examples:"
    echo "  $0 --step setup_peds --args <filepath>"
    echo "  $0 --step propose_pairs --skip_k --args <filepath>"
    exit 1
}


# parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -a|--args)
            breedhs_args="$2"
            shift 2 
            ;;
        --step)
            step="$2"
            shift 2
            ;;
        --skip_k)
            skip_k="true"
            shift 1
            ;;
        --use_pairs)
            proposed_pairs="$2"
            shift 2 
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown argument '$1'"
            usage
            ;;
    esac
done

# validate required arguments
if [[ -z "$breedhs_args" ]]; then
    echo "Error: --args is required"
    usage
fi

# convert args to absolute path if relative
if [[ ! "$breedhs_args" = /* ]]; then
    breedhs_args="$(cd "$(dirname "$breedhs_args")" && pwd)/$(basename "$breedhs_args")"
fi

# check if args file exists
if [[ ! -f "$breedhs_args" ]]; then
    echo "Error: Arguments file not found: $breedhs_args"
    exit 1
fi

if [[ -z "$step" ]]; then
    echo "Error: --step argument is required"
    usage
fi

if [[ ! "$step" =~ ^(make_id_maps|setup_peds|propose_pairs|finalize_pairs)$ ]]; then
    echo "Error: Invalid step '$step'"
    usage
fi


source ${breedhs_args}
source activate ${conda_env}

if [[ "$step" == "make_id_maps" ]]; then
echo "Running breedHS step 1: ID maps"

	logs_dir=${proj_dir}/logs/01_make_id_maps
    mkdir -p ${logs_dir}
    timestamp=$(date +"%Y%m%d-%H:%M:%S")
    log_file=${logs_dir}/01_make_id_maps_${timestamp}.o

    echo "all_wfu_shipping_sheets: ${all_wfu_shipping_sheets[@]}"
    Rscript ${code}/hsw_exchange_breedHS_01_make_id_maps.R \
        --utils ${utils} \
        --arg_parser ${arg_parser_01} \
        --wfu_raw_ped ${wfu_raw_ped} \
        --all_wfu_ss ${all_wfu_shipping_sheets[@]} \
        --hsw_raw_ped ${hsw_raw_ped} \
        --hsw_assignments ${hsw_assignments} \
        --prev_colony_df ${prev_colony_df} \
        --all_hsw_ss ${all_hsw_shipping_sheets[@]} \
        --outdir ${id_map_dir} \
        breedhs 2>&1 | tee ${log_file}

fi


if [[ "$step" == "setup_peds" ]]; then
echo "Running breedHS step 2: pedigree setup"

	logs_dir=${proj_dir}/logs/02_setup_peds
    mkdir -p ${logs_dir}
    timestamp=$(date +"%Y%m%d-%H:%M:%S")
    log_file=${logs_dir}/02_setup_peds_${timestamp}.o

    Rscript ${code}/hsw_exchange_breedHS_02_setup_peds.R \
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
        --wfu_shipping_sheet ${wfu_shipping_sheet} \
        --hsw_assignments ${hsw_assignments} \
        --prev_colony_df ${prev_colony_df} \
        --wfu_id_map ${wfu_id_map} \
        --hsw_id_map ${hsw_id_map} \
        --proj_dir ${proj_dir} \
        --arg_parser ${arg_parser_02} \
        breedhs 2>&1 | tee ${log_file}

fi


if [[ "$step" == "propose_pairs" ]]; then
echo "Running breedHS step 3: proposing breeder pairs"

    if [[ "$skip_k" == "true" ]]; then

        kinship_file=${proj_dir}/kinship/hsw_merged_gen${hsw_last_gen}_kinship_all_pairings.csv
        echo "Using kinship file: ${kinship_file}"
        echo ""
    fi

    if [[ -f "$proposed_pairs" ]]; then
    echo "Using proposed pairs file: ${proposed_pairs}"
    fi

    echo "breedhs_args: ${breedhs_args}"
    echo "proj_dir: ${proj_dir}"

	job_name=hsw_gen${hsw_last_gen}_propose_pairs
	logs_dir=${proj_dir}/logs/03_propose_pairs
	mkdir -p ${logs_dir}

	sbatch -J ${job_name} -o ${logs_dir}/${job_name}-%j.o -e ${logs_dir}/${job_name}-%j.o \
		-p ${partition} -q ${qos} --account ${allocation} \
        -N 1 -c ${cpu} --mem-per-cpu ${mem} --time 1:00:00  \
		--export breedhs_args="${breedhs_args}",skip_k="${skip_k}",proposed_pairs="${proposed_pairs}" \
        --exclude "tscc-14-10" \
		${code}/hsw_exchange_breedHS_03_propose_pairs.sh

fi



if [[ "$step" == "finalize_pairs" ]]; then
echo "Running breedHS step 4: finalizing breeder pairs + final files"
    
    logs_dir=${proj_dir}/logs/04_final_pairs
	mkdir -p ${logs_dir}
    timestamp=$(date +"%Y%m%d-%H:%M:%S")
    log_file=${logs_dir}/hsw_gen${hsw_last_gen}_final_pairs_${timestamp}.o

    proposed_pairs=$( ls ${breedhs_out}/*breeders_proposed* )
    colony_pairs=$( ls ${colony_dir}/*breeders_proposed* )
    
    Rscript ${code}/hsw_exchange_breedHS_04_finalize_pairs.R \
        --pop ${pop} \
        --generation ${hsw_last_gen} \
        --hsw_raw_stem ${hsw_raw_stem} \
        --colony_dir ${colony_dir} \
        --peds_dir ${peds_dir} \
        --pairs_proposed ${proposed_pairs} \
        --colony_pairs ${colony_pairs} \
        --colony_df ${colony_df} \
        --prev_colony_df ${prev_colony_df} \
        --utils ${utils} \
        --proj_dir ${proj_dir} \
        --arg_parser ${arg_parser_04} \
        breedhs 2>&1 | tee ${log_file}

fi


conda deactivate