#!/bin/bash

# initialize variables
skip_k="false"

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
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown argument '$1'"
            usage
            ;;
    esac
done

usage() {
    echo "Usage: $0 --step <step_name> --args <path>"
    echo ""
    echo "Arguments:"
    echo "  --step    Required. Analysis step to run: 'make_id_maps', 'setup_peds', or 'propose_pairs'"
    echo "  --args    Required. Path to the breedHS arguments file"
    echo ""
    echo "Examples:"
    echo "  $0 --step setup_peds"
    echo "  $0 --step propose_pairs --skip_k"
    exit 1
}

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

if [[ ! "$step" =~ ^(make_id_maps|setup_peds|propose_pairs|replace_breeders)$ ]]; then
    echo "Error: Invalid step '$step'"
    usage
fi


source ${breedhs_args}
source activate ${conda_env}


if [[ "$step" == "make_id_maps" ]]; then
echo "Running breedHS step 1: ID maps"

    Rscript ${code}/wfu_breedHS_01_make_id_maps.R \
        --utils ${utils} \
        --arg_parser ${arg_parser_01} \
        --wfu_raw_ped ${wfu_raw_ped} \
        --all_wfu_ss ${all_wfu_shipping_sheets[@]} \
        --hsw_raw_ped ${hsw_raw_ped} \
        --all_hsw_ss ${all_hsw_shipping_sheets[@]} \
        --outdir ${id_map_dir} \
        breedhs
fi

if [[ "$step" == "setup_peds" ]]; then
echo "Running breedHS step 2: pedigree setup"

    Rscript ${code}/wfu_breedHS_02_setup_peds.R \
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
        --proj_dir ${proj_dir} \
        --arg_parser ${arg_parser_02} \
        breedhs

fi


if [[ "$step" == "propose_pairs" ]]; then
echo "Running breedHS step 3: proposing breeder pairs"

	job_name=wfu_gen${wfu_last_gen}_propose_pairs
	logs_dir=${proj_dir}/logs/03_propose_pairs
	mkdir -p ${logs_dir}

	sbatch -J ${job_name} -o ${logs_dir}/${job_name}-%j.o -e ${logs_dir}/${job_name}-%j.o \
		-p ${partition} -q ${qos} --account ${allocation} \
        -N 1 -c ${cpu} --mem-per-cpu ${mem} --time 1:00:00  \
		--export breedhs_args="${breedhs_args}",skip_k="${skip_k}" \
		${code}/wfu_breedHS_03_propose_pairs.sh

fi


conda deactivate
