#!/usr/bin/env bash
set -euo pipefail

echo_usage() {
    echo "Usage: $0 <read_fastq> [options]"
    echo
    echo "Required:"
    echo "  <read_fastq>               Input FASTQ file for long-read assembly"
    echo
    echo "General options:"
    echo "  -o, --outdir <path>         Output directory [default: current directory]"
    echo "  -t, --threads <int>         Threads per assembly job [default: 8]"
    echo "  -r, --read-type <type>      Read type: ont_r9 | ont_r10 | pacbio_clr | pacbio_hifi [default: ont_r10]"
    echo "  -k, --keep-intermediate     Keep subsampled and filtered FASTQ files after assembly"
    echo "  -w, --overwrite             Delete output directory if it exists (use with caution!)"
    echo "  -T, --max-time <duration>   Slurm time limit per job [default: 8h]"
    echo "  -M, --max-mem <size>        Slurm memory per job [default: 32g]"
    echo "  -R, --resume-after-assembly Resume pipeline after assembly step (skip read filtering, subsampling, and assembly)"
    echo
    echo "Read filtering with filtlong:"
    echo "  -l, --min-length <int>      Filter out reads shorter than this length"
    echo "  -b, --target-bases <int>    Keep top reads until total base count is reached"
    echo "  -p, --keep-percent <float>  Keep best X%% of reads (e.g. 90)"
    echo
    echo "Help:"
    echo "  -h, --help                  Show this help message and exit"
}
# ------------------------------
# Parse args
# ------------------------------
reads=""
threads="8"
read_type="ont_r10"
keep_intermediate=false
max_time="8h"
max_mem="32g"
min_length=""
target_bases=""
keep_percent=""
outdir="."
overwrite=false
resume_after_assembly=false

while [[ $# -gt 0 ]]; do
    case "$1" in
    -o | --outdir)
        outdir="$2"
        shift 2
        ;;
    -t | --threads)
        threads="$2"
        shift 2
        ;;
    -r | --read-type)
        read_type="$2"
        shift 2
        ;;
    -R | --resume-after-assembly)
        resume_after_assembly=true
        shift
        ;;
    -k | --keep-intermediate)
        keep_intermediate=true
        shift
        ;;
    -w | --overwrite)
        overwrite=true
        shift
        ;;
    -T | --max-time)
        max_time="$2"
        shift 2
        ;;
    -M | --max-mem)
        max_mem="$2"
        shift 2
        ;;
    -l | --min-length)
        min_length="$2"
        shift 2
        ;;
    -b | --target-bases)
        target_bases="$2"
        shift 2
        ;;
    -p | --keep-percent)
        keep_percent="$2"
        shift 2
        ;;
    -h | --help)
        echo_usage
        exit 0
        ;;
    -*)
        echo "[ERROR] Unknown option: $1" >&2
        echo "Use --help to see usage." >&2
        exit 1
        ;;
    *)
        if [[ -z "$reads" ]]; then
            reads=$(realpath "$1")
        else
            echo "[ERROR] Unexpected positional argument: $1" >&2
            echo "Use --help to see usage." >&2
            exit 1
        fi
        shift
        ;;
    esac
done

if ((threads > 128)); then threads=128; fi
case $read_type in
ont_r9 | ont_r10 | pacbio_clr | pacbio_hifi) ;;
*)
    echo "[ERROR] read_type must be ont_r9, ont_r10, pacbio_clr or pacbio_hifi" >&2
    exit 1
    ;;
esac

if [[ "$overwrite" == true && "$resume_after_assembly" == true ]]; then
    echo "[ERROR] Cannot use --overwrite and --resume-after-assembly together." >&2
    echo "        Overwrite would delete assemblies/ directory needed for resume." >&2
    exit 1
fi

if [[ -d "$outdir" && "$overwrite" == true ]]; then
    echo "[WARN] Output directory '$outdir' exists and will be overwritten."
    rm -rf "$outdir"
fi

mkdir -p "$outdir"
echo "[INFO] Changing to output directory: $outdir"
cd "$outdir"

if [[ "$resume_after_assembly" == false ]]; then
    # ------------------------------
    # Optional: Filter reads with filtlong
    # ------------------------------
    if [[ -n "$min_length" || -n "$target_bases" || -n "$keep_percent" ]]; then
        echo "[INFO] Filtering reads with filtlong..."
        mkdir -p filtered

        filtlong_cmd=(filtlong)
        [[ -n "$min_length" ]] && filtlong_cmd+=(--min_length "$min_length")
        [[ -n "$target_bases" ]] && filtlong_cmd+=(--target_bases "$target_bases")
        [[ -n "$keep_percent" ]] && filtlong_cmd+=(--keep_percent "$keep_percent")

        # Strip known FASTQ extensions and gz if present
        basename_no_ext=$(basename "$reads" | sed -E 's/(\.fastq|\.fq)(\.gz)?$//')
        reads_filtered="filtered/${basename_no_ext}_filtered.fastq"

        "${filtlong_cmd[@]}" "$reads" >"$reads_filtered"
        reads="$reads_filtered"

        echo "[INFO] Filtered reads written to: $reads"
    fi

    # ------------------------------
    # Estimate genome size
    # ------------------------------
    echo "[INFO] Estimating genome size..."
    genome_size=$(lrge -t "$threads" "$reads" 2>autocycler.stderr)

    # ------------------------------
    # Subsample reads
    # ------------------------------
    echo "[INFO] Subsampling reads..."
    autocycler subsample \
        --reads "$reads" \
        --out_dir subsampled_reads \
        --genome_size "$genome_size" \
        2>>autocycler.stderr

    # ------------------------------
    # Submit assembly jobs to Slurm
    # ------------------------------
    mkdir -p assemblies/slurm_logs
    rm -f assemblies/job_ids.txt

    echo "[INFO] Submitting assembly jobs..."
    for assembler in raven myloasm miniasm flye metamdbg necat nextdenovo plassembler canu; do
        for i in 01 02 03 04; do
            sample="subsampled_reads/sample_${i}.fastq"
            out_prefix="assemblies/${assembler}_${i}"
            cmd="autocycler helper $assembler --reads $sample --out_prefix $out_prefix --threads $threads --genome_size $genome_size --read_type $read_type --min_depth_rel 0.1"
            jobname="ac_${assembler}_${i}"

            # Submit and extract job ID
            jobid=$(ssubmit \
                -t "$max_time" \
                -m "$max_mem" \
                -o "assemblies/slurm_logs/${jobname}.out" \
                -e "assemblies/slurm_logs/${jobname}.err" \
                "$jobname" \
                "$cmd" -- -c "$threads" 2>&1 |
                grep -oP 'Submitted batch job \K[0-9]+')

            if [[ -n "$jobid" ]]; then
                echo "[INFO] Submitted $jobname → JobID=$jobid"
                echo "$jobid" >>assemblies/job_ids.txt
            else
                echo "[ERROR] Failed to get job ID for $jobname" >&2
                exit 1
            fi
        done
    done

    # ------------------------------
    # Wait for all jobs to complete, then check status
    # ------------------------------
    echo "[INFO] Waiting for assembly jobs to complete..."
    sleep 10 # Give Slurm time to register jobs

    start_time=$(date +%s)
    mapfile -t job_ids <assemblies/job_ids.txt

    while true; do
        sleep 60

        still_running=()
        for jobid in "${job_ids[@]}"; do
            if squeue -j "$jobid" 2>/dev/null | grep -q "$jobid"; then
                still_running+=("$jobid")
            fi
        done

        if ((${#still_running[@]} == 0)); then
            echo "[INFO] All assembly jobs exited squeue. Verifying final status..."
            break
        fi

        elapsed=$(($(date +%s) - start_time))
        minutes=$((elapsed / 60))
        seconds=$((elapsed % 60))
        echo "[INFO] ${#still_running[@]} job(s) still running after ${minutes}m ${seconds}s: ${still_running[*]}"
    done

    # Check final state with sacct
    fail_count=0
    while read -r jobid; do
        status=$(sacct -j "$jobid" --format=State --noheader | awk '{print $1}' | head -n 1)
        if [[ "$status" != "COMPLETED" ]]; then
            echo "[ERROR] Job $jobid failed with status: $status" >&2
            ((fail_count++))
        fi
    done <assemblies/job_ids.txt

    if ((fail_count > 0)); then
        echo "[FATAL] $fail_count job(s) failed. Exiting." >&2
        exit 1
    fi

    echo "[INFO] All assembly jobs completed successfully."
else
    echo "[INFO] Resuming pipeline after assembly."
    if [[ ! -d assemblies ]]; then
        echo "[ERROR] Cannot resume: assemblies/ directory not found." >&2
        exit 1
    fi
fi # end resume_after_assembly check

# ------------------------------
# Adjust output weights
# ------------------------------
echo "[INFO] Adjusting assembly outputs..."
shopt -s nullglob
for f in assemblies/plassembler*.fasta; do
    sed -i 's/circular=True/circular=True Autocycler_cluster_weight=3/' "$f"
done
for f in assemblies/canu*.fasta assemblies/flye*.fasta; do
    sed -i 's/^>.*$/& Autocycler_consensus_weight=2/' "$f"
done
shopt -u nullglob

# ------------------------------
# Cleanup and continue pipeline
# ------------------------------
if [[ "$keep_intermediate" = false ]]; then
    echo "[INFO] Removing intermediate reads..."
    rm -rf subsampled_reads/*.fastq
    rm -f filtered/*.fastq 2>/dev/null || true
fi

echo "[INFO] Running Autocycler compression, clustering, and resolution..."
autocycler compress -i assemblies -a autocycler_out 2>>autocycler.stderr
autocycler cluster -a autocycler_out 2>>autocycler.stderr

for c in autocycler_out/clustering/qc_pass/cluster_*; do
    autocycler trim -c "$c" 2>>autocycler.stderr
    autocycler resolve -c "$c" 2>>autocycler.stderr
done

autocycler combine \
    -a autocycler_out \
    -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa \
    2>>autocycler.stderr

echo "[INFO] Autocycler pipeline completed successfully."
