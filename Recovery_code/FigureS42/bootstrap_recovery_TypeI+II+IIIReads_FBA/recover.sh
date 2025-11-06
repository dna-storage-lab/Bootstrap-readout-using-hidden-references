#!/bin/bash


echo "=================== Begin Processing ==================="
echo "Start Time: $(date)"
echo "========================================================"

start_time=$(date +%s.%3N)

# ----------------------------------------------------------
# Directory Setup
# ----------------------------------------------------------
mkdir -p ./results
result_log="./results/recovery_status.txt"
: > "$result_log"

# ----------------------------------------------------------
# Program Parameters
# ----------------------------------------------------------
WatermarkSeq="./configure/sequenceL80640NoPeriod2FILE"
ins=0.0026
del=0.0026
sub=0.0008
PayloadSize=40320
NumThreads=1 
MaxThreads=5
Job_Count=0
ExpNum=5


# ----------------------------------------------------------
# Input Files 
# ----------------------------------------------------------
InputFastq="../../../Sequencing_data/ART_simulated_data/err_0.006/DNA-40.32Kb-ER-Sim-1.fastq"       
ReadLength=150
CorrThreshold=108
NumReads=431 # 1.6x


# ----------------------------------------------------------
# Subsample Reads
# ----------------------------------------------------------
SubDir="./sequencing_data/sub_reads"
mkdir -p "$SubDir"
./bin/SubsampleFastqRandom "$InputFastq" "$NumReads" "$ExpNum" "$SubDir"

# ----------------------------------------------------------
# Run Processing for Each Experiment
# ----------------------------------------------------------
for Exp in $(seq 1 $ExpNum); do
(
    echo ""
    echo "========================================================"
    echo " Experiment #${Exp}"
    echo "--------------------------------------------------------"

    SubReadsFile="$SubDir/sub_reads_${Exp}.txt"
    type1count=0; type2count=0; type3count=0; Era=0; Err=0

    # Stage-I Recovery
    ExpDir_S1="./results/exp${Exp}/stepI_dec_result"
    mkdir -p "$ExpDir_S1"

    ./R0.5_bootstrap_recovery_step1.sh \
        "$SubReadsFile" \
        "$WatermarkSeq" \
        "$ExpDir_S1" \
        "$ReadLength" \
        "$NumReads" \
        "$ins" \
        "$del" \
        "$sub" \
        "$CorrThreshold" \
        "$NumThreads"

    err1=$(cut -d' ' -f1 "$ExpDir_S1/err_of_dec_result.txt")
    if [ "$err1" -eq 0 ]; then
        status=1
        type1count=$(wc -l < "$ExpDir_S1/Type-I_reads.txt")
        Era=$(awk 'END{print $3}' "$ExpDir_S1/BER.txt")
        Err=$(awk 'END{print $4}' "$ExpDir_S1/BER.txt")
    else
        # Stage-II Recovery
        ExpDir_S2="./results/exp${Exp}/stepII_dec_result"
        mkdir -p "$ExpDir_S2"

        ./R0.5_bootstrap_recovery_step2.sh \
            "$WatermarkSeq" \
            "$ExpDir_S2" \
            "$ReadLength" \
            "$ins" \
            "$del" \
            "$sub" \
            "$NumThreads" \
            "$ExpDir_S1"

        err2=$(cut -d' ' -f1 "$ExpDir_S2/err_of_dec_result.txt")
        if [ "$err2" -eq 0 ]; then
            status=2
            type1count=$(wc -l < "$ExpDir_S1/Type-I_reads.txt")
            type2count=$(wc -l < "$ExpDir_S2/Type-II_reads.txt")
            Era=$(awk 'END{print $3}' "$ExpDir_S2/BER.txt")
            Err=$(awk 'END{print $4}' "$ExpDir_S2/BER.txt")
        else
            # Stage-III Recovery
            ExpDir_S3="./results/exp${Exp}/stepIII_dec_result"
            mkdir -p "$ExpDir_S3"

            ./R0.5_bootstrap_recovery_step3.sh \
                "$SubReadsFile" \
                "$WatermarkSeq" \
                "$ExpDir_S3" \
                "$ReadLength" \
                "$ExpDir_S1" \
                "$ins" \
                "$del" \
                "$sub" \
                "$ExpDir_S2"

            err3=$(cut -d' ' -f1 "$ExpDir_S3/err_of_dec_result.txt")
            if [ "$err3" -eq 0 ]; then
                status=3
                type1count=$(wc -l < "$ExpDir_S1/Type-I_reads.txt")
                type2count=$(wc -l < "$ExpDir_S2/Type-II_reads.txt")
                type3count=$(wc -l < "$ExpDir_S3/Type-III_reads.txt")
                Era=$(awk 'END{print $3}' "$ExpDir_S3/BER.txt")
                Err=$(awk 'END{print $4}' "$ExpDir_S3/BER.txt")
            else
                status=0
                type1count=$(wc -l < "$ExpDir_S1/Type-I_reads.txt")
                type2count=$(wc -l < "$ExpDir_S2/Type-II_reads.txt")
                type3count=$(wc -l < "$ExpDir_S3/Type-III_reads.txt")
                Era=$(awk 'END{print $3}' "$ExpDir_S3/BER.txt")
                Err=$(awk 'END{print $4}' "$ExpDir_S3/BER.txt")
            fi
        fi
    fi

    # Write results
    case "$status" in
        1) echo "$Exp $Era $Err $status $type1count" >> "$result_log" ;;
        2) echo "$Exp $Era $Err $status $type1count $type2count" >> "$result_log" ;;
        3) echo "$Exp $Era $Err $status $type1count $type2count $type3count" >> "$result_log" ;;
        *) echo "$Exp $Era $Err $status $type1count $type2count $type3count" >> "$result_log" ;;
    esac

    # Calculate recovery rate
    awk 'BEGIN{success=0; total=0}
         {if ($4 != 0) success++; total++}
         END{printf "Recovered: %d/%d (%.2f%%)\n", success, total, (success/total)*100}' \
        "$result_log"

) &

    # Job control
    Job_Count=$((Job_Count + 1))
    if (( Job_Count >= MaxThreads )); then
        wait
        Job_Count=0
    fi
done

wait # Wait for all background jobs to complete

end_time=$(date +%s.%3N)
elapsed_time=$(echo "scale=3; $end_time - $start_time" | bc)

echo ""
echo "End Time: $(date)"
echo "Total Execution Time: ${elapsed_time} seconds"
echo "========================================================"