#!/bin/bash


echo "=================== Begin Processing ==================="
echo "Start Time: $(date)"
echo "========================================================"
echo ""

StartTime=$(date +%s.%3N)

# ----------------------------------------------------------
# Directory Setup
# ----------------------------------------------------------
mkdir -p ./results
cd ./src
SubDir="../sequencing_data/sub_reads"
mkdir -p "$SubDir"
# ----------------------------------------------------------
# Program  Parameters
# ----------------------------------------------------------
WatermarkSeq="../configure/SequenceLengthALL_FILE001R0667"            # Watermark sequence 
NumThreads=5                    
ExpNum=5                                                              # Number of replicate experiments

# ----------------------------------------------------------
# Input Files 
# ----------------------------------------------------------
# SE150
InputFastq="../../../../Sequencing_data/Pooled_sequencing_data/DNA-40.5Kb-Mix-SE150.fastq"    # 150‑bp single‑end reads 
ReadLength=150
CorrThreshold=96                                                      # Correlation threshold
NumReadsList="864"  #  3.2x



# ----------------------------------------------------------
# Run Processing for Each Coverage Level
# ----------------------------------------------------------
for NumReads in "${NumReadsList[@]}"
do
    coverage_results="../results/$NumReads/"
    mkdir -p "$coverage_results"

    ../bin/SubsampleFastqRandom  "$InputFastq" "$NumReads" "$ExpNum" "$SubDir"

    # ----------------------------------------------------------
    # Run Processing for Each Experiment
    # ----------------------------------------------------------
    for Exp in $(seq 1 $ExpNum)
    do
        echo ""
        echo "========================================================"
        echo " Experiment #${Exp}"
        echo "--------------------------------------------------------"

        ExpDir="$coverage_results/exp${Exp}"
        mkdir -p "$ExpDir"

        SubReadsFile="../sequencing_data/sub_reads/sub_reads_${Exp}.txt"

        CorrelationResultsFile="${ExpDir}/correlation_result.txt"
        SoftInfoFile="${ExpDir}/soft_info.txt"
        DecodedBitstreamFile="${ExpDir}/recovery_bitstream.txt"
        OutputImageFile="${ExpDir}/recovery_image.jpg"
        PreLDPC_Error="$coverage_results/BER.txt"
        PostLDPC_Error="$coverage_results/err_of_dec_result.txt"

        # NOTE: Only I/O file paths are listed in Inputs/Outputs for brevity.

        # ------------------------------------------------------
        # [Step 1] Sliding Correlation
        # ------------------------------------------------------
        # Inputs:
        #   • WatermarkSeq                     – known watermark sequence
        #   • SubReadsFile                     – Sequencing data
        # Outputs:
        #   • CorrelationResultsFile           – read alignment info
        echo "[Step 1] Sliding correlation"
        ../bin/SlidingCorrelation  \
            "$NumThreads" \
            "$ReadLength" \
            "$NumReads" \
            "$WatermarkSeq" \
            "$SubReadsFile" \
            "$CorrelationResultsFile"

        # ------------------------------------------------------
        # [Step 2] Bit‑wise Majority Voting
        # ------------------------------------------------------
        # Inputs:
        #   • CorrelationResultsFile          – from Step 1
        #   • WatermarkSeq                    – known watermark sequence
        # Outputs:
        #   • SoftInfoFile                    – consensus soft information
        echo "[Step 2] Majority voting"
        ../bin/BitwiseConsensusRecovery  \
            "$CorrelationResultsFile" \
            "$CorrThreshold" \
            "$WatermarkSeq" \
            "$SoftInfoFile"


        # ------------------------------------------------------
        # [Step 3] LDPC Decoding
        # ------------------------------------------------------
        # Inputs:
        #   • SoftInfoFile                    – from Step 2
        #   • EncodedFile                     – ground‑truth encoded bitstream
        # Outputs:
        #   • PreLDPC_Error                   – BER before decoding
        EncodedFile="../configure/encoded_bit.txt"
        EncodeBitLen=64800
        ../bin/CalBitError  \
            "$SoftInfoFile" \
            "$EncodedFile" \
            "$EncodeBitLen"  \
            "$PreLDPC_Error" 

        # Inputs:
        #   • SoftInfoFile                    – from Step 2
        # Outputs:
        #   • OutputImageFile                 – reconstructed image 
        #   • DecodedBitstreamFile            – binary bitstream
        echo "[Step 3] LDPC decoding"
        ../bin/LDPC_r2_3_soft_decoder \
            "$SoftInfoFile" \
            "$DecodedBitstreamFile" \
            "$OutputImageFile" \
            "S"  # Soft-Decision Decoding

        # Inputs:
        #   • DecodedBitstreamFile            – from Step 3
        #   • OriginalBitstreamFile           – ground‑truth bitstream
        # Outputs:
        #   • PostLDPC_Error                  – hamming distance
        OriginalBitstreamFile="../configure/source_bit.txt"
        BitstreamLen=43200
        ../bin/PostDecHammingDistance  \
            "$DecodedBitstreamFile" \
            "$OriginalBitstreamFile" \
            "$PostLDPC_Error" \
            "$BitstreamLen" 


        # Result statistics
        awk '{ if ($1 == 0) success++; total++ } END { printf "Recovered: %d\nRecoveryRate: %.4f\n", success, success/total }' "$PostLDPC_Error"
        echo "========================================================"
    done

done

echo ""
echo "End Time: $(date)"