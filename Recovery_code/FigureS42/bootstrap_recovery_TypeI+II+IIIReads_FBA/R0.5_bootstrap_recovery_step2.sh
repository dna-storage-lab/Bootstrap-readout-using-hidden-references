#!/bin/bash

# ==========================================================
# R0.83_bootstrap_recovery_step2
# ==========================================================

WatermarkSeq=$1
outputDir=$2
ReadLength=$3      
ins=$4
del=$5
sub=$6
NumThreads=$7
dec_results1=$8       



echo -e "\n\n========================================================"
echo "First decoding failed. Start the second recovery process!"
echo "========================================================"
# ----------------------------------------------------------
# [Step 1] Majority Voting: Generate Scaffold Reference
# ----------------------------------------------------------
echo "[Step 1] Generate Scaffold Reference"
./bin/majorityvoting "$dec_results1" "$outputDir"


# ----------------------------------------------------------
# [Step 2] Edlib Alignment for Correlation-failed Reads
# ----------------------------------------------------------
low_corr_reads="$dec_results1/lowthres_reads.txt"
TypeII_reads="$outputDir/Type-II_reads.txt"
echo "[Step 2]  Type-II reads filtering"
./bin/lowthres_pthread_edlib \
  "$NumThreads" \
  "$outputDir/scaffold_ref.txt" \
  "$low_corr_reads" \
  "$TypeII_reads"


# ----------------------------------------------------------
# [Step 3] Forward-Backward Algorithm
# ----------------------------------------------------------
echo "[Step 3]  Forward-Backward Algorithm"
FB_output_S1="$dec_results1/symbol_probability.txt"
FB_output_S2="$outputDir/indelCorrect_output_S.txt"
mode=1 # process Type-II reads  
./bin/R1_2_indel_correct \
  "$TypeII_reads" \
  "$FB_output_S2" \
  "$WatermarkSeq" \
  "$ins" \
  "$sub" \
  "$del" \
  "$mode"


# ----------------------------------------------------------
# [Step 4] Consensus soft information generation
# ----------------------------------------------------------
FB_output="$outputDir/symbol_probability.txt"
cat "$FB_output_S1" "$FB_output_S2" > "$FB_output"

echo "[Step 4]  Consensus soft information generation"
soft_info="$outputDir/soft_info.txt"
./bin/R1_2_multi-read_merging \
  "$FB_output" \
  "$soft_info" \
  "$WatermarkSeq"



# ----------------------------------------------------------
# [Step 5]  LDPC Decoding
# ----------------------------------------------------------
# BER
BitErrorRateFile="$outputDir/BER.txt"
encoded_bit="./configure/encoded_bit.txt"
EncodeBitLen=64512
./bin/CalBitError  \
    "$soft_info" \
    "$encoded_bit" \
    "$EncodeBitLen"  \
    "$BitErrorRateFile" 
  
echo "[Step 5]  LDPC decoding"
correctedBitStream="$outputDir/recovery_bitstream.txt"
recovery_image="$outputDir/recovery_image.jpg"
dec_results2="$outputDir/decodedCodeword.txt"
./bin/NLDPC_r1_2_soft_decoder \
  "$soft_info" \
  "$correctedBitStream" \
  "$recovery_image" \
  S \
  "$dec_results2"

OriginalBitstreamFile="./configure/source_bit.txt"
BitstreamLen=32256
hamming_error_after_dec="$outputDir/err_of_dec_result.txt"
./bin/post_dec_hamming_dis \
  "$correctedBitStream" \
  "$OriginalBitstreamFile" \
  "$hamming_error_after_dec" \
  "$BitstreamLen"
