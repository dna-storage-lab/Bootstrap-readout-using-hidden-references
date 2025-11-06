#!/bin/bash

# ==========================================================
# R0.5_bootstrap_recovery_step3
# ==========================================================

SubReadsFile=$1
WatermarkSeq=$2
outputDir=$3
ReadLength=$4   
dec_results1=$5                   
ins=$6
del=$7
sub=$8
dec_results2=$9 



echo -e "\n\n========================================================"
echo "Second decoding failed. Start the third recovery process!"
echo "========================================================"
# ----------------------------------------------------------
# [Step 1] Decode Result Feedback to Generate Regenerative Reference
# ----------------------------------------------------------

# Filter remaining low-threshold reads not recovered in Step II
awk 'NR==FNR{a[$1]; next} !($1 in a)' \
  $dec_results2/Type-II_reads.txt \
  $dec_results1/lowthres_reads.txt \
  > "$outputDir/scaffold_unaligned_reads.txt"

# align to regenerative ref
echo "[Step 1]  Type-III reads filtering"
./bin/decode_feedback_align  \
  "$dec_results2/decodedCodeword.txt" \
  "$outputDir/scaffold_unaligned_reads.txt" \
  "$outputDir/Type-III_reads.txt"



# ----------------------------------------------------------
# [Step 2] Forward-Backward Algorithm for Type-III Reads
# ----------------------------------------------------------
echo "[Step 2]  Forward-Backward Algorithm"
FB_output_S3="$outputDir/indelCorrect_output_S.txt"
mode=2 # process Type-III reads  
./bin/R1_2_indel_correct \
  "$outputDir/Type-III_reads.txt" \
  "$FB_output_S3" \
  "$WatermarkSeq" \
  "$ins" \
  "$sub" \
  "$del" \
  "$mode"


# ----------------------------------------------------------
# [Step 3] Consensus soft information generation
# ----------------------------------------------------------
FB_output_S2="$dec_results2/symbol_probability.txt"
FB_output="$outputDir/symbol_probability.txt"
cat "$FB_output_S2" "$FB_output_S3" > "$FB_output"

echo "[Step 3]  Consensus soft information generation"
soft_info="$outputDir/soft_info.txt"
./bin/R1_2_multi-read_merging \
  "$FB_output" \
  "$soft_info" \
  "$WatermarkSeq" \


# ----------------------------------------------------------
# [Step 4] Third LDPC Decoding
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

echo "[Step 4]  LDPC decoding"
S_correctedBitStream3="$outputDir/recovery_bitstream.txt"
recovery_image="$outputDir/recovery_image.jpg"
dec_results3="$outputDir/decodedCodeword.txt"
./bin/NLDPC_r1_2_soft_decoder \
  "$soft_info" \
  "$S_correctedBitStream3" \
  "$recovery_image" \
  S \
  "$dec_results3"

OriginalBitstreamFile="./configure/source_bit.txt"
BitstreamLen=32256
hamming_error_after_dec3="$outputDir/err_of_dec_result.txt"
./bin/post_dec_hamming_dis \
  "$S_correctedBitStream3" \
  "$OriginalBitstreamFile" \
  "$hamming_error_after_dec3" \
  "$BitstreamLen"
