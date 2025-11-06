#!/bin/bash

# ==========================================================
# R0.93_bootstrap_recovery_step3
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
  > "$outputDir/remaining_reads.txt"

regen_ref="$outputDir/regen_ref.txt"
./bin/R0.93_encode "$dec_results2/recovery_bitstream.txt" "$regen_ref"


TypeIII_reads="$outputDir/Type-III_reads.txt"
echo "[Step 1]  Type-III reads filtering"
./bin/align_regen_ref  \
  "$regen_ref" \
  "$outputDir/remaining_reads.txt" \
  "$TypeIII_reads"


# ----------------------------------------------------------
# [Step 2] Forward-Backward Algorithm for Type-III Reads
# ----------------------------------------------------------
echo "[Step 2]  Forward-Backward Algorithm"
FB_output_S3="$outputDir/indelCorrect_output_S.txt"
mode=2 
./bin/indel_correct \
  "$TypeIII_reads" \
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
soft_info="$outputDir/soft_infor.txt"
./bin/multi-read_merging \
  "$FB_output" \
  "$soft_info" \
  "$WatermarkSeq" 


# ----------------------------------------------------------
# [Step 4] Third LDPC Decoding
# ----------------------------------------------------------
BitErrorRateFile="$outputDir/BER.txt"
encoded_bit="./configure/encoded_bit.txt"
EncodeBitLen=32000
./bin/CalBitError  \
    "$soft_info" \
    "$encoded_bit" \
    "$EncodeBitLen"  \
    "$BitErrorRateFile" 

echo "[Step 4]  LDPC decoding"
correctedBitStream="$outputDir/recovery_bitstream.txt"
./bin/R0.93_decode \
  "$soft_info" \
  "$correctedBitStream" 

OriginalBitstreamFile="./configure/source_bit.txt"
BitstreamLen=29760
hamming_error_after_dec="$outputDir/err_of_dec_result.txt"
./bin/post_dec_hamming_dis \
  "$correctedBitStream" \
  "$OriginalBitstreamFile" \
  "$hamming_error_after_dec" \
  "$BitstreamLen"
