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
FB_output_error2="$outputDir/indel_corrected_error.txt"
mode=1
./bin/R5_6_indel_correct \
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
outputfilecat2="$outputDir/symbol_probability.txt"
cat "$FB_output_S1" "$FB_output_S2" > "$outputfilecat2"

echo "[Step 4]  Consensus soft information generation"
outputfile2="$outputDir/soft_info.txt"
./bin/R5_6_multi-read_merging \
  "$outputfilecat2" \
  "$outputfile2" \
  "$WatermarkSeq"



# ----------------------------------------------------------
# [Step 5]  LDPC Decoding
# ----------------------------------------------------------
# BER
BitErrorRateFile="$outputDir/BER.txt"
encoded_bit="./configure/encoded_bit.txt"
EncodeBitLen=64800
./bin/CalBitError  \
    "$outputfile2" \
    "$encoded_bit" \
    "$EncodeBitLen"  \
    "$BitErrorRateFile" 
  
echo "[Step 5]  LDPC decoding"
S_correctedBitStream2="$outputDir/recovery_bitstream.txt"
recovery_image2="$outputDir/recovery_image.jpg"
dec_results2="$outputDir/decodedCodeword.txt"
./bin/LDPC_r5_6_soft_decoder \
  "$outputfile2" \
  "$S_correctedBitStream2" \
  "$recovery_image2" \
  S \
  "$dec_results2"

OriginalBitstreamFile="./configure/source_bit.txt"
BitstreamLen=54000
hamming_error_after_dec2="$outputDir/err_of_dec_result.txt"
./bin/post_dec_hamming_dis \
  "$S_correctedBitStream2" \
  "$OriginalBitstreamFile" \
  "$hamming_error_after_dec2" \
  "$BitstreamLen"
