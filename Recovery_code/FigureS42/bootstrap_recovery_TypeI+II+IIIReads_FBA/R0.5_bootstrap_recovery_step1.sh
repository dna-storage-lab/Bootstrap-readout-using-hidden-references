#!/bin/bash


SubReadsFile=$1
WatermarkSeq=$2
outputDir=$3
ReadLength=$4   
NumReads=$5                   
ins=$6
del=$7
sub=$8
CorrThreshold=$9       
NumThreads=${10}



echo -e "\n\n========================================================"
echo "Start the first recovery process!"
echo "========================================================"

    # ----------------------------------------------------------
    # [Step 1]: Sliding Correlation
    # ----------------------------------------------------------
    # Inputs:
    #   • WatermarkSeq                  – known watermark sequence
    #   • SubReadsFile                  – Sequencing data
    # Outputs:
    #   • CorrelationResultsFile        – read alignment info
    CorrelationResultsFile="$outputDir/correlation_result.txt"
    echo "[Step 1] Sliding correlation"
    ./bin/SlidingCorrelation  \
        "$NumThreads" \
        "$ReadLength" \
        "$NumReads" \
        "$WatermarkSeq" \
        "$SubReadsFile" \
        "$CorrelationResultsFile"

    # Type-I Reads Filtering Based on Correlation Threshold
    # Inputs:
    #   • CorrelationResultsFile          – read alignment info
    # Outputs:
    #   • TypeI_reads                     – Type-I Reads
    #   • low_corr_reads                  – correlation-failed reads
    TypeI_reads="$outputDir/Type-I_reads.txt"
    low_corr_reads="$outputDir/lowthres_reads.txt"
    ./bin/getthre \
      "$CorrelationResultsFile" \
      "$NumReads" \
      "$CorrThreshold" \
      "$TypeI_reads" \
      "$low_corr_reads"


    # ----------------------------------------------------------
    # [Step 2] Forward-Backward Algorithm 
    # ----------------------------------------------------------
    # Inputs:
    #   • TypeI_reads                      –  from Step 1
    #   • WatermarkSeq                     – known watermark sequence
    # Outputs:
    #   • FB_output                        – indel corrected symbol probability
    echo "[Step 2]  Forward-Backward Algorithm "
    FB_output="$outputDir/symbol_probability.txt"
    mode=0 # process Type-I reads  
    ./bin/R1_2_indel_correct  \
      "$TypeI_reads" \
      "$FB_output" \
      "$WatermarkSeq" \
      "$ins" \
      "$sub" \
      "$del" \
      "$mode"


    # ----------------------------------------------------------
    # [Step 3] Consensus soft information generation
    # ----------------------------------------------------------
    # Inputs:
    #   • FB_output                        –  from Step 2
    #   • WatermarkSeq                     – known watermark sequence
    # Outputs:
    #   • soft_Info                        – consensus soft information
    echo "[Step 3]  Consensus soft information generation"
    soft_Info="$outputDir/soft_info.txt"
    ./bin/R1_2_multi-read_merging \
      "$FB_output" \
      "$soft_Info" \
      "$WatermarkSeq" 


    # ----------------------------------------------------------
    # [Step 4] LDPC Decoding
    # ----------------------------------------------------------
    # Inputs:
    #   • soft_Info                          –  from Step 3
    #   • encoded_bit                        –  ground‑truth encoded bitstream
    # Outputs:
    #   • BitErrorRateFile                   – BER before decoding 
    BitErrorRateFile="$outputDir/BER.txt"
    encoded_bit="./configure/encoded_bit.txt"
    EncodeBitLen=64512
    ./bin/CalBitError  \
        "$soft_Info" \
        "$encoded_bit" \
        "$EncodeBitLen"  \
        "$BitErrorRateFile" 

    # Inputs:
    #   • soft_Info                           –  from Step 3
    # Outputs:
    #   • recovery_image                      – recovery image
    #   • correctedBitStream                  – decoded bitstream
    #   • dec_cw                              – decoded codeword
    echo "[Step 4]  LDPC decoding"
    correctedBitStream="$outputDir/recovery_bitstream.txt"
    recovery_image="$outputDir/recovery_image.jpg"
    dec_cw="$outputDir/decodedCodeword.txt"
    ./bin/NLDPC_r1_2_soft_decoder \
      "$soft_Info" \
      "$correctedBitStream" \
      "$recovery_image" \
      S \
      "$dec_cw"

    # Inputs:
    #   • correctedBitStream                  – decoded bitstream
    #   • OriginalBitstream                   – ground‑truth bitstream
    # Outputs:
    #   • hamming_error_after_dec             – hamming distance
    OriginalBitstream="./configure/source_bit.txt"
    BitstreamLen=32256
    hamming_error_after_dec="$outputDir/err_of_dec_result.txt"
    ./bin/post_dec_hamming_dis \
      "$correctedBitStream" \
      "$OriginalBitstream" \
      "$hamming_error_after_dec" \
      "$BitstreamLen"

