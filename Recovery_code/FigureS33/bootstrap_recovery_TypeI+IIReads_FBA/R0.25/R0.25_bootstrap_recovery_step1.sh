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

    # NOTE: Only I/O file paths are listed in Inputs/Outputs for brevity.

    # ----------------------------------------------------------
    # [Step 1]: Sliding Correlation
    # ----------------------------------------------------------
    # Inputs:
    #   • WatermarkSeq                     – known watermark sequence
    #   • SubReadsFile                     – Sequencing data
    # Outputs:
    #   • CorrelationResultsFile           – read alignment info
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
    #   • CorrelationResultsFile           – read alignment info
    # Outputs:
    #   • TypeI_Reads                      – Type-I Reads
    #   • low_corr_reads                   – correlation-failed reads
    TypeI_Reads="$outputDir/Type-I_reads.txt"
    low_corr_reads="$outputDir/lowthres_reads.txt"
    ./bin/getthre \
      "$CorrelationResultsFile" \
      "$NumReads" \
      "$CorrThreshold" \
      "$TypeI_Reads" \
      "$low_corr_reads"



    # ----------------------------------------------------------
    # [Step 2] Forward-Backward Algorithm 
    # ----------------------------------------------------------
    # Inputs:
    #   • TypeI_Reads                      –  from Step 1
    #   • WatermarkSeq                     – known watermark sequence
    # Outputs:
    #   • FB_output                        – indel corrected symbol probability
    FB_output="$outputDir/symbol_probability.txt"
    mode=0  # process Type-I reads
    echo "[Step 2]  Forward-Backward Algorithm "
    ./bin/indel_correct  \
      "$TypeI_Reads" \
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
    #   • FB_output                         –  from Step 2
    #   • WatermarkSeq                      – known watermark sequence
    # Outputs:
    #   • soft_info                         – consensus soft information
    echo "[Step 3]  Consensus soft information generation"
    soft_info="$outputDir/soft_info.txt"
    ./bin/multi-read_merging \
      "$FB_output" \
      "$soft_info" \
      "$WatermarkSeq" 



    # ----------------------------------------------------------
    # [Step 4] LDPC Decoding
    # ----------------------------------------------------------
    # Inputs:
    #   • soft_info                          –  from Step 3
    #   • encoded_bit                        –  ground‑truth encoded bitstream
    # Outputs:
    #   • BitErrorRateFile                   – BER before decding 
    BitErrorRateFile="$outputDir/BER.txt"
    encoded_bit="./configure/encoded_bit.txt"
    EncodeBitLen=64800
    ./bin/CalBitError  \
        "$soft_info" \
        "$encoded_bit" \
        "$EncodeBitLen"  \
        "$BitErrorRateFile" 

    # Inputs:
    #   • soft_info                          –  from Step 3
    # Outputs:
    #   • correctedBitStream                 – decoded bitstream
    #   • recovery_poem                      – recovery text
    #   • dec_cw                             – decoded codeword
    echo "[Step 4]  LDPC decoding"
    correctedBitStream="$outputDir/recovery_bitstream.txt"
    recovery_poem="$outputDir/recovery_poem.txt"
    dec_cw="$outputDir/decodedCodeword.txt"
    ./bin/LDPC_r1_4_soft_decoder \
      "$soft_info" \
      "$correctedBitStream" \
      "$recovery_poem" \
      S \
      "$dec_cw"

    # Inputs:
    #   • correctedBitStream                 – decoded bitstream
    #   • OriginalBitstream                  – ground‑truth bitstream
    # Outputs:
    #   • hamming_error_after_dec            – hamming distance
    OriginalBitstream="./configure/source_bit.txt"
    BitstreamLen=16200
    hamming_error_after_dec="$outputDir/err_of_dec_result.txt"
    ./bin/post_dec_hamming_dis \
      "$correctedBitStream" \
      "$OriginalBitstream" \
      "$hamming_error_after_dec" \
      "$BitstreamLen"
