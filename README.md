# Fast bootstrap and reliable readout using hidden references for DNA data storage

![Alt Text](./image/image.png)

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Example of usage](#example-of-usage)
  - [1. Fast recovery (Illumina data, R=2/3)](#1-fast-recovery)
  - [2. Bootstrap recovery (ONT data, R=2/3)](#2-bootstrap-recovery)
- [Note](#note)
- [License](#license)

## Overview

Synthetic DNA is becoming a promising data storage medium for future large-scale data archiving. However, data readout from massive, unordered sequencing reads requires alignment based on overlapping regions and is complicated by diverse sequencing errors. We propose a multi-stage alignment and error correction strategy via multiple-fold hidden references, transforming the de novo readout into a resequencing-like workflow. We provide code for readout pipelines under different error conditions, divided into two main parts：

1. **Fast recovery**: In low-error-rate scenarios, the pipeline identifies reads via sliding correlation to watermark reference. Bit-wise consensus rapidly generates soft-decision information for LDPC decoding.
2. **Bootstrap recovery**: In the presence of indels, the pipeline progressively identifies reads with distinct features using multiple-fold references. The forward–backward algorithm (FBA) generates indel-corrected probability information for reliable readout.

The proposed scheme is compatible with next-generation sequencing (NGS) and Oxford Nanopore Technologies (ONT) sequencing platforms (see [Summary of datasets](docs/Summary%20of%20datasets.pdf)). We provide the complete source code and datasets used to generate the recovery results presented in this study
(see [Summary of experiments](docs/Summary%20of%20experiments.pdf)).

---

## Requirements

### **OS Requirements**

The program has been tested on the following operating systems:

- **Ubuntu 18.04.6 LTS**
- **Ubuntu 20.04.6 LTS**

### **Software Requirements**

The following tools and dependencies are required:

- **C compiler** — ensure `gcc` is installed
- **C++ compiler** — ensure `g++` is installed
- **Edlib** — sequence alignment library
  [https://github.com/Martinsos/edlib](https://github.com/Martinsos/edlib)

### Docker Environment

A pre-built Docker image **`bootstrap_readout_v1.0`** is provided to encapsulate the complete experimental software environment.

#### **Docker Environment Summary**

| Environment / Software | Version            |
| ---------------------- | ------------------ |
| Operating system       | Ubuntu 18.04.6 LTS |
| gcc                    | 7.5.0              |
| g++                    | 7.5.0              |
| Edlib                  | 1.2.6              |
| Velvet                 | 1.2.09             |
| ART                    | 2.5.8              |

To start the environment, run:

```bash
cd Docker_image
./run_docker.sh
```

---

## Example of usage

### 1. Fast recovery

Fast recovery was performed on Illumina sequencing data (`DNA-40.5Kb-EM-SE150.fastq`), corresponding to the results in Figure 3B.

**Command:**

```bash
cd ./Recovery_code/Figure3/R0.67/
./build.sh
./recover.sh
```

#### [Step 1] Sliding correlation

Input files:

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `DNA-40.5Kb-EM-SE150.fastq` – real sequencing data with a raw error rate of 0.2%

Output files:

- `correlation_result.txt` – read alignment information including read sequence, correlation peak, alignment position relative to the watermark, and strand orientation

#### [Step 2] Bit-wise majority voting

Input files:

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `correlation_result.txt` – read alignment information from Step 1

Output files:

- `soft_info.txt` – consensus soft information with two columns: probability of "0" and "1" at each bit position

#### [Step 3] Soft-decision LDPC decoding

Input files:

- `soft_info.txt` – consensus soft information from Step 2

Output files:

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary stream (43,200 bits)

Fast recovery workflows for both Illumina (R = 1/4, 1/2, 5/6) and ONT (R = 1/4, 2/3) data are provided, following the same structure and usage as in this example.

---

### 2. Bootstrap recovery

Bootstrap recovery was performed on ONT sequencing data (`DNA-40.5Kb-EM-ONT-1.fastq`), corresponding to the results in Figure 5E.

**Command:**

```bash
cd ./Figure5/bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/
./build.sh
./recover.sh
```

### Stage 1. Recovery with Type-I Reads

#### [Step 0] Read segmentation (only used for ONT sequencing data)

Input files:

- `DNA-40.5Kb-EM-ONT-1.fastq` – ONT sequencing data with a raw error rate of 4.7%
- `Target_length` – target fragment length (typically set to 150 in practice)

Output files:

- `DNA-40.5Kb-EM-ONT-1-segment.fastq` – segmented ONT reads

#### [Step 1] Sliding correlation

Input files:

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `DNA-40.5Kb-EM-ONT-1-segment.fastq`  – segmented ONT reads from Step 0

Output files:

- `correlation_result.txt` – read alignment information including read sequence, correlation peak, alignment position relative to the watermark, and strand orientation
- `TypeI_reads.txt` – high-correlation reads with no indels or indels near the ends
- `lowthres_reads.txt` – low-correlation reads that typically carried internal indel errors

#### [Step 2] Forward-backward algorithm (FBA)

Input files:

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `TypeI_reads.txt` – high-correlation reads from Step 1

Output files:

- `symbol_probability.txt` – indel-corrected symbol probability (from Type-I reads)

#### [Step 3] Consensus soft information generation

Input files:

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `symbol_probability.txt` – indel-corrected symbol probability from Step 2

Output files:

- `soft_info.txt` – consensus soft information (from Type-I reads)

#### [Step 4] Soft-decision LDPC decoding

Input files:

- `soft_info.txt` – consensus soft information from Step 3

Output files:

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary stream (43,200 bits)
- `decodedCodeword.txt` – decoded full LDPC codeword including systematic and parity bits (64,800 bits)

---

### Stage 2. Recovery with Type-I+II Reads

#### [Step 1–2] Scaffold reference generation and filtering Type-II reads

Input files:

- `TypeI_reads.txt` – high-correlation reads from Stage 1
- `lowthres_reads.txt` – low-correlation reads remaining from Stage 1

Output files:

- `TypeII_reads.txt` – Type-II reads that typically carried internal indel errors
- `scaffold_unaligned_reads.txt` – residual reads not aligned to the scaffold

#### [Step 3–5] FBA, consensus, and decoding

Input files:

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `TypeII_reads.txt` – Type-II reads from Step 2

Output files:

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary stream (43,200 bits)
- `decodedCodeword.txt` – decoded full LDPC codeword including systematic and parity bits (64,800 bits)

---

### Stage 3. Recovery with Type-I+II+III Reads

#### [Step 1–2] Regenerative reference generation and filtering Type-III reads

Input files:

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `decodedCodeword.txt` – decoded codeword from Stage 2
- `scaffold_unaligned_reads.txt` – residual reads not aligned to the scaffold from Stage 2

Output files:

- `TypeIII_reads.txt` – reads typically aligned to gap regions, identified during the final recovery stage

#### [Step 3–5] FBA, consensus, and decoding

Input files:

- `SequenceLengthALL_FILE001R0667` – known watermark sequence
- `TypeIII_reads.txt` – Type-III reads from Step 1

Output files:

- `recovery_image.jpg` – reconstructed image
- `recovery_bitstream.txt` – decoded binary stream (43,200 bits)
- `decodedCodeword.txt` – decoded full LDPC codeword including systematic and parity bits (64,800 bits)

Summary statistics for all independent experiments are saved in:

`./bootstrap_recovery_TypeI+II+IIIReads_FBA/R0.67/results/recovery_status.txt`

This file contains seven columns:

1. Experiment ID
2. Erasure rate
3. Substitution error rate
4. Recovery stage (0: not recovered, 1: Type-I reads only, 2: Type-I + II reads, 3: Type-I + II + III reads)
5. Number of Type-I reads
6. Number of Type-II reads
7. Number of Type-III reads

Bootstrap recovery workflow for both Illumina (R = 1/4, 1/2, 2/3, 5/6) and ONT (R = 1/4) data are provided, following the same structure and usage as in this example.

---

## Note

LDPC encoder/decoder by Radford M. Neal [https://github.com/radfordneal/LDPC-codes](https://github.com/radfordneal/LDPC-codes).

For comprehensive information on the coding scheme in this study, please refer to the following literature:

- MacKay, D.J.C., and Neal, R.M. (1996) Near Shannon limit performance of low density parity check codes. Electronics Letters, 32(18): 1645–1646.
- Davey, M.C., and MacKay, D.J.C. (2001) Reliable communication over channels with insertions, deletions, and substitutions. IEEE Transactions on Information Theory, 47(2): 687–698.
- Chen, W.G., Han, M.Z., Zhou, J.T., Ge, Q., Wang, P.P., Zhang, X.C., Zhu, S.Y., Song, L.F., and Yuan, Y.J. (2021) An artificial chromosome for data storage. National Science Review, 8: nwab028.
- Ge, Q., Qin, R., Liu, S., Guo, Q., Han, C., and Chen, W.G. (2025) Pragmatic soft-decision data readout of encoded large DNA. Briefings in Bioinformatics, 26: bbaf102.

## License

This project is licensed under the MIT License. See the [LICENSE] file for details.
