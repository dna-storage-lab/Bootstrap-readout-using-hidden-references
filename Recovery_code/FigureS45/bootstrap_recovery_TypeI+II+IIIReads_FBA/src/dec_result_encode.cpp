#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define COLUMN 180*360  // 64800

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc != 3) {
        printf("Usage: %s dec_results_64800.txt output_genome.txt\n", argv[0]);
        return 1;
    }

    FILE *fp1;

    // 内存分配
    char temp2[COLUMN];
    char *cblk = temp2;

    char temp3[COLUMN];
    char *cblk1 = temp3;

    char *sparse_code = (char *)malloc(81000 * sizeof(char));
    if (sparse_code == NULL) {
        printf("Error: malloc failed\n");
        return 1;
    }

    // 1. 读取输入比特串文件
    FILE *input_file = fopen(argv[1], "r");
    if (input_file == NULL) {
        printf("Can't open input file: %s\n", argv[1]);
        return 1;
    }

    int ch, count = 0;
    while ((ch = fgetc(input_file)) != EOF && count < COLUMN) {
        if (ch == '0' || ch == '1') {
            cblk1[count] = ch - '0';
            count++;
        }
    }
    fclose(input_file);

    if (count != COLUMN) {
        printf("Error: Expected 64800 bits, but got %d\n", count);
        return 1;
    }

    // 2. 加载 permutation 参数
    char para_file[] = "./configure/permutation64800";
    if ((fp1 = fopen(para_file, "rb")) == NULL) {
        printf("Can't open parameter file: %s\n", para_file);
        return 1;
    }

    int *permutation64800 = new int[COLUMN];
    fread(permutation64800, sizeof(int), COLUMN, fp1);
    fclose(fp1);

    // 3. 交织
    for (int j = 0; j < COLUMN; j++) {
        // cblk[j] = cblk1[permutation64800[j]];
        cblk[j] = *(cblk1 + permutation64800[j]);
        // printf("%d ", permutation64800[j]);		
    }

    // 4. 稀疏化（64800 -> 81000）
    int temp4, tk;
    for (int j = 0; j < 64800 / 4; j++) {
        temp4 = 0;
        for (tk = 0; tk < 4; tk++) {
            temp4 += cblk[tk + j * 4];
        }

        if (temp4 < 3) {
            sparse_code[j * 5 + 0] = 0;
            for (tk = 1; tk < 5; tk++) {
                sparse_code[j * 5 + tk] = cblk[tk - 1 + j * 4];
            }
        } else {
            sparse_code[j * 5 + 0] = 1;
            for (tk = 1; tk < 5; tk++) {
                sparse_code[j * 5 + tk] = 1 - cblk[tk - 1 + j * 4];
            }
        }
    }

    // 5. 输出稀疏码
    // FILE *output_file = fopen(argv[2], "w");
    // if (output_file == NULL) {
    //     printf("Can't open output file: %s\n", argv[2]);
    //     return 1;
    // }

    // for (int i = 0; i < 81000; i++) {
    //     fputc(sparse_code[i] + '0', output_file);
    // }
    // fclose(output_file);

    // 6. 加载水印并异或
    int *watermark = (int*)malloc(81000 * sizeof(int));
    FILE *watermarkFP;
    if ((watermarkFP = fopen("./configure/SequenceL81000NoPeriodOnly2ndFILE", "rb")) == NULL) {
        printf("Can't open the watermark file\n");
        return 1;
    }
    fread(watermark, sizeof(int), 81000, watermarkFP);
    fclose(watermarkFP);

    // 7. 异或操作并生成碱基序列
    FILE *genome_file = fopen(argv[2], "w");
    if (genome_file == NULL) {
        printf("Can't open genome output file: %s\n", argv[2]);
        return 1;
    }

    for (int j = 0; j < 40500; j++) {
        int bit1 = (sparse_code[j] ^ watermark[j]);
        int bit2 = (sparse_code[j + 40500] ^ watermark[j + 40500]);

        char base;
        if (bit1 == 0 && bit2 == 0) base = 'A';
        else if (bit1 == 0 && bit2 == 1) base = 'T';
        else if (bit1 == 1 && bit2 == 0) base = 'G';
        else base = 'C';

        fputc(base, genome_file);
    }
    fclose(genome_file);

    // 清理
    delete[] permutation64800;
    free(sparse_code);
    free(watermark);

    // printf("Success: base sequence written to %s\n", argv[2]);
    return 0;
}