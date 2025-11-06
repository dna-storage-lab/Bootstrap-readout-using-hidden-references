#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

int MAX_LEN = 40320; 

void major1(const char *input_dir, const char *output_dir) {
    FILE *fpeak = NULL, *fcom = NULL;

    char peakname[256], comname[256];
    snprintf(peakname, sizeof(peakname), "%s/Type-I_reads.txt", input_dir);
    snprintf(comname, sizeof(comname), "%s/scaffold_ref.txt", output_dir);

    fpeak = fopen(peakname, "rt");
    if (fpeak == NULL) {
        printf("Error: Cannot open input peak file %s\n", peakname);
        return;
    }

    fcom = fopen(comname, "wt");
    if (!fcom) {
        printf("Error: Cannot create output file %s\n", comname);
        fclose(fpeak);
        return;
    }

    char buff[500] = {0};
    char *com = (char *)malloc(sizeof(char) * MAX_LEN);
    int (*reads)[5] = (int (*)[5])calloc(MAX_LEN, sizeof(int[5]));

    if (!com || !reads) {
        printf("Memory allocation failed.\n");
        fclose(fpeak); fclose(fcom);
        free(com); free(reads);
        return;
    }

    int peak, pos, zf;

    while (fscanf(fpeak, "%s %d %d %d", 
                  buff, &peak, &pos, &zf) != EOF) {
        int seqlen = strlen(buff);
        for (int i = 0; i < seqlen && (pos + i < MAX_LEN); i++) {
            switch (buff[i]) {
                case 'A': reads[pos + i][0]++; break;
                case 'T': reads[pos + i][1]++; break;
                case 'G': reads[pos + i][2]++; break;
                case 'C': reads[pos + i][3]++; break;
                default:  reads[pos + i][4]++; break;
            }
        }
    }
    fclose(fpeak);

    for (int i = 0; i < MAX_LEN; i++) {
        int max_count = 0, best_index = 4;
        for (int j = 0; j < 4; j++) {
            if (reads[i][j] > max_count) {
                max_count = reads[i][j];
                best_index = j;
            }
        }
        com[i] = (best_index == 0) ? 'A' : (best_index == 1) ? 'T' : (best_index == 2) ? 'G' : (best_index == 3) ? 'C' : 'X';
        fprintf(fcom, "%c", com[i]);
    }
    fprintf(fcom, "\n");

    fclose(fcom);
     free(com); free(reads);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s <input_dir> <output_dir> \n", argv[0]);
        return 1;
    }

    const char *input_dir = argv[1];
    const char *output_dir = argv[2];

    MAX_LEN = 40500;
    major1(input_dir, output_dir);


    return 0;
}
