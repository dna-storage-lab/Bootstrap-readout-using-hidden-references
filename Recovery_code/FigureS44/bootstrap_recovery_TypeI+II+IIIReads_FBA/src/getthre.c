#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

void split_by_threshold(const char *input_file, int threshold, char *output_file_up, char *output_file_low) {
    int peak = 0, pos = 0, zf1 = 0;
    char buff[500] = {0};

    FILE *fpeak = NULL;
    FILE *farti_up = NULL;
    FILE *farti_low = NULL;

    fpeak = fopen(input_file, "rt");
    if (fpeak == NULL) {
        printf("Cannot open the input file: %s\n", input_file);
        return;
    }
    
    farti_up = fopen(output_file_up, "wt");
    if (farti_up == NULL) {
        printf("Cannot open the output file: %s\n", output_file_up);
        fclose(fpeak);
        return;
    }

    farti_low = fopen(output_file_low, "wt");
    if (farti_low == NULL) {
        printf("Cannot open the output file: %s\n", output_file_low);
        fclose(fpeak);
        fclose(farti_up);
        return;
    }

    while (1) {
        memset(buff, '\0', sizeof(buff));
        int ret = fscanf(fpeak, "%s %d %d %d",
                         buff, &peak, &pos, &zf1);
        if (ret == EOF) {
            break;
        }

        if (peak >= threshold) {
            fprintf(farti_up, "%s %d %d %d\n",
                    buff, peak, pos, zf1);
        } else {
            fprintf(farti_low, "%s %d %d %d\n",
                    buff, peak, pos, zf1);
        }
    }

    fclose(fpeak);
    fclose(farti_up);
    fclose(farti_low);
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Usage: ./split_and_threshold <sliding_result> <cover> <peak_threshold> <output_file_up> <output_file_low>\n");
        return 1;
    }

    char *sliding_result = argv[1];  
    int cover = atoi(argv[2]); 
    int peak_thres = atoi(argv[3]);
    char *output_file_up = argv[4];  
    char *output_file_low = argv[5]; 

    const char *input_file = sliding_result;
    split_by_threshold(input_file, peak_thres, output_file_up, output_file_low);

    return 0;
}
