#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

void split_by_threshold(const char *input_file, double threshold_ratio,
                        const char *output_file_up, const char *output_file_low) {
    int peak = 0, pos = 0, zf1 = 0;
    char buff[500] = {0};

    FILE *fpeak = fopen(input_file, "rt");
    if (fpeak == NULL) {
        printf("Cannot open the input file: %s\n", input_file);
        return;
    }

    FILE *farti_up = fopen(output_file_up, "wt");
    if (farti_up == NULL) {
        printf("Cannot open the output file: %s\n", output_file_up);
        fclose(fpeak);
        return;
    }

    FILE *farti_low = fopen(output_file_low, "wt");
    if (farti_low == NULL) {
        printf("Cannot open the output file: %s\n", output_file_low);
        fclose(fpeak);
        fclose(farti_up);
        return;
    }

    while (1) {
        memset(buff, 0, sizeof(buff));
        int ret = fscanf(fpeak, "%499s %d %d %d", buff, &peak, &pos, &zf1);
        if (ret == EOF) break;

        if (ret != 4) {
            int c;
            while ((c = fgetc(fpeak)) != '\n' && c != EOF) {}
            continue;
        }

        int read_length = (int)strlen(buff);
        double dynamic_threshold = threshold_ratio * (double)read_length * 2.0;

        if ((double)peak >= dynamic_threshold) {
            fprintf(farti_up, "%s %d %d %d\n", buff, peak, pos, zf1);
        } else {
            fprintf(farti_low, "%s %d %d %d\n", buff, peak, pos, zf1);
        }
    }

    fclose(fpeak);
    fclose(farti_up);
    fclose(farti_low);
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Usage: ./split_and_threshold <sliding_result> <cover> <threshold_ratio> <output_file_up> <output_file_low>\n");
        printf("Example: ./split_and_threshold result.txt 10 0.32 up.txt low.txt\n");
        return 1;
    }

    const char *sliding_result = argv[1];
    (void)atoi(argv[2]); 
    double threshold_ratio = atof(argv[3]);
    const char *output_file_up = argv[4];
    const char *output_file_low = argv[5];

    if (threshold_ratio < 0) {
        printf("threshold_ratio must be non-negative (e.g., 0.32).\n");
        return 1;
    }

    split_by_threshold(sliding_result, threshold_ratio, output_file_up, output_file_low);
    return 0;
}
