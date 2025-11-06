#include <stdio.h>
#include <stdlib.h>

void calculate_hamming_distance(const char *file1, const char *file2, const char *file3, int file_size) {
    FILE *f1 = fopen(file1, "r");
    FILE *f2 = fopen(file2, "r");
    FILE *f_out = fopen(file3, "a");

    if (f1 == NULL || f2 == NULL || f_out == NULL) {
        perror("Failed to open input or output files");
        exit(EXIT_FAILURE);
    }

    int hamming_distance = 0;
    char bit1, bit2;
    int count = 0;

    while (count < file_size) {
        if (fscanf(f1, " %c", &bit1) != 1 || fscanf(f2, " %c", &bit2) != 1) {
            fprintf(stderr, "Error reading files or file length mismatch (expected %d bits)\n", file_size);
            fclose(f1);
            fclose(f2);
            fclose(f_out);
            exit(EXIT_FAILURE);
        }

        if ((bit1 != '0' && bit1 != '1') || (bit2 != '0' && bit2 != '1')) {
            fprintf(stderr, "Non-binary character detected in input files\n");
            fclose(f1);
            fclose(f2);
            fclose(f_out);
            exit(EXIT_FAILURE);
        }

        if (bit1 != bit2) {
            hamming_distance++;
        }

        count++;
    }

    double error_rate = (double)hamming_distance / file_size;
    fprintf(f_out, "%d %.6f\n", hamming_distance, error_rate);

    fclose(f1);
    fclose(f2);
    fclose(f_out);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s file1.txt file2.txt output_distance.txt file_size\n", argv[0]);
        return EXIT_FAILURE;
    }

    int file_size = atoi(argv[4]);
    if (file_size <= 0) {
        fprintf(stderr, "Error: file_size must be a positive integer\n");
        return EXIT_FAILURE;
    }

    calculate_hamming_distance(argv[1], argv[2], argv[3], file_size);

    return EXIT_SUCCESS;
}