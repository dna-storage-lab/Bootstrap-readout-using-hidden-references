#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_LEN 100000 
#define EPSILON 1e-6 

void calculate_and_write_error(const char *input_file, const char *reference_file, int num_lines, FILE *output_file) {
    FILE *input_fp = fopen(input_file, "r");
    FILE *ref_fp = fopen(reference_file, "r");

    if (input_fp == NULL || ref_fp == NULL) {
        fprintf(stderr, "Error opening files\n");
        return;
    }

    double input_prob[num_lines][2]; 
    int reference[num_lines];         

    for (int i = 0; i < num_lines; i++) {
        fscanf(input_fp, "%lf %lf", &input_prob[i][0], &input_prob[i][1]);
    }

    for (int i = 0; i < num_lines; i++) {
        char ch = fgetc(ref_fp); 
        if (ch == '0') {
            reference[i] = 0; 
        } else if (ch == '1') {
            reference[i] = 1; 
        } else if (ch == EOF) {
            if (i < num_lines) {
                fprintf(stderr, "Error: File contains fewer than %d bits.\n", num_lines);
            }
            break;
        } else {
            fprintf(stderr, "Error: Invalid character in reference file at line %d\n", i);
            fclose(input_fp);
            fclose(ref_fp);
            return;
        }
    }


    if (fgetc(ref_fp) != EOF) {
        fprintf(stderr, "Error: File contains more than %d bits.\n", num_lines);
        fclose(input_fp);
        fclose(ref_fp);
        return;
    }


    int erase_count = 0;   
    int replace_count = 0; 
    int predicted[num_lines];

    for (int i = 0; i < num_lines; i++) {
        if (fabs(input_prob[i][0] - input_prob[i][1]) < EPSILON) {
            predicted[i] = 2;  
            erase_count++;  
        } else if (input_prob[i][0] > input_prob[i][1]) {
            predicted[i] = 0; 
        } else {
            predicted[i] = 1;  
        }

        if (predicted[i] != reference[i]) {
            if (predicted[i] != 2) {
                replace_count++;  
            }
        }
    }

    double erase_rate = (double)erase_count / num_lines;
    double replace_rate = (double)replace_count / (num_lines-erase_count);

    double total_error_rate = (erase_rate / 2) + replace_rate;

    fprintf(output_file, "%d %d %.6f %.6f %.6f\n", erase_count, replace_count, erase_rate, replace_rate, total_error_rate);

    fclose(input_fp);
    fclose(ref_fp);
}


int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: ./cal_error <input_file> <decoded_file> <NUM_LINES> <output_file>\n");
        return 1;
    }

    const char *reference_file = argv[2];  
    int num_lines = atoi(argv[3]);
    const char *output_file_name = argv[4];  

    FILE *output_fp = fopen(output_file_name, "a");
    if (output_fp == NULL) {
        fprintf(stderr, "Error opening output file\n");
        return 1;
    }

    calculate_and_write_error(argv[1], reference_file, num_lines, output_fp);

    fclose(output_fp);
    return 0;
}