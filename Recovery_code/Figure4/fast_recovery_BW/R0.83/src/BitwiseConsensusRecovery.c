#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define MAXLEN 3000000
#define MAX_READS 254886
#define MAX_LINE_LENGTH 500
#define SEQ_LEN 40500
#define EPSILON 1e-6

typedef struct {
    char vote;
    double prob;
} PositionResult;

typedef struct {
    double p0;
    double p1;
} ProbPair;


void sparsify_group(PositionResult group[5], ProbPair output[4]) {
    const char combinations[16][5] = {
        {'0','0','0','0','0'}, {'0','0','0','0','1'}, {'0','0','0','1','0'}, {'0','0','0','1','1'},
        {'0','0','1','0','0'}, {'0','0','1','0','1'}, {'0','0','1','1','0'}, {'1','1','0','0','0'},
        {'0','1','0','0','0'}, {'0','1','0','0','1'}, {'0','1','0','1','0'}, {'1','0','1','0','0'},
        {'0','1','1','0','0'}, {'1','0','0','1','0'}, {'1','0','0','0','1'}, {'1','0','0','0','0'}
    };

    double p[16] = {0};
    for (int i = 0; i < 16; i++) {
        double prob = 1.0;
        for (int j = 0; j < 5; j++) {
            double adjusted_prob = group[j].prob;
            if (adjusted_prob == 0.000000) adjusted_prob = 0.000100;
            else if (adjusted_prob == 1.000000) adjusted_prob = 0.999900;

            if (combinations[i][j] == '0') {
                if (group[j].vote == '0') prob *= adjusted_prob;
                else if (group[j].vote == '1') prob *= (1 - adjusted_prob);
                else prob *= 0.5; 
            } else if (combinations[i][j] == '1') {
                if (group[j].vote == '0') prob *= (1 - adjusted_prob);
                else if (group[j].vote == '1') prob *= adjusted_prob;
                else prob *= 0.5; 
            }
        }
        p[i] = prob;
    }

    double sum_p = 0.0;
    for (int i = 0; i < 16; i++) sum_p += p[i];
    for (int i = 0; i < 16; i++) p[i] /= sum_p;

    output[0].p0 = p[0]+p[1]+p[2]+p[3]+p[4]+p[5]+p[6]+p[7];
    output[0].p1 = 1 - output[0].p0;
    
    output[1].p0 = p[0]+p[1]+p[2]+p[3]+p[8]+p[9]+p[10]+p[11];
    output[1].p1 = 1 - output[1].p0;
    
    output[2].p0 = p[0]+p[1]+p[4]+p[5]+p[8]+p[9]+p[12]+p[13];
    output[2].p1 = 1 - output[2].p0;
    
    output[3].p0 = p[0]+p[2]+p[4]+p[6]+p[8]+p[10]+p[12]+p[14];
    output[3].p1 = 1 - output[3].p0;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <filtered_in> <peak_threshold> <wmkfile> <output_softfile>\n", argv[0]);
        return 1;
    }
    
    const char *filtered_in = argv[1];
    int peak_threshold = atoi(argv[2]);
    const char *wmkfile = argv[3];
    const char *softfile = argv[4];

    FILE *fpeak = fopen(filtered_in, "rt");
    if (!fpeak) {
        fprintf(stderr, "Cannot open %s\n", filtered_in);
        return 1;
    }

    char seq[MAX_LINE_LENGTH];
    int peak_value, pos, zf;
    int (*counts)[4] = calloc(SEQ_LEN, sizeof(int[4])); 

    while (fscanf(fpeak, "%s %d %d %d", seq, &peak_value, &pos, &zf) == 4) {
        if (peak_value < peak_threshold) continue;
        
        int len = strlen(seq);
        for (int i = 0; i < len; i++) {
            int idx = (pos + i) ;
            switch (seq[i]) {
                case 'A': counts[idx][0]++; break;
                case 'T': counts[idx][1]++; break;
                case 'G': counts[idx][2]++; break;
                case 'C': counts[idx][3]++; break;
            }
        }
    }
    fclose(fpeak);

    PositionResult *upper_results = malloc(SEQ_LEN * sizeof(PositionResult));
    PositionResult *lower_results = malloc(SEQ_LEN * sizeof(PositionResult));
    
    for (int i = 0; i < SEQ_LEN; i++) {
        int up0 = counts[i][0] + counts[i][1]; 
        int up1 = counts[i][2] + counts[i][3]; 
        int up_total = up0 + up1;
        
        upper_results[i].vote = (up0 > up1) ? '0' : (up0 < up1) ? '1' : 'X';
        upper_results[i].prob = (up_total == 0) ? 0.5 : 
                              (up0 > up1) ? (double)up0 / up_total : (double)up1 / up_total;

        int lo0 = counts[i][0] + counts[i][2]; 
        int lo1 = counts[i][1] + counts[i][3]; 
        int lo_total = lo0 + lo1;
        
        lower_results[i].vote = (lo0 > lo1) ? '0' : (lo0 < lo1) ? '1' : 'X';
        lower_results[i].prob = (lo_total == 0) ? 0.5 : 
                              (lo0 > lo1) ? (double)lo0 / lo_total : (double)lo1 / lo_total;
    }
    free(counts);

    int *watermark = malloc(2 * SEQ_LEN * sizeof(int));
    FILE *fw = fopen(wmkfile, "rb");
    if (!fw || fread(watermark, sizeof(int), 2*SEQ_LEN, fw) != 2*SEQ_LEN) {
        fprintf(stderr, "Error reading watermark\n");
        free(upper_results); free(lower_results); free(watermark);
        if (fw) fclose(fw);
        return 1;
    }
    fclose(fw);
    
    for (int i = 0; i < SEQ_LEN; i++) {
        if (watermark[i]) {
            if (upper_results[i].vote == '0') upper_results[i].vote = '1';
            else if (upper_results[i].vote == '1') upper_results[i].vote = '0';
        }
        if (watermark[SEQ_LEN + i]) { 
            if (lower_results[i].vote == '0') lower_results[i].vote = '1';
            else if (lower_results[i].vote == '1') lower_results[i].vote = '0';
        }
    }
    free(watermark);

    int num_groups = SEQ_LEN / 5;
    ProbPair *sparse_upper = malloc(num_groups * 4 * sizeof(ProbPair));
    ProbPair *sparse_lower = malloc(num_groups * 4 * sizeof(ProbPair));
    
    for (int g = 0; g < num_groups; g++) {
        PositionResult upper_group[5], lower_group[5];
        for (int i = 0; i < 5; i++) {
            int idx = g * 5 + i;
            upper_group[i] = upper_results[idx];
            lower_group[i] = lower_results[idx];
        }
        
        sparsify_group(upper_group, &sparse_upper[g*4]);
        sparsify_group(lower_group, &sparse_lower[g*4]);
    }
    
    free(upper_results);
    free(lower_results);

    FILE *fo = fopen(softfile, "wt");
    if (!fo) {
        fprintf(stderr, "Cannot create %s\n", softfile);
        free(sparse_upper); free(sparse_lower);
        return 1;
    }
    
    for (int i = 0; i < num_groups * 4; i++) {
        double p0 = sparse_upper[i].p0;
        double p1 = sparse_upper[i].p1;
        p0 = (p0 < EPSILON) ? 0.0001 : ((1.0 - p0) < EPSILON ? 0.9999 : p0);
        p1 = (p1 < EPSILON) ? 0.0001 : ((1.0 - p1) < EPSILON ? 0.9999 : p1);
        fprintf(fo, "%.6f %.6f\n", p0, p1);
    }
    
    for (int i = 0; i < num_groups * 4; i++) {
        double p0 = sparse_lower[i].p0;
        double p1 = sparse_lower[i].p1;
        p0 = (p0 < EPSILON) ? 0.0001 : ((1.0 - p0) < EPSILON ? 0.9999 : p0);
        p1 = (p1 < EPSILON) ? 0.0001 : ((1.0 - p1) < EPSILON ? 0.9999 : p1);
        fprintf(fo, "%.6f %.6f\n", p0, p1);
    }

    fclose(fo);
    
    free(sparse_upper);
    free(sparse_lower);
    return 0;
}