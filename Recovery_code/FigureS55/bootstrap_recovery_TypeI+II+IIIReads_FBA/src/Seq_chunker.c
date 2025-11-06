#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define BUFFER_SIZE 1048576

int read_fastq_record(FILE *in, char *header, char *seq, char *plus, char *qual) {
    if (fgets(header, BUFFER_SIZE, in) == NULL) return 0;
    if (header[0] != '@') return 0;
    
    if (fgets(seq, BUFFER_SIZE, in) == NULL) return 0;
    
    if (fgets(plus, BUFFER_SIZE, in) == NULL) return 0;
    if (plus[0] != '+') return 0;
    
    if (fgets(qual, BUFFER_SIZE, in) == NULL) return 0;
    
    header[strcspn(header, "\r\n")] = '\0';
    seq[strcspn(seq, "\r\n")] = '\0';
    plus[strcspn(plus, "\r\n")] = '\0';
    qual[strcspn(qual, "\r\n")] = '\0';
    
    return 1;
}

void write_fastq_record(FILE *out, const char *header, const char *seq, 
                       const char *plus, const char *qual) {
    fprintf(out, "%s\n", header);
    fprintf(out, "%s\n", seq);
    fprintf(out, "%s\n", plus);
    fprintf(out, "%s\n", qual);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input FASTQ file> <output FASTQ file> <split length L>\n", argv[0]);
        return 1;
    }
    
    char *endptr;
    long L = strtol(argv[3], &endptr, 10);
    if (*endptr != '\0' || L <= 0 || L > INT_MAX) {
        fprintf(stderr, "Error: Split length must be a positive integer\n");
        return 1;
    }
    
    FILE *in = fopen(argv[1], "r");
    if (!in) {
        perror("Failed to open input file");
        return 1;
    }
    
    FILE *out = fopen(argv[2], "w");
    if (!out) {
        perror("Failed to open output file");
        fclose(in);
        return 1;
    }
    
    char *header = malloc(BUFFER_SIZE);
    char *seq = malloc(BUFFER_SIZE);
    char *plus = malloc(BUFFER_SIZE);
    char *qual = malloc(BUFFER_SIZE);
    char *new_header = malloc(BUFFER_SIZE);
    char *sub_seq = malloc(BUFFER_SIZE);
    char *sub_qual = malloc(BUFFER_SIZE);
    
    if (!header || !seq || !plus || !qual || !new_header || !sub_seq || !sub_qual) {
        fprintf(stderr, "Memory allocation failed\n");
        free(header);
        free(seq);
        free(plus);
        free(qual);
        free(new_header);
        free(sub_seq);
        free(sub_qual);
        fclose(in);
        fclose(out);
        return 1;
    }
    
    int record_count = 0;
    while (read_fastq_record(in, header, seq, plus, qual)) {
        record_count++;
        int seq_len = strlen(seq);
        
        if (seq_len != strlen(qual)) {
            fprintf(stderr, "Warning: Record %d has inconsistent sequence and quality lengths, skipped\n", record_count);
            continue;
        }
        
        if (seq_len < L) {
            write_fastq_record(out, header, seq, plus, qual);
        } else {
            int n_total = (seq_len + L - 1) / L;
            
            for (int i = 0; i < n_total; i++) {
                int start = i * L;
                int end = start + L;
                if (end > seq_len) end = seq_len;
                int chunk_len = end - start;
                
                strncpy(sub_seq, seq + start, chunk_len);
                sub_seq[chunk_len] = '\0';
                strncpy(sub_qual, qual + start, chunk_len);
                sub_qual[chunk_len] = '\0';
                
                snprintf(new_header, BUFFER_SIZE, "%s_chunk%d/%d", header, i + 1, n_total);
                
                write_fastq_record(out, new_header, sub_seq, plus, sub_qual);
            }
        }
    }
    
    free(header);
    free(seq);
    free(plus);
    free(qual);
    free(new_header);
    free(sub_seq);
    free(sub_qual);
    fclose(in);
    fclose(out);
    
    return 0;
}
    