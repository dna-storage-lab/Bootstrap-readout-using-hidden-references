#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#define PATH_MAX 4096 
#define MAXLEN 100000
#define MAX_LINE_LENGTH 1000

int num_threads = 0;  
int watermark_sequence[2][40500] = {0}; 
int seq[200][254886] = {0};  
int segment_length = 0; 
int peak[3000000][300] = {0};  
int read_length = 150;

struct sig {
    int num_threads;      
    int cover;  
    const char *reads_path;
};

void invert_sequence(char *x) {
    int n = strlen(x);
    int m = (n - 1) / 2;
    for (int i = 0; i <= m; i++) {
        char tmp = x[i];
        x[i] = x[n - 1 - i];
        x[n - 1 - i] = tmp;
    }
}

void complement(const char *input, char *output) {
    int n = strlen(input);
    for (int i = 0; i < n; i++) {
        switch (input[i]) {
            case 'A': output[i] = 'T'; break;
            case 'T': output[i] = 'A'; break;
            case 'G': output[i] = 'C'; break;
            case 'C': output[i] = 'G'; break;
            default : output[i] = input[i]; break;
        }
    }
    output[n] = '\0';
}

void load_watermark_reference(const char *path_watermark) {
    FILE *fwm = fopen(path_watermark, "rb");
    if (!fwm) {
        fprintf(stderr, "Can't open watermark file: %s\n", path_watermark);
        exit(EXIT_FAILURE);
    }
    fseek(fwm, 0, SEEK_END);
    long bytes = ftell(fwm);
    rewind(fwm);
    int n_ints = bytes / sizeof(int);
    int *buf = malloc(n_ints * sizeof(int));
    fread(buf, sizeof(int), n_ints, fwm);
    fclose(fwm);
    int half = n_ints / 2;

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < half && j < 40500; j++) {
            int v = buf[i * half + j];
            watermark_sequence[i][j] = (v == 0 ? -1 : v);
        }
    }
    free(buf);
}

void split_reference_segments(int n) {
    segment_length = (40500 + (n - 1) * (read_length+1)) / n;
    for (int i = 0; i < n; i++) {
        int offset = i * (segment_length - read_length - 1);
        for (int j = 0; j < segment_length; j++) {
            seq[2 * i][j]     = watermark_sequence[0][offset + j];
            seq[2 * i + 1][j] = watermark_sequence[1][offset + j];
        }
    }
}

void select_optimal_alignments(int threads, int cover, const char *outpath, const char *reads_path) {
    FILE *fout = fopen(outpath, "w");
    if (!fout) { perror("fopen"); return; }
    FILE *f = fopen(reads_path, "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n", reads_path); fclose(fout); return; }
    char seqstr[502];

    for (int i = 0; i < cover; i++) {
        fscanf(f, "%s", seqstr);
        int best_v = peak[i][0], best_p = peak[i][1], best_z = peak[i][2];
        for (int th = 1; th < threads; th++) {
            int v=peak[i][3*th], p=peak[i][3*th+1], z=peak[i][3*th+2];
            if (v > best_v) { best_v=v; best_p=p; best_z=z; }
        }

        // if (best_z) { invert_sequence(seqstr); complement(seqstr, seqstr); }
        // fprintf(fout, "%s %d %d %d\n", seqstr, best_v, best_p, best_z);

        
        if (best_z) 
        { 
            invert_sequence(seqstr); complement(seqstr, seqstr); 
        }



        size_t len = strlen(seqstr);
        if (len <= 40) {        
            continue;
        } else if (len <= 100) { 
            
            int peak_value = 0;
            fprintf(fout, "%s %d %d %d\n", seqstr, peak_value, best_p, best_z);
        } else {                 
            fprintf(fout, "%s %d %d %d\n", seqstr, best_v, best_p, best_z);
        }


    }
    fclose(f); fclose(fout);
}

void *thread_func(void *arg) {
    struct sig *s = (struct sig *)arg;
    int tid = s->num_threads, cover = s->cover;
    const char *reads_path = s->reads_path; 
    FILE *f = fopen(reads_path, "r");
    if (!f) { 
        fprintf(stderr, "Cannot open %s\n", reads_path); 
        return NULL; 
    }

    char seqstr[502], rcstr[502];
    int reads[2][502], reads_rc[2][502];
    int *corr = malloc(segment_length * sizeof(int)), *corr_rc = malloc(segment_length * sizeof(int));
    for (int n = 0; n < cover; n++) {
        fscanf(f, "%s", seqstr);
        int seqlen = strlen(seqstr);
        strcpy(rcstr, seqstr); invert_sequence(rcstr); complement(rcstr, rcstr);

        for (int i = 0; i < seqlen; i++) {
            char c = seqstr[i]; reads[0][i] = (c=='G'||c=='C')?1:-1; reads[1][i] = (c=='T'||c=='C')?1:-1;
            c = rcstr[i]; reads_rc[0][i] = (c=='G'||c=='C')?1:-1; reads_rc[1][i] = (c=='T'||c=='C')?1:-1;
        }
        for (int i = 0; i < segment_length; i++) corr[i] = corr_rc[i] = 0;

        for (int i = 0; i < segment_length; i++) for (int j = 0; j < seqlen; j++) {
            corr[i] += reads[0][j]*seq[2*tid][i+j] + reads[1][j]*seq[2*tid+1][i+j];
            corr_rc[i] += reads_rc[0][j]*seq[2*tid][i+j] + reads_rc[1][j]*seq[2*tid+1][i+j];
        }
        int pv1 = corr[0], pp1 = 0, pv2 = corr_rc[0], pp2 = 0;

        for (int i = 1; i < segment_length; i++) {
            if (corr[i] > pv1) { pv1 = corr[i]; pp1 = i; }
            if (corr_rc[i] > pv2) { pv2 = corr_rc[i]; pp2 = i; }
        }

        if (pv1 >= pv2) { peak[n][3*tid]=pv1; peak[n][3*tid+1]=pp1+tid*(segment_length-read_length-1); peak[n][3*tid+2]=0; }
        else             { peak[n][3*tid]=pv2; peak[n][3*tid+1]=pp2+tid*(segment_length-read_length-1); peak[n][3*tid+2]=1; }
    }
    free(corr); free(corr_rc); fclose(f);
    return NULL;
}

int main(int argc, char *argv[]) {
    if (argc < 7) {
        fprintf(stderr, "Usage: %s <threads> <length> <coverage>  <watermark_path> <reads_path> <output_path>\n", argv[0]);
        return EXIT_FAILURE;
    }

    num_threads = atoi(argv[1]);
    read_length = atoi(argv[2]);
    int coverage = atoi(argv[3]);
    const char *wmk = argv[4];
    const char *reads_path = argv[5]; 
    const char *outpath = argv[6];

    load_watermark_reference(wmk);
    split_reference_segments(num_threads);

    pthread_t threads_id[num_threads];
    struct sig args[num_threads];

    for (int i = 0; i < num_threads; i++) { 
        args[i].num_threads = i; 
        args[i].cover = coverage; 
        args[i].reads_path = reads_path;
    }
    for (int i = 0; i < num_threads; i++) pthread_create(&threads_id[i], NULL, thread_func, &args[i]);
    for (int i = 0; i < num_threads; i++) pthread_join(threads_id[i], NULL);

    select_optimal_alignments(num_threads, coverage, outpath, reads_path);
    return 0;
}
