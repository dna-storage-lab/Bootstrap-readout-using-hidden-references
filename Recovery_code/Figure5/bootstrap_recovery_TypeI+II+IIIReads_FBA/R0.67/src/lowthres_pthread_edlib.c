#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "edlib.h"

#define MAXSEQ 500
#define MAXBUF 40500

typedef struct {
    char  **lines;
    size_t  size;
    size_t  capacity;
} ThreadResult;

typedef struct {
    char *seq;          
    int  seq_len;    
    int  pea, pos, zf;
} PeakRecord;

typedef struct {
    PeakRecord   *records;  
    size_t        start;  
    size_t        end;      
    const char   *ref_seq; 
    int           ref_len;
    ThreadResult *res;    
} ThreadArg;

void inv(char* x, int n) {
    for (int i = 0, j = n - 1; i < j; i++, j--) {
        char t = x[i]; x[i] = x[j]; x[j] = t;
    }
}

void complement(const char *in, char *out, int len) {
    for (int i = 0; i < len; i++) {
        switch (in[i]) {
            case 'A': out[i] = 'T'; break;
            case 'T': out[i] = 'A'; break;
            case 'G': out[i] = 'C'; break;
            case 'C': out[i] = 'G'; break;
            default:  out[i] = 'N'; break;
        }
    }
}

void* thread_func(void *vp) {
    ThreadArg *A = (ThreadArg*)vp;
    ThreadResult *res = A->res;
    char linebuf[512];

    for (size_t i = A->start; i < A->end; i++) {
        PeakRecord *r = &A->records[i];
        int read_length = r->seq_len;
        
        int half_len = read_length / 2;
        int first_half_len = half_len;
        int second_half_len = read_length - half_len;

        char *before_half = malloc(first_half_len + 1);
        char *after_half = malloc(second_half_len + 1);
        char *tmp = malloc(second_half_len + 1); 
        
        if (!before_half || !after_half || !tmp) {
            if (before_half) free(before_half);
            if (after_half) free(after_half);
            if (tmp) free(tmp);
            continue;
        }
        

        strncpy(before_half, r->seq, first_half_len);
        before_half[first_half_len] = '\0';
        
        strncpy(after_half, r->seq + first_half_len, second_half_len);
        after_half[second_half_len] = '\0';

        int bf_ins=0, bf_del=0, bf_sub=0, ed_bf, pos_bf;
        int br_ins=0, br_del=0, br_sub=0, ed_br, pos_br;
        int af_ins=0, af_del=0, af_sub=0, ed_af, pos_af;
        int ar_ins=0, ar_del=0, ar_sub=0, ed_ar, pos_ar;
        EdlibAlignResult R;
        

        R = edlibAlign(before_half, first_half_len, A->ref_seq, A->ref_len,
                       edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        ed_bf = R.editDistance; 
        pos_bf = R.numLocations > 0 ? R.startLocations[0] : -1;

        for (int k = 0; k < R.alignmentLength; k++) {
            if (R.alignment[k] == 1) bf_del++;
            else if (R.alignment[k] == 2) bf_ins++;
            else if (R.alignment[k] == 3) bf_sub++;
        }
        edlibFreeAlignResult(R);

        strcpy(tmp, before_half);
        inv(tmp, first_half_len);
        complement(tmp, tmp, first_half_len);
        R = edlibAlign(tmp, first_half_len, A->ref_seq, A->ref_len,
                       edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        ed_br = R.editDistance; 
        pos_br = R.numLocations > 0 ? R.startLocations[0] : -1;
        for (int k = 0; k < R.alignmentLength; k++) {
            if (R.alignment[k] == 1) br_del++;
            else if (R.alignment[k] == 2) br_ins++;
            else if (R.alignment[k] == 3) br_sub++;
        }
        edlibFreeAlignResult(R);
        
        
        R = edlibAlign(after_half, second_half_len, A->ref_seq, A->ref_len,
                       edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        ed_af = R.editDistance; 
        pos_af = R.numLocations > 0 ? R.startLocations[0] : -1;
        for (int k = 0; k < R.alignmentLength; k++) {
            if (R.alignment[k] == 1) af_del++;
            else if (R.alignment[k] == 2) af_ins++;
            else if (R.alignment[k] == 3) af_sub++;
        }
        edlibFreeAlignResult(R);
        
       
        strcpy(tmp, after_half);
        inv(tmp, second_half_len);
        complement(tmp, tmp, second_half_len);
        R = edlibAlign(tmp, second_half_len, A->ref_seq, A->ref_len,
                       edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        ed_ar = R.editDistance; 
        pos_ar = R.numLocations > 0 ? R.startLocations[0] : -1;
        for (int k = 0; k < R.alignmentLength; k++) {
            if (R.alignment[k] == 1) ar_del++;
            else if (R.alignment[k] == 2) ar_ins++;
            else if (R.alignment[k] == 3) ar_sub++;
        }
        edlibFreeAlignResult(R);
        
       
        int min_ed = ed_bf;
        int min_pos = pos_bf;
        int choice = 0;
        
        if (ed_br < min_ed) { 
            min_ed = ed_br; 
            min_pos = pos_br; 
            choice = 1; 
        }
        if (ed_af < min_ed) { 
            min_ed = ed_af; 
            min_pos = pos_af; 
            choice = 2; 
        }
        if (ed_ar < min_ed) { 
            min_ed = ed_ar; 
            min_pos = pos_ar; 
            choice = 3; 
        }
        
       
        double ins, del, sub, err;
        int orient, outpos;
        int used_length;
        
        switch(choice) {
            case 0: 
                ins = bf_ins / (double)first_half_len;
                del = bf_del / (double)first_half_len;
                sub = bf_sub / (double)first_half_len;
                orient = 0;
                outpos = min_pos;
                used_length = first_half_len;
                break;
            case 1: 
                ins = br_ins / (double)first_half_len;
                del = br_del / (double)first_half_len;
                sub = br_sub / (double)first_half_len;
                orient = 1;
                outpos = min_pos - first_half_len;
                used_length = first_half_len;
                break;
            case 2: 
                ins = af_ins / (double)second_half_len;
                del = af_del / (double)second_half_len;
                sub = af_sub / (double)second_half_len;
                orient = 0;
                outpos = min_pos - first_half_len;
                used_length = second_half_len;
                break;
            default: 
                ins = ar_ins / (double)second_half_len;
                del = ar_del / (double)second_half_len;
                sub = ar_sub / (double)second_half_len;
                orient = 1;
                outpos = min_pos;
                used_length = second_half_len;
                break;
        }
        
        err = min_ed / (double)used_length;
        
        if (err < 0.2 && min_pos >= 0) {
            int L = snprintf(linebuf, sizeof(linebuf),
                           "%s %d %d %d %d %d %.6f %d\n",
                           r->seq, r->pea, r->pos, r->zf,
                           min_ed, outpos, err, orient);
            if (L > 0 && L < (int)sizeof(linebuf)) {
                if (res->size + 1 >= res->capacity) {
                    res->capacity = (res->capacity == 0) ? 1024 : res->capacity * 2;
                    res->lines = realloc(res->lines, res->capacity * sizeof(char*));
                }
                char *dup = malloc(L + 1);
                memcpy(dup, linebuf, L + 1);
                res->lines[res->size++] = dup;
            }
        }
        

        free(before_half);
        free(after_half);
        free(tmp);
    }
    return res;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <threads> <reference_file> <peak_file> <output_file>\n", argv[0]);
        return 1;
    }
    int t = atoi(argv[1]);
    const char *ref_file  = argv[2];
    const char *peak_file = argv[3];
    const char *out_file  = argv[4];


    FILE *fp = fopen(peak_file, "r");
    if (!fp) { 
        perror("open peak"); 
        return 1; 
    }
    

    size_t total = 0;
    char buffer[MAXBUF];
    while (fgets(buffer, sizeof(buffer), fp)) {
        total++;
    }
    rewind(fp);


    PeakRecord *all = malloc(total * sizeof(*all));
    if (!all) {
        perror("malloc all");
        fclose(fp);
        return 1;
    }
    
    for (size_t i = 0; i < total; i++) {
        if (fgets(buffer, sizeof(buffer), fp)) {
            char seq_buf[MAXSEQ];
            int pea, pos, zf;
            if (sscanf(buffer, "%s %d %d %d", seq_buf, &pea, &pos, &zf) == 4) {
                int len = strlen(seq_buf);
                all[i].seq_len = len;
                all[i].pea = pea;
                all[i].pos = pos;
                all[i].zf = zf;
                all[i].seq = malloc(len + 1);
                if (all[i].seq) {
                    strcpy(all[i].seq, seq_buf);
                } else {
                    all[i].seq_len = 0;
                }
            } else {
                all[i].seq_len = 0;
                all[i].seq = NULL;
            }
        } else {
            all[i].seq_len = 0;
            all[i].seq = NULL;
        }
    }
    fclose(fp);


    FILE *fr = fopen(ref_file, "r");
    if (!fr) { 
        perror("open ref"); 
        for (size_t i = 0; i < total; i++) {
            if (all[i].seq) free(all[i].seq);
        }
        free(all);
        return 1; 
    }
    
    char *ref_seq = malloc(MAXBUF);
    if (!ref_seq) {
        perror("malloc ref_seq");
        fclose(fr);
        for (size_t i = 0; i < total; i++) {
            if (all[i].seq) free(all[i].seq);
        }
        free(all);
        return 1;
    }
    
    if (fscanf(fr, "%s", ref_seq) != 1) {
        fprintf(stderr, "Error reading reference sequence\n");
        fclose(fr);
        free(ref_seq);
        for (size_t i = 0; i < total; i++) {
            if (all[i].seq) free(all[i].seq);
        }
        free(all);
        return 1;
    }
    fclose(fr);
    
    int ref_len = strlen(ref_seq);

    pthread_t *threads = malloc(t * sizeof(*threads));
    ThreadArg *args = malloc(t * sizeof(*args));
    
    if (!threads || !args) {
        perror("malloc threads or args");
        if (threads) free(threads);
        if (args) free(args);
        free(ref_seq);
        for (size_t i = 0; i < total; i++) {
            if (all[i].seq) free(all[i].seq);
        }
        free(all);
        return 1;
    }

    size_t per = (total + t - 1) / t;
    for (int i = 0; i < t; i++) {
        size_t st = i * per;
        size_t ed = (st + per > total) ? total : st + per;

        args[i].records = all;
        args[i].start   = st;
        args[i].end     = ed;
        args[i].ref_seq = ref_seq;
        args[i].ref_len = ref_len;
        args[i].res     = malloc(sizeof(ThreadResult));
        
        if (args[i].res) {
            args[i].res->capacity = 1024;
            args[i].res->size     = 0;
            args[i].res->lines    = malloc(1024 * sizeof(char*));
        }
        
        if (pthread_create(&threads[i], NULL, thread_func, &args[i]) != 0) {
            fprintf(stderr, "Error creating thread %d\n", i);
        }
    }


    FILE *fo = fopen(out_file, "w");
    if (!fo) { 
        perror("open output"); 
        for (int i = 0; i < t; i++) {
            pthread_join(threads[i], NULL);
        }
        free(threads);
        free(args);
        free(ref_seq);
        for (size_t i = 0; i < total; i++) {
            if (all[i].seq) free(all[i].seq);
        }
        free(all);
        return 1; 
    }

    for (int i = 0; i < t; i++) {
        void *ret;
        pthread_join(threads[i], &ret);
        if (ret) {
            ThreadResult *R = (ThreadResult*)ret;
            for (size_t j = 0; j < R->size; j++) {
                fputs(R->lines[j], fo);
                free(R->lines[j]);
            }
            free(R->lines);
            free(R);
        }
    }
    fclose(fo);


    free(ref_seq);
    free(threads);
    for (int i = 0; i < t; i++) {
        if (args[i].res) {
          
        }
    }
    free(args);
    
    for (size_t i = 0; i < total; i++) {
        if (all[i].seq) {
            free(all[i].seq);
        }
    }
    free(all);

    return 0;
}