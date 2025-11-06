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
    char seq[151];  
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
    char before75[76], after75[76], tmp[76], linebuf[256];

    for (size_t i = A->start; i < A->end; i++) {
        PeakRecord *r = &A->records[i];
        if ((int)strlen(r->seq) != 150) continue;


        strncpy(before75, r->seq, 75); before75[75] = '\0';
        strncpy(after75,  r->seq + 75, 75); after75[75] = '\0';


        int bf_ins=0,bf_del=0,bf_sub=0, ed_bf, pos_bf;
        int br_ins=0,br_del=0,br_sub=0, ed_br, pos_br;
        int af_ins=0,af_del=0,af_sub=0, ed_af, pos_af;
        int ar_ins=0,ar_del=0,ar_sub=0, ed_ar, pos_ar;
        EdlibAlignResult R;


        R = edlibAlign(before75,75, A->ref_seq,A->ref_len,
                       edlibNewAlignConfig(-1,EDLIB_MODE_HW,EDLIB_TASK_PATH,NULL,0));
        ed_bf = R.editDistance; pos_bf = *R.startLocations;
        for (int k=0; k<R.alignmentLength; k++){
            if (R.alignment[k]==1) bf_del++;
            else if (R.alignment[k]==2) bf_ins++;
            else if (R.alignment[k]==3) bf_sub++;
        }
        edlibFreeAlignResult(R);


        strcpy(tmp,before75);
        inv(tmp,75); complement(tmp,tmp,75);
        R = edlibAlign(tmp,75, A->ref_seq,A->ref_len,
                       edlibNewAlignConfig(-1,EDLIB_MODE_HW,EDLIB_TASK_PATH,NULL,0));
        ed_br = R.editDistance; pos_br = *R.startLocations;
        for (int k=0; k<R.alignmentLength; k++){
            if (R.alignment[k]==1) br_del++;
            else if (R.alignment[k]==2) br_ins++;
            else if (R.alignment[k]==3) br_sub++;
        }
        edlibFreeAlignResult(R);


        R = edlibAlign(after75,75, A->ref_seq,A->ref_len,
                       edlibNewAlignConfig(-1,EDLIB_MODE_HW,EDLIB_TASK_PATH,NULL,0));
        ed_af = R.editDistance; pos_af = *R.startLocations;
        for (int k=0; k<R.alignmentLength; k++){
            if (R.alignment[k]==1) af_del++;
            else if (R.alignment[k]==2) af_ins++;
            else if (R.alignment[k]==3) af_sub++;
        }
        edlibFreeAlignResult(R);


        strcpy(tmp,after75);
        inv(tmp,75); complement(tmp,tmp,75);
        R = edlibAlign(tmp,75, A->ref_seq,A->ref_len,
                       edlibNewAlignConfig(-1,EDLIB_MODE_HW,EDLIB_TASK_PATH,NULL,0));
        ed_ar = R.editDistance; pos_ar = *R.startLocations;
        for (int k=0; k<R.alignmentLength; k++){
            if (R.alignment[k]==1) ar_del++;
            else if (R.alignment[k]==2) ar_ins++;
            else if (R.alignment[k]==3) ar_sub++;
        }
        edlibFreeAlignResult(R);


        int min_ed = ed_bf, min_pos = pos_bf, choice = 0;
        if (ed_br < min_ed) { min_ed = ed_br; min_pos = pos_br; choice = 1; }
        if (ed_af < min_ed) { min_ed = ed_af; min_pos = pos_af; choice = 2; }
        if (ed_ar < min_ed) { min_ed = ed_ar; min_pos = pos_ar; choice = 3; }

        double ins, del, sub, err;
        int orient, outpos;
        switch(choice) {
            case 0: ins = bf_ins/75.0; del = bf_del/75.0; sub = bf_sub/75.0; orient = 0; outpos = min_pos;       break;
            case 1: ins = br_ins/75.0; del = br_del/75.0; sub = br_sub/75.0; orient = 1; outpos = min_pos-75; break;
            case 2: ins = af_ins/75.0; del = af_del/75.0; sub = af_sub/75.0; orient = 0; outpos = min_pos-75; break;
            default:ins = ar_ins/75.0; del = ar_del/75.0; sub = ar_sub/75.0; orient = 1; outpos = min_pos;     break;
        }
        err = min_ed / 75.0;

        if (err < 0.2) {
            int L = snprintf(linebuf, sizeof linebuf,
                             "%s %d %d %d %d %d %.6f %d\n",
                             r->seq, r->pea, r->pos, r->zf,
                             min_ed, outpos, err, orient);
            if (res->size + 1 >= res->capacity) {
                res->capacity *= 2;
                res->lines = realloc(res->lines, res->capacity * sizeof(char*));
            }
            char *dup = malloc(L + 1);
            memcpy(dup, linebuf, L+1);
            res->lines[res->size++] = dup;
        }
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
    if (!fp) { perror("open peak"); return 1; }
    size_t total = 0; char tmpseq[MAXSEQ]; int pea, pos, zf;
    while (fscanf(fp, "%s %d %d %d", tmpseq,&pea,&pos,&zf)==4) total++;
    rewind(fp);

    PeakRecord *all = malloc(total * sizeof *all);
    for (size_t i = 0; i < total; i++) {
        fscanf(fp, "%s %d %d %d", all[i].seq, &all[i].pea, &all[i].pos, &all[i].zf);
    }
    fclose(fp);


    FILE *fr = fopen(ref_file, "r");
    if (!fr) { perror("open ref"); return 1; }
    char *ref_seq = malloc(MAXBUF);
    fscanf(fr, "%s", ref_seq);
    fclose(fr);
    int ref_len = strlen(ref_seq);


    pthread_t *threads = malloc(t * sizeof *threads);
    ThreadArg *args   = malloc(t * sizeof *args);

    size_t per = (total + t - 1) / t;
    for (int i = 0; i < t; i++) {
        size_t st = i * per;
        size_t ed = (st + per > total ? total : st + per);

        args[i].records = all;
        args[i].start   = st;
        args[i].end     = ed;
        args[i].ref_seq = ref_seq;
        args[i].ref_len = ref_len;
        args[i].res     = malloc(sizeof(ThreadResult));
        args[i].res->capacity = 1024;
        args[i].res->size     = 0;
        args[i].res->lines    = malloc(1024 * sizeof(char*));

        pthread_create(&threads[i], NULL, thread_func, &args[i]);
    }

    FILE *fo = fopen(out_file, "w");
    if (!fo) { perror("open output"); return 1; }

    for (int i = 0; i < t; i++) {
        void *ret;
        pthread_join(threads[i], &ret);
        ThreadResult *R = (ThreadResult*)ret;
        for (size_t j = 0; j < R->size; j++) {
            fputs(R->lines[j], fo);
            free(R->lines[j]);
        }
        free(R->lines);
        free(R);
    }
    fclose(fo);

    free(all);
    free(ref_seq);
    free(threads);
    free(args);

    return 0;
}
