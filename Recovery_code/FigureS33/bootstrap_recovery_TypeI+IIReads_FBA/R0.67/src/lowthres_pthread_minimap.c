#define _XOPEN_SOURCE 700
#define _POSIX_C_SOURCE 200809L


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>     // mkstemp, unlink, close
#include <errno.h>

#define MAXSEQ 500
#define LINEBUF 1024

typedef struct {
    char seq[151];  // 150bp + NUL
    int  pea, pos, zf;
} PeakRecord;

/* 反转字符串（长度已知，不含末尾NUL） */
static void inv(char* x, int n) {
    for (int i = 0, j = n - 1; i < j; i++, j--) {
        char t = x[i]; x[i] = x[j]; x[j] = t;
    }
}

/* 互补碱基（不负责添加NUL） */
static void complement(const char *in, char *out, int len) {
    for (int i = 0; i < len; i++) {
        switch (in[i]) {
            case 'A': out[i] = 'T'; break;
            case 'T': out[i] = 'A'; break;
            case 'G': out[i] = 'C'; break;
            case 'C': out[i] = 'G'; break;
            case 'a': out[i] = 't'; break;
            case 't': out[i] = 'a'; break;
            case 'g': out[i] = 'c'; break;
            case 'c': out[i] = 'g'; break;
            default:  out[i] = 'N'; break;
        }
    }
}

/*
 * 用 minimap2 比对一个 75bp 片段。
 * 输入：
 *   ref_path: 参考 .mmi 或 .fasta 路径
 *   qseq/qlen: 查询序列与长度（75）
 * 输出：
 *   *editDist: 编辑距离（优先 NM:i:，否则用 aln_len-nmatch 估算）
 *   *tstart  : 目标起点（0-based）
 *   *strand  : '+' 或 '-'
 * 返回：
 *   1 = 有命中并成功解析；0 = 无命中或失败
 */
static int align_with_mm2(const char *ref_path,
                          const char *qseq, int qlen,
                          int *editDist, int *tstart, char *strand)
{
    *editDist = 1000000000;
    *tstart = -1;
    *strand = '+';

    // 1) 写临时 query FASTA
    char qtmp[] = "/tmp/mm2qXXXXXX";
    int qfd = mkstemp(qtmp);
    if (qfd == -1) {
        perror("mkstemp(query)");
        return 0;
    }
    FILE *fq = fdopen(qfd, "w");
    if (!fq) {
        perror("fdopen(query)");
        close(qfd);
        unlink(qtmp);
        return 0;
    }
    // 序列可能不是以NUL结尾，这里用长度限定的写法
    fprintf(fq, ">q\n%.*s\n", qlen, qseq);
    fclose(fq); // 也会关闭 qfd

    // 2) 输出临时文件
    char otmp[] = "/tmp/mm2oXXXXXX";
    int ofd = mkstemp(otmp);
    if (ofd == -1) {
        perror("mkstemp(output)");
        unlink(qtmp);
        return 0;
    }
    close(ofd);

    // 3) 调 minimap2 产生 PAF
    //    -x sr: 短读段预设；--secondary=no: 去次要；-N 1: 只要一个最佳；-c: 输出CIGAR/NM
    char cmd[2048];
    snprintf(cmd, sizeof(cmd),
             "minimap2 -x sr --secondary=no -N 1 -c \"%s\" \"%s\" > \"%s\" 2>/dev/null",
             ref_path, qtmp, otmp);
    int rc = system(cmd);
    if (rc == -1) {
        perror("system(minimap2)");
        unlink(qtmp);
        unlink(otmp);
        return 0;
    }

    // 4) 解析 PAF 第一行
    FILE *fo = fopen(otmp, "r");
    if (!fo) {
        perror("open mm2 output");
        unlink(qtmp);
        unlink(otmp);
        return 0;
    }
    char line[4096];
    int got = 0;
    if (fgets(line, sizeof(line), fo) && line[0] != '\n' && line[0] != '\0') {
        // PAF: 0 qname 1 qlen 2 qstart 3 qend 4 strand 5 tname 6 tlen 7 tstart 8 tend
        //      9 nmatch 10 alnlen 11 mapq ...
        char *tok;
        char *saveptr = NULL;
        int col = 0;
        long nmatch = -1, alnlen = -1;
        int local_tstart = -1;
        char local_strand = '+';
        int local_nm = -1;

        for (tok = strtok_r(line, "\t\n", &saveptr);
             tok;
             tok = strtok_r(NULL, "\t\n", &saveptr), col++)
        {
            if (col == 4) {
                local_strand = tok[0];
            } else if (col == 7) {
                local_tstart = atoi(tok);
            } else if (col == 9) {
                nmatch = atol(tok);
            } else if (col == 10) {
                alnlen = atol(tok);
            } else if (col >= 12) {
                if (strncmp(tok, "NM:i:", 5) == 0) {
                    local_nm = atoi(tok + 5);
                }
            }
        }

        if (local_tstart >= 0 && ((nmatch >= 0 && alnlen > 0) || local_nm >= 0)) {
            *tstart = local_tstart;
            *strand = local_strand;
            if (local_nm >= 0) {
                *editDist = local_nm;
            } else {
                long diff = alnlen - nmatch;
                if (diff < 0) diff = 0;
                *editDist = (int)diff;
            }
            got = 1;
        }
    }
    fclose(fo);

    // 5) 清理
    unlink(qtmp);
    unlink(otmp);

    return got;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr,
            "Usage: %s <reference_mmi_or_fasta> <peak_file> <output_file>\n"
            "Example:\n"
            "  %s scaffold_ref.mmi peaks.txt output.txt\n", argv[0], argv[0]);
        return 1;
    }
    const char *ref_path  = argv[1];  // .mmi 或 .fasta
    const char *peak_file = argv[2];
    const char *out_file  = argv[3];

    // 读取 peak_file（seq 150bp + pea + pos + zf）
    FILE *fp = fopen(peak_file, "r");
    if (!fp) { perror("open peak_file"); return 1; }

    // 先统计条数
    size_t total = 0;
    char tmpseq[MAXSEQ];
    int pea, pos, zf;
    while (fscanf(fp, "%s %d %d %d", tmpseq, &pea, &pos, &zf) == 4) total++;
    rewind(fp);

    PeakRecord *all = (PeakRecord*)malloc(total * sizeof *all);
    if (!all) { fprintf(stderr, "malloc failed for %zu records\n", total); fclose(fp); return 1; }

    for (size_t i = 0; i < total; i++) {
        if (fscanf(fp, "%150s %d %d %d", all[i].seq, &all[i].pea, &all[i].pos, &all[i].zf) != 4) {
            fprintf(stderr, "Input parse error at line %zu\n", i+1);
            free(all); fclose(fp); return 1;
        }
    }
    fclose(fp);

    FILE *fo = fopen(out_file, "w");
    if (!fo) { perror("open output"); free(all); return 1; }

    char before75[76], after75[76], tmp[76], linebuf[LINEBUF];

    for (size_t i = 0; i < total; i++) {
        PeakRecord *r = &all[i];
        if ((int)strlen(r->seq) != 150) continue;  // 跳过非150bp

        // 切成两段
        strncpy(before75, r->seq, 75); before75[75] = '\0';
        strncpy(after75,  r->seq + 75, 75); after75[75] = '\0';

        int ed_bf=1000000000, pos_bf=-1;
        int ed_br=1000000000, pos_br=-1;
        int ed_af=1000000000, pos_af=-1;
        int ed_ar=1000000000, pos_ar=-1;
        char st_bf='+', st_br='+', st_af='+', st_ar='+';

        // 1) before75 正向
        align_with_mm2(ref_path, before75, 75, &ed_bf, &pos_bf, &st_bf);

        // 2) before75 反向互补
        memcpy(tmp, before75, 75);
        inv(tmp, 75); complement(tmp, tmp, 75); tmp[75] = '\0';
        align_with_mm2(ref_path, tmp, 75, &ed_br, &pos_br, &st_br);

        // 3) after75 正向
        align_with_mm2(ref_path, after75, 75, &ed_af, &pos_af, &st_af);

        // 4) after75 反向互补
        memcpy(tmp, after75, 75);
        inv(tmp, 75); complement(tmp, tmp, 75); tmp[75] = '\0';
        align_with_mm2(ref_path, tmp, 75, &ed_ar, &pos_ar, &st_ar);

        // 选择最小编辑距离
        int min_ed = ed_bf, min_pos = pos_bf, choice = 0;
        char min_strand = st_bf;
        if (ed_br < min_ed) { min_ed = ed_br; min_pos = pos_br; min_strand = st_br; choice = 1; }
        if (ed_af < min_ed) { min_ed = ed_af; min_pos = pos_af; min_strand = st_af; choice = 2; }
        if (ed_ar < min_ed) { min_ed = ed_ar; min_pos = pos_ar; min_strand = st_ar; choice = 3; }

        // 误差率（minimap2 无法直接分出 ins/del/sub，保持为0）
        double ins=0.0, del=0.0, sub=0.0;
        double err = (min_ed >= 0 && min_ed < 1000000000) ? (min_ed / 75.0) : 1.0;

        // orient/outpos 与原逻辑一致
        int orient, outpos;
        switch (choice) {
            case 0: orient = (min_strand == '+') ? 0 : 1; outpos = min_pos;         break;
            case 1: orient = 1;                         outpos = min_pos - 75;    break;
            case 2: orient = (min_strand == '+') ? 0 : 1; outpos = min_pos - 75;    break;
            default:orient = 1;                         outpos = min_pos;          break;
        }

        if (min_pos >= 0 && min_pos < 100000000 && err >= 0.0 && err < 0.2) {
            // 输出格式保持不变：seq pea pos zf min_ed outpos err orient
            int L = snprintf(linebuf, sizeof linebuf,
                             "%s %d %d %d %d %d %.6f %d\n",
                             r->seq, r->pea, r->pos, r->zf,
                             min_ed, outpos, err, orient);
            if (L > 0) fputs(linebuf, fo);
        }
    }

    fclose(fo);
    free(all);
    return 0;
}
