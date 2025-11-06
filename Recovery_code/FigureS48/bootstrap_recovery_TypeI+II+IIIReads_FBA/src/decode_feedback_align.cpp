#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "edlib.h"

int process_r083(const char* input_file, std::string &out_seq) {
    const int COLUMN = 180 * 360; 
    const int BASE   = 40500;     


    std::vector<char> cblk1(COLUMN);
    {
        std::ifstream fin(input_file);
        if (!fin) {
            std::perror("Error opening dec_results file");
            return -1;
        }
        int ch, cnt = 0;
        while ((ch = fin.get()) != EOF && cnt < COLUMN) {
            if (ch=='0' || ch=='1') {
                cblk1[cnt++] = ch - '0';
            }
        }
        fin.close();
        if (cnt != COLUMN) {
            std::fprintf(stderr, "Error: expected %d bits but got %d\n", COLUMN, cnt);
            return -2;
        }
    }


    std::vector<int> perm(COLUMN);
    {
        std::ifstream fperm("./configure/permutation64800", std::ios::binary);
        if (!fperm) {
            std::perror("Error opening permutation file");
            return -3;
        }
        fperm.read(reinterpret_cast<char*>(perm.data()), COLUMN * sizeof(int));
        fperm.close();
    }


    std::vector<char> cblk(COLUMN);
    for (int i = 0; i < COLUMN; i++) {
        cblk[i] = cblk1[ perm[i] ];
    }


    std::vector<char> sparse_code(BASE * 2);
    for (int j = 0; j < COLUMN / 4; j++) {
        int sum = 0;
        for (int k = 0; k < 4; k++) sum += cblk[j*4 + k];
        if (sum < 3) {
            sparse_code[j*5 + 0] = 0;
            for (int k = 1; k < 5; k++)
                sparse_code[j*5 + k] = cblk[j*4 + (k-1)];
        } else {
            sparse_code[j*5 + 0] = 1;
            for (int k = 1; k < 5; k++)
                sparse_code[j*5 + k] = 1 - cblk[j*4 + (k-1)];
        }
    }


    std::vector<int> watermark(BASE * 2);
    {
        std::ifstream fwm("./configure/SequenceL81000NoPeriodOnly2ndFILE", std::ios::binary);
        if (!fwm) {
            std::perror("Error opening watermark file");
            return -4;
        }
        fwm.read(reinterpret_cast<char*>(watermark.data()), watermark.size() * sizeof(int));
        fwm.close();
    }


    out_seq.resize(BASE);
    for (int j = 0; j < BASE; j++) {
        int b1 = sparse_code[j]         ^ watermark[j];
        int b2 = sparse_code[j + BASE]  ^ watermark[j + BASE];
        char base;
        if (b1==0 && b2==0)      base = 'A';
        else if (b1==0 && b2==1) base = 'T';
        else if (b1==1 && b2==0) base = 'G';
        else                     base = 'C';
        out_seq[j] = base;
    }
    return 0;
}

void encodeLayer(const std::string &in, std::string &out, int mode) {
    out.resize(in.size());
    for (size_t i = 0; i < in.size(); i++) {
        char c = in[i];
        if (mode == 0) { 
            out[i] = (c=='A'||c=='T') ? c :
                     (c=='G')         ? 'A' :
                     (c=='C')         ? 'T' : 'N';
        } else {     
            out[i] = (c=='A'||c=='T') ? 'A' :
                     (c=='G'||c=='C') ? 'T' : 'N';
        }
    }
}


std::pair<int,int> alignHW(const std::string &query, const std::string &ref) {
    EdlibAlignResult r = edlibAlign(
        query.c_str(), query.size(),
        ref.c_str(),   ref.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0)
    );
    int st = r.startLocations[0];
    int ed = r.editDistance;
    edlibFreeAlignResult(r);
    return std::make_pair(st, ed);
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "Usage: %s dec_results.txt stage2_residual_reads output.txt\n", argv[0]);
        return 1;
    }


    std::string enc_seq;
    int ret = process_r083(argv[1], enc_seq);
    if (ret != 0) {
        return 10 + ret;
    }


    std::string ref_up, ref_lo;
    encodeLayer(enc_seq, ref_up, 1);
    encodeLayer(enc_seq, ref_lo, 0);


    std::ifstream fin_reads(argv[2]);
    if (!fin_reads) {
        std::perror("Error opening remaining_reads");
        return 2;
    }
    std::ofstream fout(argv[3]);
    if (!fout) {
        std::perror("Error opening output file");
        return 3;
    }


    std::string line;
    while (std::getline(fin_reads, line)) {
        if (line.empty()) continue;
        std::istringstream is(line);
        std::vector<std::string> cols;
        std::string tok;

        for (int i = 0; i < 4 && (is >> tok); i++) {
            cols.push_back(tok);
        }

        std::string rest;
        std::getline(is, rest);
        if (!rest.empty() && (rest[0]==' '||rest[0]=='\t')) rest.erase(0,1);
        if (!rest.empty()) cols.push_back(rest);
        if (cols.empty()) continue;


        const std::string &seq = cols[0];
        std::string q_up, q_lo;
        encodeLayer(seq, q_up, 1);
        encodeLayer(seq, q_lo, 0);

        auto ru = alignHW(q_up, ref_up);
        auto rl = alignHW(q_lo, ref_lo);
        int su = ru.first, sl = rl.first;


        if (std::abs(su - sl) >= 10) continue;


        for (size_t i = 0; i < cols.size(); ++i) {
            fout << cols[i] << ' ';
        }
        fout << su << ' ' << sl << '\n';
    }

    fin_reads.close();
    fout.close();
    return 0;
}
