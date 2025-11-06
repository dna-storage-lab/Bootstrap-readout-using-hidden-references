#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include "edlib.h"

void encodeLayer(const std::string &in, std::string &out, int mode) {
    out.resize(in.size());
    for (size_t i = 0; i < in.size(); i++) {
        char c = in[i];
        if (mode == 0) { 
            out[i] = (c=='A' || c=='T') ? c :
                     (c=='G') ? 'A' :
                     (c=='C') ? 'T' : 'N';
        } else {       
            out[i] = (c=='A' || c=='T') ? 'A' :
                     (c=='G' || c=='C') ? 'T' : 'N';
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

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <enc_base_seq.txt> <reads.txt> <out.txt>\n", argv[0]);
        return 1;
    }
    const char* enc_path   = argv[1];
    const char* reads_path = argv[2];
    const char* out_path   = argv[3];


    std::ifstream fin(enc_path);
    if (!fin) {
        perror("Error opening enc_base_seq file");
        return 1;
    }
    std::string enc_seq;
    std::getline(fin, enc_seq);
    fin.close();


    std::string ref_upper, ref_lower;
    encodeLayer(enc_seq, ref_upper, 1);
    encodeLayer(enc_seq, ref_lower, 0);


    std::ifstream fr(reads_path);
    if (!fr) {
        perror("Error opening reads file");
        return 1;
    }
    std::ofstream fo(out_path);
    if (!fo) {
        perror("Error opening output file");
        return 1;
    }

    std::string line;
    while (std::getline(fr, line)) {
        if (line.empty()) continue;

        std::istringstream is(line);
        std::vector<std::string> cols;
        std::string tok;
        for (int i = 0; i < 4 && (is >> tok); i++) {
            cols.push_back(tok);
        }

        std::string rest;
        std::getline(is, rest);
        if (!rest.empty()) {
            if (rest[0]=='\t' || rest[0]==' ') rest.erase(0,1);
            cols.push_back(rest);
        }
        if (cols.empty()) continue;

        const std::string &seq = cols[0];

 
        std::string q_upper, q_lower;
        encodeLayer(seq, q_upper, 1);
        encodeLayer(seq, q_lower, 0);


        auto ru = alignHW(q_upper, ref_upper);
        int su = ru.first;
        auto rl = alignHW(q_lower, ref_lower);
        int sl = rl.first;


        if (std::abs(su - sl) >= 10) continue;

    
        for (size_t i = 0; i < cols.size(); ++i) {
            fo << cols[i] << ' ';
        }
        fo << su << ' ' << sl << '\n';
    }

    fr.close();
    fo.close();
    return 0;
}
