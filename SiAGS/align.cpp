#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void profile(char* seq_f, std::vector<std::string>& chr_n, std::vector<uint32_t>& chr_c, std::uint32_t& len) {
    std::ifstream infile;
    infile.open(std::string(seq_f)+".chr");
    std::string line;
    while (std::getline(infile, line)) {
        chr_n.push_back(line);
        std::getline(infile, line);
        chr_c.push_back((std::uint32_t) std::stol(line));
        len += (std::uint32_t) std::stol(line);
    }
    len = (((len>>6)+1)<<6);
    infile.close();
}

void load(char* seq_f, std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ, std::uint32_t* cco) {
    std::ifstream infile;
    infile.open(std::string(seq_f)+".seq", std::ios::binary);
    infile.read((char*) seq, len>>2);
    infile.close();
    infile.open(std::string(seq_f)+".sfa", std::ios::binary);
    infile.read((char*) sfa, len>>2);
    infile.close();
    infile.open(std::string(seq_f)+".bwt", std::ios::binary);
    infile.read((char*) bwt, len>>2);
    infile.close();
    infile.open(std::string(seq_f)+".occ", std::ios::binary);
    infile.read((char*) occ, len>>2);
    infile.close();
    for (int i=1; i<4; ++i)
        cco[i] = cco[i-1] + occ[(((len>>6)-1)<<2)+i-1];
}

std::uint32_t nucs(std::uint32_t* seq, std::uint32_t len, std::uint32_t p, std::uint32_t r) {
    std::uint32_t b = p >> 4;
    std::uint32_t h = p & 0b1111;
    std::uint32_t t = h + r;
    if (p >= len) {
        return 0;
    }
    else if (p+r >= len) {
        t -= 16;
        return (((seq[b]<<(h<<1))>>(h<<1))<<(t<<1));
    }
    else if (t <= 16) {
        return ((seq[b]<<(h<<1))>>((16-r)<<1));
    }
    else {
        t -= 16;
        return  ((((seq[b]<<(h<<1))>>(h<<1))<<(t<<1))+(seq[b+1]>>((16-t)<<1)));
    }
}

std::uint32_t lfm(std::uint32_t r, std::uint32_t c, std::uint32_t* bwt, std::uint32_t* occ, std::uint32_t* cco) {
    std::uint32_t ans{};
    if ((r>>5)&1) {
        ans = cco[c] + occ[((r>>6)<<2)+c];
        // for
    }
    else {
        if (r<64)
            ans = cco[c] + ((c==3)?1:0);
        else
            ans = cco[c] + occ[(((r>>6)-1)<<2)+c];
        // for 
    }
    return ans;
}

void search(std::string qry, std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ, std::uint32_t* cco) {
    std::uint32_t head{0};
    std::uint32_t tail{len};
    for (std::int32_t i=qry.size()-1; i>=0; --i) {
        return;
    }
}

void align(char** argv) {
    std::vector<std::string> chr_n;
    std::vector<std::uint32_t> chr_c;
    std::uint32_t len{};
    profile(argv[1], chr_n, chr_c, len);

    // for (int i=0; i<chr_n.size(); ++i) {
    //     std::cout << chr_n[i] << ' ' << chr_c[i] << '\n';
    // }
    // std::cout << len << '\n';
    // return;

    std::uint32_t* seq{new std::uint32_t[len>>4]{}};
    std::uint32_t* sfa{new std::uint32_t[len>>4]{}};
    std::uint32_t* bwt{new std::uint32_t[len>>4]{}};
    std::uint32_t* occ{new std::uint32_t[len>>4]{}};
    std::uint32_t cco[4]{};
    load(argv[1], len, seq, sfa, bwt, occ, cco);
    search("ATATTGTGATA", len, seq, sfa, bwt, occ, cco);

    // for (int i=0; i<len; ++i) {
    //     if (nucs(seq, len, i, 1)==0)
    //         std::cout << 'A';
    //     else if (nucs(seq, len, i, 1)==1)
    //         std::cout << 'C';
    //     else if (nucs(seq, len, i, 1)==2)
    //         std::cout << 'G';
    //     else if (nucs(seq, len, i, 1)==3)
    //         std::cout << 'T';
    // }
    // std::cout << '\n';
    // for (int i=0; i<len/16; ++i) {
    //     std::cout << sfa[i] << ' ';
    // }
    // std::cout << '\n';
    // for (int i=0; i<len; ++i) {
    //     if (nucs(bwt, len, i, 1)==0)
    //         std::cout << 'A';
    //     else if (nucs(bwt, len, i, 1)==1)
    //         std::cout << 'C';
    //     else if (nucs(bwt, len, i, 1)==2)
    //         std::cout << 'G';
    //     else if (nucs(bwt, len, i, 1)==3)
    //         std::cout << 'T';
    // }
    // std::cout << '\n';
    for (int i=0; i<len/16; ++i) {
        std::cout << occ[i] << ' ';
    }
    std::cout << '\n';
    for (int i=0; i<4; ++i) {
        std::cout << cco[i] << ' ';
    }
    std::cout << '\n';

    delete[] seq;
    delete[] sfa;
    delete[] bwt;
    delete[] occ;
}

int main(int argc, char** argv) {
    align(argv);
    return 0;
}