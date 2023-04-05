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

void load(char* seq_f, std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::ifstream infile;
    infile.open(std::string(seq_f)+".seq", std::ios::binary);
    infile.read((char*) seq, len>>2);
    infile.close();
    infile.open(std::string(seq_f)+".sfa", std::ios::binary);
    infile.read((char*) sfa, (len>>2)+4);
    infile.close();
    infile.open(std::string(seq_f)+".bwt", std::ios::binary);
    infile.read((char*) bwt, len>>2);
    infile.close();
    infile.open(std::string(seq_f)+".occ", std::ios::binary);
    infile.read((char*) occ, (len>>2)+16);
    infile.close();
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

std::uint32_t cti(char c) {
    if (c=='A' || c=='a')
        return 0;
    else if (c=='C' || c=='c')
        return 1;
    else if (c=='G' || c=='g')
        return 2;
    else if (c=='T' || c=='t')
        return 3;
    else
        return 0;
}

char itc(std::uint32_t i) {
    if (i==0)
        return 'A';
    else if (i==1)
        return 'C';
    else if (i==2)
        return 'G';
    else if (i==3)
        return 'T';
    else
        return 'A';
}

std::uint32_t lfm(std::uint32_t r, std::uint32_t c, std::uint32_t len, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::uint32_t ans{};
    if ((r>>5)&1) {
        ans = occ[(len>>4)+c] + occ[((r>>6)<<2)+c];
        for (std::uint32_t i=r; i<(((r>>6)+1)<<6); ++i) {
            if (nucs(bwt, len, i, 1)==c)
                ans -= 1;
        }
    }
    else {
        if (r<64)
            ans = occ[(len>>4)+c];
        else
            ans = occ[(len>>4)+c] + occ[(((r>>6)-1)<<2)+c];
        // std::cout << "Hello\n";
        // std::cout << "Hello\n";
        for (std::uint32_t i=((r>>6)<<6); i<r; ++i) {
            // if (r==19)
            //     std::cout << i << ' ' << nucs(bwt, len, i, 1) << ' ' << ans << '\n';
            if (nucs(bwt, len, i, 1)==c)
                ans += 1;
        }
    }
    // if (r==15) {
    //     std::cout << ans << '\n';
    //     std::cout << "sfa: " << sfa[len>>16] << '\n';
    // }
    if (c==3 && r<=sfa[len/16])
        ans += 1;
    // if (r==15)
    //     std::cout << ans << '\n';
    return ans;
}

std::uint32_t rpm(std::uint32_t r, std::uint32_t len, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::uint32_t row{r};
    std::uint32_t ans{};
    std::uint32_t off{};
    for (off=0; off<len; ++off) {
        if (row%16==0) {
            ans = sfa[row/16];
            break;
        }
        else if (row==sfa[len/16]) {
            ans = 0;
            break;
        }
        // std::cout << row << ' ' << bwt[row] << '\n';
        row = lfm(row, nucs(bwt, len, row, 1), len, sfa, bwt, occ);
        // std::cout << row << '\n';
        // std::cout << row << '\n';
    }
    ans += off;
    return ans;
}

void search(std::string qry, std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::uint32_t head{0};
    std::uint32_t tail{len};
    for (std::int64_t i=qry.size()-1; i>=0; --i) {
        head = lfm(head, cti(qry[i]), len, sfa, bwt, occ);
        tail = lfm(tail, cti(qry[i]), len, sfa, bwt, occ);
        // std::cout << head << ' ' << tail << '\n';
    }
    for (std::uint32_t i=head; i<tail; ++i) {
        std::cout << rpm(i, len, sfa, bwt, occ) << ' ';
    }
    std::cout << '\n';
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
    std::uint32_t* sfa{new std::uint32_t[(len>>4)+1]{}};
    std::uint32_t* bwt{new std::uint32_t[len>>4]{}};
    std::uint32_t* occ{new std::uint32_t[(len>>4)+4]{}};
    load(argv[1], len, seq, sfa, bwt, occ);
    search(argv[2], len, seq, sfa, bwt, occ);

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
    // for (int i=0; i<len/16+1; ++i) {
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
    // for (int i=0; i<len/16+4; ++i) {
    //     std::cout << occ[i] << ' ';
    // }
    // std::cout << '\n';

    // std::cout << lfm(1, 1, len, sfa, bwt, occ);
    // return;

    // for (int i=0; i<len; ++i) {
    //     int pos = rpm(i, len, sfa, bwt, occ);
    //     for (int j=pos; j<len; ++j) {
    //         std::cout << itc(nucs(seq, len, j, 1));
    //     }
    //     std::cout << '\n';
    // }

    // std::cout << rpm(15, len, sfa, bwt, occ) << '\n';

    delete[] seq;
    delete[] sfa;
    delete[] bwt;
    delete[] occ;
}

int main(int argc, char** argv) {
    align(argv);
    return 0;
}