#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <string>

#include <bitset>

void write_acgt(std::uint32_t* ref, int pos, char nuc) {
    int box = pos >> 4;
    int pnt = (0b1111 - pos & 0b1111) << 1;
    if (nuc == 'A' || nuc == 'a') {
        *(ref+box) &= ~(0b01<<pnt);
        *(ref+box) &= ~(0b10<<pnt);
    }
    else if (nuc == 'C' || nuc == 'c') {
        *(ref+box) |= (0b01<<pnt);
        *(ref+box) &= ~(0b10<<pnt);
    }
    else if (nuc == 'G' || nuc == 'g') {
        *(ref+box) &= ~(0b01<<pnt);
        *(ref+box) |= (0b10<<pnt);
    }
    else if (nuc == 'T' || nuc == 't') {
        *(ref+box) |= (0b01<<pnt);
        *(ref+box) |= (0b10<<pnt);
    }
    else {
        if (std::rand()&1)
            *(ref+box) &= ~(0b01<<pnt);
        else
            *(ref+box) |= (0b01<<pnt);
        if (std::rand()&1)
            *(ref+box) &= ~(0b10<<pnt);
        else
            *(ref+box) |= (0b10<<pnt);
    }
}

void pack(std::string ref_f, std::uint32_t* ref) {
    char* read{new char[4294967296]{}};   // +4GiB
    std::ifstream infile;
    infile.open(ref_f, std::ios::binary);
    infile.read((char*)read, 4294967296);
    infile.close();
    std::ofstream outfile;
    outfile.open(ref_f+".chr", std::ios::trunc);
    int read_f{0}, nuc_s{0}, nuc_e{0};
    std::string chr_name{};
    for (int i=0; i<4294967296; ++i) {
        if (!(*(read+i))) {
            outfile << chr_name << '\n' << nuc_e-nuc_s << '\n';
            break;
        }
        if (*(read+i)=='>') {
            if (nuc_e-nuc_s)
                outfile << chr_name << '\n' << nuc_e-nuc_s << '\n';
            read_f = 1;
            nuc_s = nuc_e;
            chr_name = "";
        }
        else if (*(read+i)=='\n' && read_f==1)
            read_f = 2;
        else if (read_f==1)
            chr_name += *(read+i);
        else if (read_f==2 && *(read+i)!='\n') {
            write_acgt(ref, nuc_e, *(read+i));
            nuc_e += 1;
        }
    }
    outfile.close();
    outfile.open(ref_f+".seq", std::ios::binary);
    outfile.write((char*)ref, sizeof(*ref) * ((nuc_e>>4)+1));
    outfile.close();
    delete[] read;   // -4GiB
}

void index(std::string ref_f) {
    std::uint32_t* ref{new std::uint32_t[268435456]{}};   // +1GiB
    pack(ref_f, ref);
    // for (int i=0; i<20; ++i) {
    //     std::cout << std::bitset<32>(*(ref+i)) << '\n';
    // }
    delete[] ref;   // -1GiB;
}

int main(int argc, char** argv) {
    std::string ref_f{argv[1]};
    index(ref_f);
    return 0;
}