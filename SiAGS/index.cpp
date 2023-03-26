#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>

void write_actg(std::uint32_t* ref, int pos, char nuc) {
    return;
}

void pack(std::string ref_f, std::uint32_t* ref) {
    char* read{new char[4294967296]{}};   // +4GiB
    std::ifstream infile(ref_f, std::ios::binary);
    infile.read((char*)read, 4294967296);
    infile.close();
    std::ofstream outfile(ref_f+".chr", std::ios::trunc);
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
            write_actg(ref, nuc_e, *(read+i));
            nuc_e += 1;
        }
    }
    outfile.close();
    delete[] read;   // -4GiB
}

void index(std::string ref_f) {
    std::uint32_t* ref{new std::uint32_t[268435456]{}};   // +1GiB
    pack(ref_f, ref);
    delete[] ref;   // -1GiB;
}

int main(int argc, char** argv) {
    std::string ref_f{argv[1]};
    index(ref_f);
    return 0;
}