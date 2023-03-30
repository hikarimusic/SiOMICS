#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <string>


#include <iostream>
#include <bitset>

// void pack(std::string ref_f, std::uint32_t* ref) {
//     char* read{new char[4294967296]{}};   // +4GiB
//     std::ifstream infile;
//     infile.open(ref_f, std::ios::binary);
//     infile.read((char*)read, 4294967296);
//     infile.close();
//     std::ofstream outfile;
//     outfile.open(ref_f+".chr", std::ios::trunc);
//     int read_f{0}, nuc_s{0}, nuc_e{0};
//     std::string chr_name{};
//     for (int i=0; i<4294967296; ++i) {
//         if (!(*(read+i))) {
//             outfile << chr_name << '\n' << nuc_e-nuc_s << '\n';
//             break;
//         }
//         if (*(read+i)=='>') {
//             if (nuc_e-nuc_s)
//                 outfile << chr_name << '\n' << nuc_e-nuc_s << '\n';
//             read_f = 1;
//             nuc_s = nuc_e;
//             chr_name = "";
//         }
//         else if (*(read+i)=='\n' && read_f==1)
//             read_f = 2;
//         else if (read_f==1)
//             chr_name += *(read+i);
//         else if (read_f==2 && *(read+i)!='\n') {
//             write_acgt(ref, nuc_e, *(read+i));
//             nuc_e += 1;
//         }
//     }
//     outfile.close();
//     // outfile.open(ref_f+".seq", std::ios::binary);
//     // outfile.write((char*)ref, sizeof(*ref) * ((nuc_e>>4)+1));
//     // outfile.close();
//     delete[] read;   // -4GiB
// }

void write_acgt(std::uint32_t* seq, uint32_t pos, char nuc) {
    uint32_t box = pos >> 4;
    uint32_t pnt = (0b1111 - pos & 0b1111) << 1;
    if (nuc == 'A' || nuc == 'a') {
        seq[box] &= ~(0b01<<pnt);
        seq[box] &= ~(0b10<<pnt);
    }
    else if (nuc == 'C' || nuc == 'c') {
        seq[box] |= (0b01<<pnt);
        seq[box] &= ~(0b10<<pnt);
    }
    else if (nuc == 'G' || nuc == 'g') {
        seq[box] &= ~(0b01<<pnt);
        seq[box] |= (0b10<<pnt);
    }
    else if (nuc == 'T' || nuc == 't') {
        seq[box] |= (0b01<<pnt);
        seq[box] |= (0b10<<pnt);
    }
    else {
        if (std::rand()&1)
            seq[box] &= ~(0b01<<pnt);
        else
            seq[box] |= (0b01<<pnt);
        if (std::rand()&1)
            seq[box] &= ~(0b10<<pnt);
        else
            seq[box] |= (0b10<<pnt);
    }
}

void read(char* seq_f, std::string& chr, std::uint32_t* seq, std::uint32_t& len) {
    std::ifstream infile;
    infile.open(seq_f, std::ios::binary);
    infile.seekg(0, std::ios::end);
    std::uint32_t seq_l = infile.tellg();
    char* read(new char[seq_l]{});
    infile.seekg(0, std::ios::beg);
    infile.read((char*)read, seq_l);
    infile.close();
    uint32_t read_f{0}, nuc_s{0}, nuc_e{0};
    std::string chr_name{};
    for (uint32_t i=0; i<seq_l; ++i) {
        if (!(read[i])) {
            if (nuc_e-nuc_s)
                chr += chr_name + '\n' + std::to_string(nuc_e-nuc_s) + '\n';
            break;
        }
        if (read[i]=='>') {
            if (nuc_e-nuc_s)
                chr += chr_name + '\n' + std::to_string(nuc_e-nuc_s) + '\n';
            read_f = 1;
            nuc_s = nuc_e;
            chr_name = "";
        }
        else if (read[i]=='\n' && read_f==1)
            read_f = 2;
        else if (read_f==1)
            chr_name += read[i];
        else if (read_f==2 && read[i]!='\n') {
            write_acgt(seq, nuc_e, read[i]);
            nuc_e += 1;
        }
    }
    len = (((nuc_e>>4)+2)<<4);
    for (uint32_t i=nuc_e; i<len; ++i)
        write_acgt(seq, i, 't');
    delete[] read;

} 

std::uint32_t nucs(std::uint32_t* seq, std::uint32_t len, std::uint32_t p, std::uint32_t r) {
    uint32_t b = p >> 4;
    uint32_t h = p & 0b1111;
    uint32_t t = h + r;
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

void build(std::uint32_t* seq, std::uint32_t len, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    uint32_t grp[256]{};
    for (uint32_t i=0; i<len; ++i)
        grp[nucs(seq, len, i, 4)] += 1;
    for (uint32_t g=0; g<256; ++g) {
        uint32_t* pmt{new uint32_t[grp[g]]{}};
        uint32_t* tmp{new uint32_t[grp[g]]{}};
        uint32_t* pos{new uint32_t[grp[g]]{}};
        uint32_t* cnt{new uint32_t[256]{}};
        uint32_t com{0};
        uint32_t moc{0};
        uint32_t j{};
        for (uint32_t i=0; i<len; ++i) {
            if (nucs(seq, len, i, 4)==g) {
                pmt[j] = i;
                pos[j] = 4;
                j += 1;
            }
        }
        while (com < grp[g]) {
            for (uint32_t i=0; i<256; ++i)
                cnt[i] = 0;
            for (moc=com; moc<grp[g]; ++moc) {
                if (pos[moc]!=pos[com])
                    break;
                if (nucs(seq, len, pmt[moc]+pos[moc]-4, 4)!=nucs(seq, len, pmt[com]+pos[com]-4, 4))
                    break;
                cnt[nucs(seq, len, pmt[moc]+pos[moc], 4)] += 1;
            }
            bool sorted{1};
            if (cnt[0] > 1)
                sorted = 0;
            for (uint32_t i=1; i<256; ++i) {
                if (cnt[i] > 1)
                    sorted = 0;
                cnt[i] += cnt[i-1];
            }
            for (int64_t i=moc-1; i>=com; --i) 
                tmp[com+(--cnt[nucs(seq, len, pmt[i]+pos[i], 4)])] = pmt[i];
            for (int64_t i=moc-1; i>=com; --i) {
                pmt[i] = tmp[i];
                pos[i] += 4;
            }
            if (sorted)
                com = moc;
        }
        for (int i=0; i<grp[g]; ++i) {
            for (int j=pmt[i]; j<pmt[i]+100 && j<len; ++j) {
                if (nucs(seq, len, j, 1)==0)
                    std::cout << 'A';
                else if (nucs(seq, len, j, 1)==1)
                    std::cout << 'C';
                else if (nucs(seq, len, j, 1)==2)
                    std::cout << 'G';
                else if (nucs(seq, len, j, 1)==3)
                    std::cout << 'T';
            }
            std::cout << '\n';
        }
        delete[] pmt;
        delete[] tmp;
        delete[] pos;
        delete[] cnt;
    }
}

void index(char** argv) {
    std::string chr{};
    std::uint32_t* seq{new std::uint32_t[268435456]{}};   // +1GiB
    std::uint32_t len{};
    read(argv[1], chr, seq, len);
    std::uint32_t* sfa{new std::uint32_t[268435456]{}};   // +1GiB
    std::uint32_t* bwt{new std::uint32_t[268435456]{}};   // +1GiB
    std::uint32_t* occ{new std::uint32_t[268435456]{}};   // +1GiB
    build(seq, len, sfa, bwt, occ);
    // save(chr, seq, bwt, sfa, occ);
    delete[] seq;   // -1GiB
    delete[] sfa;   // -1GiB
    delete[] bwt;   // -1GiB
    delete[] occ;   // -1GiB
}

int main(int argc, char** argv) {
    index(argv);
    return 0;
}