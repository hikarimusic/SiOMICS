#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <string>
#include <ctime>

#include <iostream>
#include <bitset>

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
    char* read(new char[seq_l+1]{});
    infile.seekg(0, std::ios::beg);
    infile.read((char*)read, seq_l);
    infile.close();
    uint32_t read_f{0}, nuc_s{0}, nuc_e{0};
    std::string chr_name{};
    for (uint32_t i=0; i<=seq_l; ++i) {
        if (!(read[i])) {
            if (nuc_e-nuc_s)
                chr += chr_name + '\n' + std::to_string(nuc_e-nuc_s) + '\n';
            break;
        }
        if (read[i]=='>') {
            if (nuc_e-nuc_s)
                chr += chr_name + '\n' + std::to_string(nuc_e-nuc_s) + '\n';
                nuc_s = nuc_e;
                chr_name = "";
            if (read[i+1]=='N' && read[i+2]=='C') {
                read_f = 1;
            }
            else
                read_f = 0;
        }
        else if (read[i]=='\n' && read_f==1)
            read_f = 2;
        else if (read_f==1)
            chr_name += read[i];
        else if (read_f==2 && read[i]!='\n') {
            write_acgt(seq, nuc_e, read[i]);
            nuc_e += 1;
            if ((nuc_e % 1024)==0)
                std::cout << '\r' << std::flush << "[Read Sequence] " << nuc_e << "/" << seq_l << "                    ";
        }
    }
    len = (((nuc_e>>4)+2)<<4);
    for (uint32_t i=nuc_e; i<len; ++i)
        write_acgt(seq, i, 't');
    delete[] read;
    std::cout << '\r' << std::flush << "[Read Sequence] " << "Complete" << "                    " << '\n';
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
    const uint32_t dig{4};
    const uint32_t gid{256};
    uint32_t grp[256]{};
    for (uint32_t i=0; i<len; ++i)
        grp[nucs(seq, len, i, 4)] += 1;
    uint32_t pgs{0};
    for (uint32_t g=0; g<256; ++g) {
        uint32_t* pmt{new uint32_t[grp[g]]{}};
        uint32_t* tmp{new uint32_t[grp[g]]{}};
        uint32_t* pos{new uint32_t[grp[g]]{}};
        uint32_t* nxt{new uint32_t[grp[g]]{}};
        uint32_t* cnt{new uint32_t[gid]{}};
        uint32_t com{0};
        uint32_t moc{0};
        uint32_t id{};
        for (uint32_t i=0; i<len; ++i) {
            if (nucs(seq, len, i, 4)==g) {
                pmt[id] = i;
                pos[id] = 4;
                nxt[id] = grp[g];
                id += 1;
            }
        }
        while (com < grp[g]) {
            if (nxt[com]-com>=2 && nxt[com]-com<=16) {
                uint32_t head{pos[com]};
                uint32_t num{nxt[com]-com};
                bool flag{1};
                while (flag) {
                    for (uint32_t i=1; i<num; ++i) {
                        if (nucs(seq, len, pmt[com]+head, 16)!=nucs(seq, len, pmt[com+i]+head, 16))
                            flag = 0;
                    }
                    if (flag)
                        head += 16;
                }
                for (uint32_t i=0; i<num; ++i)
                    pos[com+i] = head;
            }
            for (uint32_t i=0; i<gid; ++i)
                cnt[i] = 0;
            for (moc=com; moc<nxt[com]; ++moc) {
                cnt[nucs(seq, len, pmt[moc]+pos[moc], dig)] += 1; 
            }
            for (uint32_t j=0; j<cnt[0]; ++j)
                nxt[com+j] = com + cnt[0];
            for (uint32_t i=1; i<gid; ++i) {
                cnt[i] += cnt[i-1];
                for (uint32_t j=cnt[i-1]; j<cnt[i]; ++j)
                    nxt[com+j] = com + cnt[i];
            }
            for (int64_t i=moc-1; i>=com; --i) 
                tmp[com+(--cnt[nucs(seq, len, pmt[i]+pos[i], dig)])] = pmt[i];
            for (int64_t i=moc-1; i>=com; --i) {
                pmt[i] = tmp[i];
                pos[i] += dig;
            }
            while (com+1 == nxt[com])
                com += 1;
            if (com % 1024 == 0)
                std::cout <<'\r' << std::flush << "[Build Index] " << com+pgs << "/" << len << "                    "; 
        }
        delete[] pmt;
        delete[] tmp;
        delete[] pos;
        delete[] cnt;
        pgs += com;
    }
    std::cout << '\r' << std::flush << "[Build Index] " << "Complete" << "                    " << '\n'; 
}

void index(char** argv) {
    std::time_t start, finish;
    time(&start);
    std::string chr{};
    std::uint32_t* seq{new std::uint32_t[268435456]{}};   // +1GiB
    std::uint32_t len{};
    read(argv[1], chr, seq, len);
    // std::uint32_t* sfa{new std::uint32_t[268435456]{}};   // +1GiB
    // std::uint32_t* bwt{new std::uint32_t[268435456]{}};   // +1GiB
    // std::uint32_t* occ{new std::uint32_t[268435456]{}};   // +1GiB
    std::uint32_t* sfa{new std::uint32_t[1]{}};   // +1GiB
    std::uint32_t* bwt{new std::uint32_t[1]{}};   // +1GiB
    std::uint32_t* occ{new std::uint32_t[1]{}};   // +1GiB
    build(seq, len, sfa, bwt, occ);
    // save(chr, seq, bwt, sfa, occ);
    delete[] seq;   // -1GiB
    delete[] sfa;   // -1GiB
    delete[] bwt;   // -1GiB
    delete[] occ;   // -1GiB
    time(&finish);
    std::cout << "[Finish] Total time: " << difftime(finish, start) << " seconds\n";
}

int main(int argc, char** argv) {
    index(argv);
    return 0;
}