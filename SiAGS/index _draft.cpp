#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <string>
#include <ctime>

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
            if ((nuc_e % 65536)==0)
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
            // std::cout << '\r' << "[Build Index] " << "Com: " << com << " / Num: " << nxt[com]-com << "                    ";
            // std::cout << '\n' << "[Build Index] " << "Pos: " << pos[com] << " / Num: " << nxt[com]-com << "                    ";
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
                // while (nucs(seq, len, pmt[com]+head, 16)==nucs(seq, len, pmt[com+1]+head, 16))
                //     head += 16;
                // while (nucs(seq, len, pmt[com]+head, 16)==nucs(seq, len, pmt[com+1]+head, 16)) {
                //     head += 16;
                // }
                // std::cout << nucs(seq, len, pmt[com]+head, 16) << " " << nucs(seq, len, pmt[com+1]+head, 16) << '\n';
                // std::cout << "Hello\n" ;
                // if (nucs(seq, len, pmt[com]+head, 16)>nucs(seq, len, pmt[com+1]+head, 16)) {
                //     tmp[com] = pmt[com];
                //     pmt[com] = pmt[com+1];
                //     pmt[com+1] = tmp[com];
                // }
                // std::cout << '\n';
                // for (int j=pmt[com]; j<pmt[com]+100 && j<len; ++j) {
                //     if (nucs(seq, len, j, 1)==0)
                //         std::cout << 'A';
                //     else if (nucs(seq, len, j, 1)==1)
                //         std::cout << 'C';
                //     else if (nucs(seq, len, j, 1)==2)
                //         std::cout << 'G';
                //     else if (nucs(seq, len, j, 1)==3)
                //         std::cout << 'T';
                // }
                // std::cout << '\n';
                // // std::cout << pmt[com+1] << '\n';
                // for (int j=pmt[com+1]; j<pmt[com+1]+100 && j<len; ++j) {
                //     if (nucs(seq, len, j, 1)==0)
                //         std::cout << 'A';
                //     else if (nucs(seq, len, j, 1)==1)
                //         std::cout << 'C';
                //     else if (nucs(seq, len, j, 1)==2)
                //         std::cout << 'G';
                //     else if (nucs(seq, len, j, 1)==3)
                //         std::cout << 'T';
                // }
                // std::cout << pmt[com+1] << '\n';
                // com += 2;
                // while (nxt[com]-com==1)
                //     com += 1;
                // while (com+1 == nxt[com]) {
                //     std::cout << '\n';
                //     for (int j=pmt[com]; j<pmt[com]+100 && j<len; ++j) {
                //         if (nucs(seq, len, j, 1)==0)
                //             std::cout << 'A';
                //         else if (nucs(seq, len, j, 1)==1)
                //             std::cout << 'C';
                //         else if (nucs(seq, len, j, 1)==2)
                //             std::cout << 'G';
                //         else if (nucs(seq, len, j, 1)==3)
                //             std::cout << 'T';
                //     }
                //     com += 1;
                // }
                // continue;
            }
            // std::cout <<'\r' << std::flush << "[Build Index] Counting" << "                    "; 
            for (uint32_t i=0; i<gid; ++i)
                cnt[i] = 0;
            // std::cout <<'\r' << std::flush << "[Build Index] Com Moc" << "                    "; 
            // for (moc=com; moc<grp[g]; ++moc) {
            //     if (pos[moc]!=pos[com])
            //         break;
            //     if (nucs(seq, len, pmt[moc]+pos[moc]-4, 4)!=nucs(seq, len, pmt[com]+pos[com]-4, 4))
            //         break;
            //     cnt[nucs(seq, len, pmt[moc]+pos[moc], 4)] += 1;
            // }
            // for (moc=com; moc<nxt[com]; ++moc)
            //     cnt[nucs(seq, len, pmt[moc]+pos[moc], 4)] += 1;
            for (moc=com; moc<nxt[com]; ++moc) {
                cnt[nucs(seq, len, pmt[moc]+pos[moc], dig)] += 1; 
            }
            // std::cout <<'\r' << std::flush << "[Build Index] " << "Com: " << com << " / Num: " << nxt[com]-com << "                    ";
            // bool sorted{1};
            // std::cout <<'\r' << std::flush << "[Build Index] Culmulate" << "                    "; 
            // if (cnt[0] > 1)
            //     sorted = 0;
            for (uint32_t j=0; j<cnt[0]; ++j)
                nxt[com+j] = com + cnt[0];
            for (uint32_t i=1; i<gid; ++i) {
                // if (cnt[i] > 1)
                //     sorted = 0;
                cnt[i] += cnt[i-1];
                for (uint32_t j=cnt[i-1]; j<cnt[i]; ++j)
                    nxt[com+j] = com + cnt[i];
            }
            // std::cout <<'\r' << std::flush << "[Build Index] Build Tmp" << "                    "; 
            for (int64_t i=moc-1; i>=com; --i) 
                tmp[com+(--cnt[nucs(seq, len, pmt[i]+pos[i], dig)])] = pmt[i];
            // std::cout <<'\r' << std::flush << "[Build Index] Build Pmt" << "                    "; 
            for (int64_t i=moc-1; i>=com; --i) {
                pmt[i] = tmp[i];
                pos[i] += dig;
            }
            // if (sorted)
            //     com = moc;
            //     std::cout <<'\r' << std::flush << "[Build Index] (" << g << "/256) " << com << "/" << grp[g] << "                    \n"; 
            while (com+1 == nxt[com])
                com += 1;
            // while (com+1 == nxt[com]) {
            //     std::cout << '\n';
            //     for (int j=pmt[com]; j<pmt[com]+100 && j<len; ++j) {
            //         if (nucs(seq, len, j, 1)==0)
            //             std::cout << 'A';
            //         else if (nucs(seq, len, j, 1)==1)
            //             std::cout << 'C';
            //         else if (nucs(seq, len, j, 1)==2)
            //             std::cout << 'G';
            //         else if (nucs(seq, len, j, 1)==3)
            //             std::cout << 'T';
            //     }
            //     com += 1;
            // }
            // std::cout << '\n';
            if (com % 65536 == 0)
                std::cout <<'\r' << std::flush << "[Build Index] (" << g << "/256) " << com << "/" << grp[g] << "                    "; 
        }
        // for (int i=0; i<grp[g]; ++i) {
        //     for (int j=pmt[i]; j<pmt[i]+100 && j<len; ++j) {
        //         if (nucs(seq, len, j, 1)==0)
        //             std::cout << 'A';
        //         else if (nucs(seq, len, j, 1)==1)
        //             std::cout << 'C';
        //         else if (nucs(seq, len, j, 1)==2)
        //             std::cout << 'G';
        //         else if (nucs(seq, len, j, 1)==3)
        //             std::cout << 'T';
        //     }
        //     std::cout << '\n';
        // }
        delete[] pmt;
        delete[] tmp;
        delete[] pos;
        delete[] cnt;
        std::cout <<'\r' << std::flush << "[Build Index] (" << g << "/256) Complete" << "                    " << '\n'; 
        // std::cout <<'\n' << "[Build Index] (" << g << "/256) Complete" << "                    " ; 
    }
}

void index(char** argv) {
    std::time_t start, finish;
    time(&start);
    std::string chr{};
    std::uint32_t* seq{new std::uint32_t[268435456]{}};   // +1GiB
    std::uint32_t len{};
    read(argv[1], chr, seq, len);
    // std::cout << chr;
    // return;
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