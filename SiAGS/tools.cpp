#include <algorithm>
#include <cstdint>
#include <ctime>
#include <deque>
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
    if (c=='A')
        return 0;
    else if (c=='C')
        return 1;
    else if (c=='G')
        return 2;
    else if (c=='T')
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
        for (std::uint32_t i=((r>>6)<<6); i<r; ++i) {
            if (nucs(bwt, len, i, 1)==c)
                ans += 1;
        }
    }
    if (c==3 && r<=sfa[len/16])
        ans += 1;
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
        row = lfm(row, nucs(bwt, len, row, 1), len, sfa, bwt, occ);
    }
    ans += off;
    return ans;
}

std::string rcseq(std::string s) {
    std::string res;
    for (std::int32_t i=s.size()-1; i>=0; --i) {
        if (s[i]=='A')
            res += 'T';
        else if (s[i]=='C')
            res += 'G';
        else if (s[i]=='G')
            res += 'C';
        else if (s[i]=='T')
            res += 'A';
    }
    return res;
}

void capitalize(std::string& s) {
    for (std::int32_t i=0; i<s.size(); ++i) {
        if (s[i]=='a')
            s[i] = 'A';
        else if (s[i]=='c')
            s[i] = 'C';
        else if (s[i]=='g')
            s[i] = 'G';
        else if (s[i]=='t')
            s[i] = 'T';
    }
}

struct seed {
    std::int64_t qs{};
    std::int64_t qe{};
    std::int64_t rs{};
    std::int64_t re{};
    std::int64_t fg{};
};

struct cluster {
    std::int64_t fg{};
    std::int64_t nc{};
    std::vector<seed> seds;
};

void map(std::string seq1, std::string seq2, std::uint32_t len, std::uint32_t* seq, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::cout << "[Map] ";
    capitalize(seq1);
    capitalize(seq2);
    const std::int64_t s_wid{10};
    const std::int64_t gap_1{50};
    const std::int64_t gap_2{1000};
    std::vector<std::string> seqs{seq1, rcseq(seq1), seq2, rcseq(seq2)};
    std::vector<seed> seds;
    std::deque<std::int64_t> sedt(s_wid);
    for (std::int64_t fg=0; fg<4; ++fg) {
        std::string qry = seqs[fg];
        for (std::int64_t qe=seqs[fg].size(); qe>0; --qe) {
            // std::cout << "qe: " << qe << '\n';
            std::int64_t head = 0;
            std::int64_t tail = len;
            for (std::int64_t qp=qe-1; qp>=0; --qp) {
                head = lfm(head, cti(qry[qp]), len, sfa, bwt, occ);
                tail = lfm(tail, cti(qry[qp]), len, sfa, bwt, occ);
                if ((tail-head)<=s_wid) {
                    // std::cout << "qp: " << qp << '\n';
                    for (std::int64_t rps=head; rps<tail; ++rps) {
                        std::int64_t rp = rpm(rps, len, sfa, bwt, occ);
                        if (std::find(sedt.begin(), sedt.end(), rp+qry.size()-qp)==sedt.end()) {
                            std::int64_t qs{}; 
                            for (qs=qp-1; qs>=0; --qs) {
                                if (qry[qs]!=itc(nucs(seq, len, rp-qp+qs, 1)))
                                    break;
                            }
                            qs += 1;
                            // std::cout << "qs: " << qs << '\n';
                            seds.push_back({qs, qe, rp-qp+qs, rp+qe-qp, fg});
                            // std::cout << qs << ' ' << qe << ' ' << rp-qp+qs << ' ' << rp+qe-qp << ' ' << fg << '\n';
                        }
                        sedt.push_back(rp+qry.size()-qp);
                    }
                    for (std::int64_t i=0; i<s_wid-(tail-head); ++i)
                        sedt.push_back(0);
                    for (std::int64_t i=0; i<s_wid; ++i)
                        sedt.pop_front();
                    // for (auto x : sedt)
                    //     std::cout << x << ' ';
                    // std::cout << '\n';
                    break;
                }
            }
        }
    }
    // std::cout << "seds size: " << seds.size() << '\n';
    std::sort(seds.begin(), seds.end(), [](seed a, seed b){return (a.rs<b.rs) ? 1 : 0;});
    // for (int i=1; i<seds.size(); ++i) {
    //     if (seds[i].rs-seds[i-1].re < 100)
    //         std::cout << seds[i].rs-seds[i-1].re << ' ';
    // }
    // std::cout << '\n';
    // std::int64_t cum{seds[0].re-seds[0].rs};
    // for (std::int64_t i=1; i<seds.size(); ++i) {
    //     if (seds[i].rs<seds[i-1].rs+1000) {
    //         cum += (std::max(seds[i].re-seds[i-1].re, (std::int64_t) 0) - std::max(seds[i].rs-seds[i-1].re, (std::int64_t) 0));
    //         // std::cout << seds[i-1].rs << ' ' << seds[i-1].re << ' ' << seds[i].rs << ' ' << seds[i].re << " cum: " << cum << '\n';
    //     }
    //     else {
    //         std::cout << cum << '\n';
    //         cum = seds[i].re - seds[i].rs;
    //     }
    // }
    std::vector<cluster> cls;
    int now[4]{-1, -1, -1, -1};
    int chs[4]{-1, -1, -1, -1};
    int val[4]{};
    int l{}, fg{};
    for (int i=0; i<seds.size(); ++i) {
        fg = seds[i].fg;
        if ((now[fg])==-1) {
            cls.push_back({seds[i].fg, seds[i].qe-seds[i].qs, {{seds[i].qs, seds[i].qe, seds[i].rs, seds[i].re, seds[i].fg}}});
            now[fg] = cls.size()-1;
            continue;
        }
        if ((seds[i].rs > cls[now[fg]].seds.back().re + gap_1)) {
            if (cls[now[fg]].nc > val[fg]) {
                chs[fg] = now[fg];
                val[fg] = cls[now[fg]].nc;
            }
            for (int j=now[fg]-1; j>=0; --j) {
                if (cls[now[fg]].seds.back().rs>cls[j].seds.back().rs+gap_2)
                    break;
                if (fg>>1 != cls[j].fg>>1) {
                    if (cls[now[fg]].nc+cls[j].nc>val[fg]) {
                        chs[fg] = now[fg];
                        val[fg] = cls[now[fg]].nc+cls[j].nc;
                    }
                    if (cls[now[fg]].nc+cls[j].nc>val[cls[j].fg]) {
                        chs[cls[j].fg] = j;
                        val[cls[j].fg] = cls[now[fg]].nc+cls[j].nc;
                    }
                }
            }
            cls.push_back({seds[i].fg, seds[i].qe-seds[i].qs, {{seds[i].qs, seds[i].qe, seds[i].rs, seds[i].re, seds[i].fg}}});
            now[fg] = cls.size()-1;
        }
        else {
            l = std::min(std::min(seds[i].qe-cls[now[fg]].seds.back().qe, seds[i].re-cls[now[fg]].seds.back().re), seds[i].qe-seds[i].qs);
            if (l>0) {
                cls[now[fg]].nc += l;
                cls[now[fg]].seds.push_back({seds[i].qe-l, seds[i].qe, seds[i].re-l, seds[i].re, seds[i].fg});
            }
        }
    }
    for (int i=0; i<4; ++i) {
        if (cls[now[i]].nc > val[i]) {
            chs[i] = now[i];
            val[i] = cls[now[i]].nc;
        }
        for (int j=now[i]-1; j>=0; --j) {
            if (cls[now[i]].seds.back().rs>cls[j].seds.back().rs+gap_2)
                break;
            if (i>>1 != cls[j].fg>>1) {
                if (cls[now[i]].nc+cls[j].nc>val[i]) {
                    chs[i] = now[i];
                    val[i] = cls[now[i]].nc+cls[j].nc;
                }
                if (cls[now[i]].nc+cls[j].nc>val[cls[j].fg]) {
                    chs[cls[j].fg] = j;
                    val[cls[j].fg] = cls[now[i]].nc+cls[j].nc;
                }
            }
        }
    }
    // if (cls[now[seds.back().fg]].nc > val[seds.back().fg]) {
    //     chs[seds.back().fg] = now[seds.back().fg];
    //     val[seds.back().fg] = cls[now[seds.back().fg]].nc;
    // }
    // for (int i=0; i<4; ++i)
    //     std::cout << now[i] << ' ' << chs[i] << ' ' << val[i] << ' ' << cls[chs[i]].seds.size() << ' ' << cls[chs[i]].seds[0].rs << '\n';
    for (int i=0; i<4; ++i) {
        for (auto x : cls[chs[i]].seds)
            std::cout << x.rs << ' ' << x.re << " / ";
        std::cout << '\n';
    }
    for (int i=0; i<4; ++i) {
        for (auto x : cls[chs[i]].seds)
            std::cout << x.qs << ' ' << x.qe << " / ";
        std::cout << '\n';
    }
}

void search(std::string qry, std::uint32_t len, std::uint32_t* sfa, std::uint32_t* bwt, std::uint32_t* occ) {
    std::cout << "[Search] ";
    capitalize(qry);
    std::uint32_t head{0};
    std::uint32_t tail{len};
    for (std::int32_t i=qry.size()-1; i>=0; --i) {
        head = lfm(head, cti(qry[i]), len, sfa, bwt, occ);
        tail = lfm(tail, cti(qry[i]), len, sfa, bwt, occ);
        if (head==tail)
            break;
    }
    for (std::uint32_t i=head; i<tail; ++i) {
        std::cout << rpm(i, len, sfa, bwt, occ) << ' ';
    }
}

void print(std::uint32_t pos, std::uint32_t _len, std::uint32_t len, std::uint32_t* seq) {
    std::cout << "[Print] ";
    for (std::uint32_t i=pos; i<pos+_len; ++i)
        std::cout << itc(nucs(seq, len, i, 1));
}

void revcom(std::string str) {
    capitalize(str);
    std::cout << "[Revcom] ";
    for (int i=str.size()-1; i>=0; --i) {
        std::cout << itc(3-cti(str[i]));
    }
}

void tools(char** argv) {
    std::cout <<'\r' << std::flush << "[Load Index] " << "Profiling" << "                    ";
    std::vector<std::string> chr_n;
    std::vector<std::uint32_t> chr_c;
    std::uint32_t len{};
    profile(argv[1], chr_n, chr_c, len);
    std::cout <<'\r' << std::flush << "[Load Index] " << "Loading" << "                    ";
    std::uint32_t* seq{new std::uint32_t[len>>4]{}};
    std::uint32_t* sfa{new std::uint32_t[(len>>4)+1]{}};
    std::uint32_t* bwt{new std::uint32_t[len>>4]{}};
    std::uint32_t* occ{new std::uint32_t[(len>>4)+4]{}};
    load(argv[1], len, seq, sfa, bwt, occ);
    std::cout << '\r' << std::flush << "[Load Index] " << "Complete" << "                    " << '\n'; 
    std::string cmd;
    std::cout << '\r' << std::flush << "[Start] " << "Available commands:" << "                    " << '\n';
    std::cout << '\r' << std::flush << "    map [seq1] [seq2] " << "                    " << '\n';
    std::cout << '\r' << std::flush << "    search [str] " << "                    " << '\n';
    std::cout << '\r' << std::flush << "    print [pos] [len] " << "                    " << '\n';
    std::cout << '\r' << std::flush << "    revcom [str] " << "                    " << '\n';
    std::cout << '\r' << std::flush << "    end          " << "                    " << '\n';
    std::cout << '\r' << std::flush << "[Command] ";
    while (std::cin>>cmd) {
        if (cmd=="map" || cmd=="m") {
            std::string seq1;
            std::string seq2;
            std::cin >> seq1 >> seq2;
            map(seq1, seq2, len, seq, sfa, bwt, occ);
        }
        else if (cmd=="search" || cmd=="s") {
            std::string str;
            std::cin >> str;
            search(str, len, sfa, bwt, occ);
        }
        else if (cmd=="print" || cmd=="p") {
            std::uint32_t pos;
            std::uint32_t _len;
            std::cin >> pos >> _len;
            print(pos, _len, len, seq);
        }
        else if (cmd=="revcom" || cmd=="r") {
            std::string str;
            std::cin >> str;
            revcom(str);
        }
        else if (cmd=="end" || cmd=="e")
            break;
        else
            std::cout << "Command not found!";
        std::cout << "\n[Command] ";
    }
    delete[] seq;
    delete[] sfa;
    delete[] bwt;
    delete[] occ;
}

int main(int argc, char** argv) {
    tools(argv);
    return 0;
}