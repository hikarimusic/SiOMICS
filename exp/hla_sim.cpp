#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

std::string revcom(std::string s) {
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

void generate(std::string gen, std::string seq, std::ofstream& outfile1, std::ofstream& outfile2) {
    std::string sim1;
    std::string sim2;
    std::int32_t inl{}, pos1{}, pos2{};
    for (std::int32_t i=0; i<100; ++i) {
        inl = std::min((std::rand()%600)+201, (int) seq.length());
        pos1 = std::rand()%(seq.size()-inl+1);
        pos2 = pos1+inl-100;
        sim1 += ("@"+gen+" "+std::to_string(i)+"/1\n");
        sim2 += ("@"+gen+" "+std::to_string(i)+"/2\n");
        if (std::rand()&1) {
            sim1 += (seq.substr(pos1, 100)+'\n');
            sim2 += (revcom(seq.substr(pos2, 100))+'\n');
        }
        else {
            sim1 += (revcom(seq.substr(pos2, 100))+'\n');;
            sim2 += (seq.substr(pos1, 100)+'\n');
        }
        sim1 += "+\nEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n";
        sim2 += "+\nEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n";
    }
    outfile1 << sim1;
    outfile2 << sim2;
}

void hla_sim(char** argv) {
    std::ifstream infile;
    infile.open(argv[1], std::ios::binary);
    infile.seekg(0, std::ios::end);
    std::uint32_t seq_l = infile.tellg();
    char* read(new char[seq_l+1]{});
    infile.seekg(0, std::ios::beg);
    infile.read((char*)read, seq_l);
    infile.close();
    std::ofstream outfile1;
    std::ofstream outfile2;
    outfile1.open(std::string(argv[1])+".1.fastq", std::ios::trunc);
    outfile2.open(std::string(argv[1])+".2.fastq", std::ios::trunc);
    std::uint32_t read_f{0}, nuc_s{0}, nuc_e{0};
    std::string gen;
    std::string seq;
    for (std::uint32_t i=0; i<=seq_l; ++i) {
        if (!(read[i])) {
            if (nuc_e-nuc_s)
                generate(gen, seq, outfile1, outfile2);
            break;
        }
        if (read[i]=='>') {
            if (nuc_e-nuc_s) 
                generate(gen, seq, outfile1, outfile2);
            nuc_s = nuc_e;
            read_f = 1;
            gen = "";
            seq = "";
        }
        else if (read[i]=='\n' && read_f==1)
            read_f = 2;
        else if (read_f==1)
            gen += read[i];
        else if (read_f==2 && read[i]!='\n') {
            seq += read[i];
            nuc_e += 1;
            if ((nuc_e % 1024)==0)
                std::cout << '\r' << std::flush << "[Process] " << nuc_e << "/" << seq_l << "                    ";
        }
    }
    outfile1.close();
    outfile2.close();
    std::cout << '\r' << std::flush << "[Process] " << "Complete" << "                    " << '\n';
}

int main(int argc, char** argv) {
    hla_sim(argv);
    return 0;
}