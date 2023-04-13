#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>

void load_dict(std::map<std::string, std::pair<std::uint32_t, std::uint32_t>>& dict) {
    // GRCh38.p14
    dict = {
        {"A", {29942532, 29945870}},
        {"B", {31353875, 31357179}},
        {"C", {31268749, 31272092}},
        {"DMA", {32948618, 32953097}},
        {"DMB", {32934636, 32941028}},
        {"DOA", {33004182, 33009591}},
        {"DOB", {32812763, 32817002}},
        {"DPA1", {33064569, 33080748}},
        {"DPA2", {33091482, 33093314}},
        {"DPB1", {33075990, 33089696}},
        {"DPB2", {33112516, 33129113}},
        {"DQA1", {32637406, 32655272}},
        {"DQA2", {32741391, 32747198}},
        {"DQB1", {32659467, 32666657}},
        {"DQB2", {32756098, 32763532}},
        {"DRA", {32439887, 32445046}},
        {"DRB1", {32578775, 32589848}},
        {"DRB3", {10000, 10000}}, // not in nc
        {"DRB4", {10000, 10000}}, // not in nc
        {"DRB5", {32517353, 32530287}},
        {"E", {30489509, 30494194}},
        {"F", {29723434, 29738532}},
        {"G", {29826474, 29831021}},
        {"H", {29887573, 29891079}},
        {"HFE", {26087429, 26098343}},
        {"J", {30005971, 30009956}},
        {"K", {29926659, 29929825}},
        {"L", {30259562, 30266951}},
        {"MICA", {31400711, 31415315}},
        {"MICB", {31494918, 31511124}},
        {"N", {30351074, 30352038}},
        {"P", {29800044, 29803079}},
        {"S", {31381569, 31382487}},
        {"T", {29896443, 29898947}},
        {"TAP1", {32845209, 32853704}},
        {"TAP2", {32821831, 32838739}},
        {"U", {29933764, 29934880}},
        {"V", {29791906, 29797807}},
        {"W", {29955834, 29959058}},
        {"Y", {10000, 10000}} // not in nc
    };
}

void hla_val(char** argv) {
    std::ifstream infile;
    infile.open(argv[1], std::ios::binary);
    infile.seekg(0, std::ios::end);
    std::uint32_t seq_l = infile.tellg();
    char* read(new char[seq_l+1]{});
    infile.seekg(0, std::ios::beg);
    infile.read((char*)read, seq_l);
    infile.close();
    std::map<std::string, std::pair<std::uint32_t, std::uint32_t>> dict;
    load_dict(dict);
    std::map<std::string, std::pair<std::uint32_t, std::uint32_t>> val;
    for (auto x : dict)
        val[x.first] = {0, 0};
    std::int32_t read_f{0}, star{0}, tab{0};
    std::string cls;
    std::string pos;
    for (std::uint32_t i=0; i<=seq_l; ++i) {
        if (read[i]=='\n') {
            if (cls!="" && pos!="") {
                val[cls].second += 1;
                if (std::stol(pos)>dict[cls].first-5000 && std::stol(pos)<dict[cls].second+5000)
                    val[cls].first += 1;
            }
            read_f = 1;
            star = 0;
            tab = 0;
            cls = "";
            pos = "";
        }
        else if (read[i]=='@')
            read_f = 0;
        else if (read_f==1 && read[i]=='*')
            star += 1;
        else if (read_f==1 && read[i]=='\t')
            tab += 1;
        else if (read_f==1) {
            if (star<1)
                cls += read[i];
            if (tab==3)
                pos += read[i];
        }
        if ((i % 1024)==0)
            std::cout << '\r' << std::flush << "[Process] " << i << "/" << seq_l << "                    ";
    }

    for (auto x : val)
        std::cout << '\n' << x.first << ':' << x.second.first << ' ' << x.second.second << '\n';

    delete[] read;
    std::cout << '\r' << std::flush << "[Process] " << "Complete" << "                    " << '\n';
}

int main(int argc, char** argv) {
    hla_val(argv);
    return 0;
}