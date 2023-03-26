#include <iostream>
#include <fstream>
#include <cstdint>

int main() {
    int n = 2000000000;
    std::uint32_t* arr{new std::uint32_t[n]{1, 2, 3, 4, 5, 6, 7, 8}};
    std::ofstream o("d", std::ios::binary);
    o.write((char*)arr, sizeof(*arr) * n);
    o.close();
    delete[] arr;
    std::uint32_t* res{new std::uint32_t[n]{}};
    std::ifstream i("d", std::ios::binary);
    i.read((char*)res, sizeof(*res) * n);
    i.close();
    // for (int x=0; x<n; ++x)
    //     std::cout << *(res+x) << ' ';
    delete[] res;
    return 0;
}