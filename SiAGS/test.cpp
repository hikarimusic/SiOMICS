#include <iostream>
#include <fstream>
#include <cstdint>

int main() {
    uint32_t a{1299936};
    if (!(a & (uint32_t) 8192))
        std::cout << (a % 8192);
    return 0;
}