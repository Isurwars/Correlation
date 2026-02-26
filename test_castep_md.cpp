#include "readers/CastepMdReader.hpp"
#include <iostream>

int main() {
    auto frames = CastepMdReader::read("/tmp/nonexistent");
    return 0;
}
