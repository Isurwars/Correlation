#include <gtest/gtest.h>
#include "readers/CP2KReader.hpp"
#include <fstream>
#include <cstdio>

using namespace correlation::readers;

class CP2KReaderTests : public ::testing::Test {
protected:
    void SetUp() override {
        std::ofstream out("test_cp2k.restart");
        out << "&CELL\n";
        out << "  ABC 12.0 12.0 12.0\n";
        out << "&END CELL\n";
        out << "&COORD\n";
        out << "  H 0.0 0.0 0.0\n";
        out << "  O 1.0 1.0 1.0\n";
        out << "&END COORD\n";
        out.close();
    }
    void TearDown() override {
        std::remove("test_cp2k.restart");
    }
};

TEST_F(CP2KReaderTests, ReadsSingleFrame) {
    CP2KReader reader;
    auto cell = reader.readStructure("test_cp2k.restart");
    EXPECT_EQ(cell.atomCount(), 2);
    EXPECT_EQ(cell.atoms()[0].element().symbol, "H");
    EXPECT_EQ(cell.lattice_parameters()[0], 12.0);
}
