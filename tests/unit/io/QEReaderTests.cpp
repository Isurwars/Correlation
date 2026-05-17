#include <gtest/gtest.h>
#include "readers/QEReader.hpp"
#include <fstream>
#include <cstdio>

using namespace correlation::readers;

class QEReaderTests : public ::testing::Test {
protected:
    void SetUp() override {
        std::ofstream out("test_qe.pwo");
        out << "CELL_PARAMETERS (angstrom)\n";
        out << "10.0 0.0 0.0\n";
        out << "0.0 10.0 0.0\n";
        out << "0.0 0.0 10.0\n";
        out << "ATOMIC_POSITIONS (angstrom)\n";
        out << "C 1.0 2.0 3.0\n";
        out << "O 4.0 5.0 6.0\n";
        out.close();
    }
    void TearDown() override {
        std::remove("test_qe.pwo");
    }
};

TEST_F(QEReaderTests, ReadsSingleFrame) {
    QEReader reader;
    auto cell = reader.readStructure("test_qe.pwo");
    EXPECT_EQ(cell.atomCount(), 2);
    EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
    EXPECT_EQ(cell.lattice_parameters()[0], 10.0);
}
