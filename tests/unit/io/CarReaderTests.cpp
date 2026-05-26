#include <gtest/gtest.h>
#include "readers/CarReader.hpp"
#include <fstream>
#include <cstdio>

using namespace correlation::readers;

TEST(CarReaderTests, Properties) {
    CarReader reader;
    EXPECT_EQ(reader.getName(), "Accelrys CAR");
    EXPECT_FALSE(reader.isTrajectory());
    auto exts = reader.getExtensions();
    EXPECT_EQ(exts.size(), 1);
    EXPECT_EQ(exts[0], "car");
}

TEST(CarReaderTests, ReadsStructure) {
    std::ofstream out("temp_car.car");
    out << "!BIOSYM archive 3\n"
        << "!Comment line\n"
        << "PBC=ON\n"
        << "PBC 10.5 11.5 12.5 90.0 90.0 90.0\n"
        << "C1      1.00  2.00  3.00 XXXX 1      xx      C    0.000\n"
        << "H1      1.50  2.50  3.50 XXXX 1      xx      H    0.000\n"
        << "end\n"
        << "end\n";
    out.close();

    CarReader reader;
    auto cell = reader.readStructure("temp_car.car");
    std::remove("temp_car.car");

    EXPECT_EQ(cell.atomCount(), 2);
    EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 10.5);
    EXPECT_DOUBLE_EQ(cell.lattice_parameters()[1], 11.5);
    EXPECT_DOUBLE_EQ(cell.lattice_parameters()[2], 12.5);
    EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
    EXPECT_EQ(cell.atoms()[1].element().symbol, "H");
}

TEST(CarReaderTests, ReadsStructurePbcOff) {
    std::ofstream out("temp_car_off.car");
    out << "!BIOSYM archive 3\n"
        << "PBC=OFF\n"
        << "O1      1.00  2.00  3.00 XXXX 1      xx      O    0.000\n"
        << "end\n";
    out.close();

    CarReader reader;
    auto cell = reader.readStructure("temp_car_off.car");
    std::remove("temp_car_off.car");

    EXPECT_EQ(cell.atomCount(), 1);
    EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 100.0); // PBC=OFF sets 100.0
    EXPECT_EQ(cell.atoms()[0].element().symbol, "O");
}

TEST(CarReaderTests, ThrowsOnReadTrajectory) {
    CarReader reader;
    EXPECT_THROW(reader.readTrajectory("dummy.car"), std::runtime_error);
}
