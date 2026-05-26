#include <gtest/gtest.h>
#include "readers/OnetepDatReader.hpp"
#include <fstream>
#include <cstdio>

using namespace correlation::readers;

TEST(OnetepDatReaderTests, Properties) {
    OnetepDatReader reader;
    EXPECT_EQ(reader.getName(), "ONETEP DAT");
    EXPECT_FALSE(reader.isTrajectory());
    auto exts = reader.getExtensions();
    EXPECT_EQ(exts.size(), 1);
    EXPECT_EQ(exts[0], "dat");
}

TEST(OnetepDatReaderTests, ReadsStructureCartesianAndBohr) {
    std::ofstream out("temp_dat.dat");
    out << "# Comment at the start\n"
        << "%BLOCK LATTICE_CART\n"
        << "bohr\n" // test unit token
        << "10.0 0.0 0.0\n"
        << "0.0 10.0 0.0\n"
        << "0.0 0.0 10.0\n"
        << "%ENDBLOCK LATTICE_CART\n\n"
        << "! Comment before positions\n"
        << "%BLOCK POSITIONS_ABS\n"
        << "angstrom\n" // test unit token
        << "C 1.0 2.0 3.0 # inline comment\n"
        << "%ENDBLOCK POSITIONS_ABS\n";
    out.close();

    OnetepDatReader reader;
    auto cell = reader.readStructure("temp_dat.dat");
    std::remove("temp_dat.dat");

    EXPECT_EQ(cell.atomCount(), 1);
    // Since unit bohr was used for lattice: 10.0 bohr (meaning the cell has those parameters directly, wait, does OnetepDatReader convert bohr to angstrom? No, OnetepDatReader doesn't do a bohr-to-angstrom conversion for lattice/positions, it just skips the token "bohr" or "angstrom"). Let's verify.
    // Yes: in OnetepDatReader.cpp:
    // if (lower_first == "angstrom" || lower_first == "bohr") { continue; }
    // So it skips them. Thus cell parameters remain 10.0.
    EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 10.0);
    EXPECT_EQ(cell.atoms()[0].element().symbol, "C");
    EXPECT_DOUBLE_EQ(cell.atoms()[0].position().x(), 1.0);
}

TEST(OnetepDatReaderTests, ReadsStructureLatticeAbcAndPositionsFrac) {
    std::ofstream out("temp_dat_abc.dat");
    out << "%BLOCK LATTICE_ABC\n"
        << "15.0 15.0 15.0\n"
        << "90.0 90.0 90.0\n"
        << "%ENDBLOCK LATTICE_ABC\n"
        << "%BLOCK POSITIONS_FRAC\n"
        << "H 0.1 0.1 0.1\n"
        << "O 1.2 0.2 -0.1\n" // should wrap to (0.2, 0.2, 0.9) * 15.0 -> (3.0, 3.0, 13.5)
        << "%ENDBLOCK POSITIONS_FRAC\n";
    out.close();

    OnetepDatReader reader;
    auto cell = reader.readStructure("temp_dat_abc.dat");
    std::remove("temp_dat_abc.dat");

    EXPECT_EQ(cell.atomCount(), 2);
    EXPECT_DOUBLE_EQ(cell.lattice_parameters()[0], 15.0);
    EXPECT_EQ(cell.atoms()[0].element().symbol, "H");
    EXPECT_EQ(cell.atoms()[1].element().symbol, "O");

    const auto &p1 = cell.atoms()[0].position();
    EXPECT_NEAR(p1.x(), 1.5, 1e-5);
    EXPECT_NEAR(p1.y(), 1.5, 1e-5);
    EXPECT_NEAR(p1.z(), 1.5, 1e-5);

    const auto &p2 = cell.atoms()[1].position();
    EXPECT_NEAR(p2.x(), 3.0, 1e-5);
    EXPECT_NEAR(p2.y(), 3.0, 1e-5);
    EXPECT_NEAR(p2.z(), 13.5, 1e-5);
}

TEST(OnetepDatReaderTests, ThrowsOnReadTrajectory) {
    OnetepDatReader reader;
    EXPECT_THROW(reader.readTrajectory("dummy.dat"), std::runtime_error);
}
