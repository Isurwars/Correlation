#include <gtest/gtest.h>
#include "readers/LammpsDumpReader.hpp"
#include <fstream>
#include <cstdio>

using namespace correlation::readers;

TEST(LammpsDumpReaderTests, Properties) {
    LammpsDumpReader reader;
    EXPECT_EQ(reader.getName(), "LAMMPS Dump");
    EXPECT_TRUE(reader.isTrajectory());
    auto exts = reader.getExtensions();
    EXPECT_EQ(exts.size(), 2);
    EXPECT_EQ(exts[0], "dump");
    EXPECT_EQ(exts[1], "lammpstrj");
}

TEST(LammpsDumpReaderTests, ReadsTrajectoryOrtho) {
    std::ofstream out("temp_dump.dump");
    out << "ITEM: TIMESTEP\n"
        << "100\n"
        << "ITEM: NUMBER OF ATOMS\n"
        << "1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n"
        << "0.0 10.0\n"
        << "0.0 11.0\n"
        << "0.0 12.0\n"
        << "ITEM: ATOMS id type x y z\n"
        << "1 2 1.0 2.0 3.0\n";
    out.close();

    LammpsDumpReader reader;
    auto traj = reader.readTrajectory("temp_dump.dump");
    std::remove("temp_dump.dump");

    EXPECT_EQ(traj.getFrameCount(), 1);
    const auto &f = traj.getFrame(0);
    EXPECT_EQ(f.atomCount(), 1);
    EXPECT_DOUBLE_EQ(f.lattice_parameters()[0], 10.0);
    EXPECT_DOUBLE_EQ(f.lattice_parameters()[1], 11.0);
    EXPECT_DOUBLE_EQ(f.lattice_parameters()[2], 12.0);
    // Standard atom mapping should fall back to type since no element was specified
    EXPECT_EQ(f.atoms()[0].element().symbol, "2");
}

TEST(LammpsDumpReaderTests, ReadsTrajectoryTriclinicScaledAndElement) {
    std::ofstream out("temp_dump_tri.dump");
    out << "ITEM: TIMESTEP\n"
        << "200\n"
        << "ITEM: NUMBER OF ATOMS\n"
        << "2\n"
        << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n"
        << "0.0 10.0 1.0\n"  // xlo xhi xy
        << "0.0 10.0 2.0\n"  // ylo yhi xz
        << "0.0 10.0 3.0\n"  // zlo zhi yz
        << "ITEM: ATOMS id type element xs ys zs\n"
        << "1 1 C 0.5 0.5 0.5\n"
        << "2 2 H 0.1 0.2 0.3\n";
    out.close();

    LammpsDumpReader reader;
    auto traj = reader.readTrajectory("temp_dump_tri.dump");
    std::remove("temp_dump_tri.dump");

    EXPECT_EQ(traj.getFrameCount(), 1);
    const auto &f = traj.getFrame(0);
    EXPECT_EQ(f.atomCount(), 2);
    
    // Check elements
    EXPECT_EQ(f.atoms()[0].element().symbol, "C");
    EXPECT_EQ(f.atoms()[1].element().symbol, "H");

    // Triclinic lattice:
    // lx = 10, ly = 10, lz = 10
    // a = (lx, 0, 0) = (10, 0, 0)
    // b = (xy, ly, 0) = (1, 10, 0)
    // c = (xz, yz, lz) = (2, 3, 10)
    // For Atom 1 (xs=0.5, ys=0.5, zs=0.5):
    // pos = 0.5*a + 0.5*b + 0.5*c = (5.0, 0.0, 0.0) + (0.5, 5.0, 0.0) + (1.0, 1.5, 5.0) = (6.5, 6.5, 5.0)
    const auto &pos1 = f.atoms()[0].position();
    EXPECT_NEAR(pos1.x(), 6.5, 1e-5);
    EXPECT_NEAR(pos1.y(), 6.5, 1e-5);
    EXPECT_NEAR(pos1.z(), 5.0, 1e-5);
}

TEST(LammpsDumpReaderTests, ThrowsOnInvalidFile) {
    LammpsDumpReader reader;
    EXPECT_THROW(reader.readTrajectory("nonexistent_file.dump"), std::runtime_error);
}
