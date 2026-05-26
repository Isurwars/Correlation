#include <gtest/gtest.h>
#include "readers/ArcReader.hpp"
#include <fstream>
#include <cstdio>

using namespace correlation::readers;

TEST(ArcReaderTests, Properties) {
    ArcReader reader;
    EXPECT_EQ(reader.getName(), "Accelrys ARC");
    EXPECT_TRUE(reader.isTrajectory());
    auto exts = reader.getExtensions();
    EXPECT_EQ(exts.size(), 1);
    EXPECT_EQ(exts[0], "arc");
}

TEST(ArcReaderTests, ReadsTrajectory) {
    std::ofstream out("temp_arc.arc");
    out << "!BIOSYM archive 3\n"
        << "PBC=ON\n"
        << "PBC 10.0 10.0 10.0 90.0 90.0 90.0\n"
        << " -12.345\n" // Energy value
        << "C1      1.00  1.00  1.00 XXXX 1      xx      C    0.000\n"
        << "end\n"
        << "PBC=OFF\n"
        << " 42.0\n" // Energy value
        << "C1      2.00  2.00  2.00 XXXX 1      xx      C    0.000\n"
        << "end\n";
    out.close();

    ArcReader reader;
    auto traj = reader.readTrajectory("temp_arc.arc");
    std::remove("temp_arc.arc");

    EXPECT_EQ(traj.getFrameCount(), 2);

    // Frame 1 check
    const auto &f1 = traj.getFrame(0);
    EXPECT_EQ(f1.atomCount(), 1);
    EXPECT_DOUBLE_EQ(f1.lattice_parameters()[0], 10.0);
    EXPECT_DOUBLE_EQ(f1.getEnergy(), -12.345);
    EXPECT_EQ(f1.atoms()[0].element().symbol, "C");

    // Frame 2 check
    const auto &f2 = traj.getFrame(1);
    EXPECT_EQ(f2.atomCount(), 1);
    EXPECT_DOUBLE_EQ(f2.lattice_parameters()[0], 100.0); // PBC=OFF sets 100.0
    EXPECT_DOUBLE_EQ(f2.getEnergy(), 42.0);
    EXPECT_EQ(f2.atoms()[0].element().symbol, "C");
}

TEST(ArcReaderTests, ThrowsOnReadStructure) {
    ArcReader reader;
    EXPECT_THROW(reader.readStructure("dummy.arc"), std::runtime_error);
}
