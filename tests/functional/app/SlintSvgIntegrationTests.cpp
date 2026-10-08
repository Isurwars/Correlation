#include <gtest/gtest.h>
#include <slint.h>
#include <span>
#include <string>
#include <vector>

namespace {

class SlintSvgIntegrationTests : public ::testing::Test {
protected:
  void SetUp() override {}
  void TearDown() override {}
};

TEST_F(SlintSvgIntegrationTests, LoadsBasicSvgFromEmbeddedData) {
  std::string svg =
      R"(<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100"><circle cx="50" cy="50" r="40" stroke="green" stroke-width="4" fill="yellow" /></svg>)";
  std::vector<uint8_t> data(svg.begin(), svg.end());
  auto img = slint::private_api::load_image_from_embedded_data(data, "svg");

  EXPECT_EQ(img.size().width, 100);
  EXPECT_EQ(img.size().height, 100);
}

TEST_F(SlintSvgIntegrationTests, LoadsSvgWithGradientFromEmbeddedData) {
  std::string svg =
      R"lit(<svg xmlns="http://www.w3.org/2000/svg" width="120" height="120"><defs><linearGradient id="g1"><stop offset="0%" stop-color="red"/><stop offset="100%" stop-color="blue"/></linearGradient></defs><rect width="120" height="120" fill="url(#g1)"/></svg>)lit";
  std::vector<uint8_t> data(svg.begin(), svg.end());
  auto img = slint::private_api::load_image_from_embedded_data(data, "svg");

  EXPECT_EQ(img.size().width, 120);
  EXPECT_EQ(img.size().height, 120);
}

TEST_F(SlintSvgIntegrationTests, HandlesInvalidSvgFromEmbeddedData) {
  std::string bad_svg = "not_an_svg";
  std::vector<uint8_t> data(bad_svg.begin(), bad_svg.end());
  auto img = slint::private_api::load_image_from_embedded_data(data, "svg");

  // Slint should return a 0x0 image when parsing fails.
  EXPECT_EQ(img.size().width, 0);
  EXPECT_EQ(img.size().height, 0);
}

TEST_F(SlintSvgIntegrationTests, FailsWithDotExtension) {
  std::string svg =
      R"(<svg xmlns="http://www.w3.org/2000/svg" width="50" height="50"><circle cx="25" cy="25" r="20"/></svg>)";
  std::vector<uint8_t> data(svg.begin(), svg.end());
  auto img = slint::private_api::load_image_from_embedded_data(data, ".svg");

  // ".svg" is an invalid extension format string for Slint (must be "svg")
  EXPECT_EQ(img.size().width, 0);
  EXPECT_EQ(img.size().height, 0);
}

TEST_F(SlintSvgIntegrationTests, LoadsSvgFromDataWithoutCacheCollision) {
  std::string svg1 =
      R"(<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100"><rect width="100" height="100" fill="red"/></svg>)";
  std::vector<uint8_t> data1(svg1.begin(), svg1.end());
  auto img1 = slint::Image::load_from_data(data1, "svg");
  EXPECT_EQ(img1.size().width, 100);
  EXPECT_EQ(img1.size().height, 100);

  std::string svg2 =
      R"(<svg xmlns="http://www.w3.org/2000/svg" width="200" height="150"><rect width="200" height="150" fill="blue"/></svg>)";
  std::vector<uint8_t> data2(svg2.begin(), svg2.end());
  auto img2 = slint::Image::load_from_data(data2, "svg");
  EXPECT_EQ(img2.size().width, 200);
  EXPECT_EQ(img2.size().height, 150);
}

TEST_F(SlintSvgIntegrationTests, LoadsSvgFromDataInSeparateThread) {
  std::string svg =
      R"(<svg xmlns="http://www.w3.org/2000/svg" width="300" height="200"><circle cx="150" cy="100" r="50" fill="green"/></svg>)";
  std::vector<uint8_t> data(svg.begin(), svg.end());
  slint::Image img;

  std::jthread worker([&]() { img = slint::Image::load_from_data(data, "svg"); });
  worker.join();

  EXPECT_EQ(img.size().width, 300);
  EXPECT_EQ(img.size().height, 200);
}

} // namespace
