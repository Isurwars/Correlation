#include <gtest/gtest.h>
#include <slint.h>
#include <string>
#include <span>

class SlintSvgIntegrationTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(SlintSvgIntegrationTests, LoadsBasicSvgFromEmbeddedData) {
    std::string svg = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"100\" height=\"100\"><circle cx=\"50\" cy=\"50\" r=\"40\" stroke=\"green\" stroke-width=\"4\" fill=\"yellow\" /></svg>";
    auto img = slint::private_api::load_image_from_embedded_data(
          std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(svg.data()), svg.size()), "svg");
    
    EXPECT_EQ(img.size().width, 100);
    EXPECT_EQ(img.size().height, 100);
}

TEST_F(SlintSvgIntegrationTests, LoadsSvgWithGradientFromEmbeddedData) {
    std::string svg = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"120\" height=\"120\"><defs><linearGradient id=\"g1\"><stop offset=\"0%\" stop-color=\"red\"/><stop offset=\"100%\" stop-color=\"blue\"/></linearGradient></defs><rect width=\"120\" height=\"120\" fill=\"url(#g1)\"/></svg>";
    auto img = slint::private_api::load_image_from_embedded_data(
          std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(svg.data()), svg.size()), "svg");
    
    EXPECT_EQ(img.size().width, 120);
    EXPECT_EQ(img.size().height, 120);
}

TEST_F(SlintSvgIntegrationTests, HandlesInvalidSvgFromEmbeddedData) {
    std::string bad_svg = "not_an_svg";
    auto img = slint::private_api::load_image_from_embedded_data(
          std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(bad_svg.data()), bad_svg.size()), "svg");
    
    // Slint should return a 0x0 image when parsing fails.
    EXPECT_EQ(img.size().width, 0);
    EXPECT_EQ(img.size().height, 0);
}

TEST_F(SlintSvgIntegrationTests, FailsWithDotExtension) {
    std::string svg = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"50\" height=\"50\"><circle cx=\"25\" cy=\"25\" r=\"20\"/></svg>";
    auto img = slint::private_api::load_image_from_embedded_data(
          std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(svg.data()), svg.size()), ".svg");
    
    // ".svg" is an invalid extension format string for Slint (must be "svg")
    EXPECT_EQ(img.size().width, 0);
    EXPECT_EQ(img.size().height, 0);
}

TEST_F(SlintSvgIntegrationTests, DataUriLoadingFailsAsExpected) {
    // Data URIs are typically not supported natively by slint_image_load_from_path in this version
    std::string b64 = "PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIxMDAiIGhlaWdodD0iMTAwIj48Y2lyY2xlIGN4PSI1MCIgY3k9IjUwIiByPSI0MCIgc3Ryb2tlPSJncmVlbiIgc3Ryb2tlLXdpZHRoPSI0IiBmaWxsPSJ5ZWxsb3ciIC8+PC9zdmc+";
    std::string uri = "data:image/svg+xml;base64," + b64;
    auto img = slint::Image::load_from_path(slint::SharedString(uri));
    
    EXPECT_EQ(img.size().width, 0);
    EXPECT_EQ(img.size().height, 0);
}
