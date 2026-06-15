#pragma once

#include <cstddef>

namespace correlation::analysis {

struct MaxFrames {
  int value;
};

struct StartFrame {
  size_t value;
};

struct EndFrame {
  size_t value;
};

} // namespace correlation::analysis
