#pragma once

namespace correlation {

#if defined(CORRELATION_USE_SINGLE_PRECISION) && CORRELATION_USE_SINGLE_PRECISION
using real_t = float;
inline constexpr const char *REAL_TYPE_NAME = "float";
inline constexpr bool is_single_precision = true;
#else
using real_t = double;
inline constexpr const char *REAL_TYPE_NAME = "double";
inline constexpr bool is_single_precision = false;
#endif

} // namespace correlation
