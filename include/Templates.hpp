#ifndef INCLUDE_TEMPLATES_HPP_
#define INCLUDE_TEMPLATES_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <algorithm>
#include <numeric>
#include <optional>
#include <vector>

//---------------------------------------------------------------------------//
//------------------------- Functions over Vectors --------------------------//
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//---------------------------- Divide by Scalar -----------------------------//
//---------------------------------------------------------------------------//
template <typename T>
std::vector<T> DivideVectorByScalar(const std::vector<T> &v, T k) {
  std::vector<T> aux(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    aux[i] = v[i] / k;
  }
  return aux;
} // DivideVectorByScalar

//---------------------------------------------------------------------------//
//----------------------------- Sum over Vector -----------------------------//
//---------------------------------------------------------------------------//
template <typename T> T VectorSum(const std::vector<T> &v) {
  return std::accumulate(v.begin(), v.end(), T());
} // VectorSum

//---------------------------------------------------------------------------//
//----------------------------- Normalize Vector ----------------------------//
//---------------------------------------------------------------------------//
template <typename T> std::vector<T> NormalizeVector(const std::vector<T> &v) {
  return DivideVectorByScalar(v, VectorSum(v));
} // NormalizeVector

///--------------------------------------------------------------------------//
//------------------------------ Find in Vector -----------------------------//
//---------------------------------------------------------------------------//
template <typename T>
std::pair<bool, int> findInVector(const std::vector<T> &vec, const T &element) {
  auto it = std::find(vec.begin(), vec.end(), element);
  if (it != vec.end()) {
    return {true, std::distance(vec.begin(), it)};
  } else {
    return {false, -1};
  }
} // findInVector

///--------------------------------------------------------------------------//
//--------------------------- Find Index in Vector --------------------------//
//---------------------------------------------------------------------------//
template <typename T>
constexpr std::optional<std::size_t> findIndex(const std::vector<T> &vec,
                                               const T &value) {
  const auto it = std::find(vec.cbegin(), vec.cend(), value);
  return (it != vec.cend()) ? std::optional(std::distance(vec.cbegin(), it))
                            : std::nullopt;
} // findIndex//---------------------------------------------------------------------------//
//-------------------------------- Methods ----------------------------------//
//---------------------------------------------------------------------------//

template <typename T>
std::vector<T> DivideVectorByScalar(const std::vector<T> &v, T k);

template <typename T> T VectorSum(const std::vector<T> &v);

template <typename T> std::vector<T> NormalizeVector(const std::vector<T> &v);

template <typename T>
std::pair<bool, int> findInVector(const std::vector<T> &vec, const T &element);

template <typename T>
constexpr std::optional<std::size_t> findIndex(const std::vector<T> &vec,
                                               const T &value);

#endif // INCLUDE_TEMPLATES_HPP_
