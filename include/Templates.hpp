#ifndef SRC_TEMPLATES_H_
#define SRC_TEMPLATES_H_
/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2025 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 * ----------------------------------------------------------------------------
 */
#include <algorithm> // For std::find
#include <cmath>     // For exp and pow
#include <numeric>   // For std::accumulate
#include <optional>  // For std::optional
#include <vector>    // For std::vector

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

#endif // SRC_TEMPLATES_H_
