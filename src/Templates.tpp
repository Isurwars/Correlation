/* ----------------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2024 Isaías Rodríguez <isurwars@gmail.com>
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
#include "Templates.hpp"

#include <algorithm> // For std::find
#include <cmath>     // For exp and pow
#include <numeric>   // For std::accumulate

/*
 * Functions over vectors
 */
template <typename T>
std::vector<T> DivideVectorByScalar(const std::vector<T> &v, T k) {
  std::vector<T> aux(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    aux[i] = v[i] / k;
  }
  return aux;
} // DivideVectorByScalar

template <typename T> T VectorSum(const std::vector<T> &v) {
  return std::accumulate(v.begin(), v.end(), T());
} // VectorSum

template <typename T> std::vector<T> NormalizeVector(const std::vector<T> &v) {
  return DivideVectorByScalar(v, VectorSum(v));
} // NormalizeVector

/*
 * Generic function to find if an element of any type exists in a vector,
 * if true, then returns the index.
 */
template <typename T>
std::pair<bool, int> findInVector(const std::vector<T> &vecOfElements,
                                  const T &element) {
  std::pair<bool, int> result;

  // Find given element in vector
  auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

  if (it != vecOfElements.end()) {
    result.second = std::distance(vecOfElements.begin(), it);
    result.first = true;
  } else {
    result.first = false;
    result.second = -1;
  }
  return result;
} // findInVector
