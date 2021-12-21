/* ---------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2021 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * ----------------------------------------------------------------------
 */
#include "Smoothing.h"

#include <cmath>

/*
 * Functions over vectors
 */
template<typename T> std::vector<T> DivideVectorByScalar(const std::vector<T>& v, T k) {
  std::vector<T> aux(v.size(), 0);
  for (int i = 0; i < v.size(); i++) {
    aux[i] = v[i] / k;
  }
  return aux;
}  // DivideVectorByScalar

template<typename T> T VectorSum(const std::vector<T>& v) {
  T sum = 0.0;
  for (auto k : v) {
    sum += k;
  }
  return sum;
}  // VectorSum

template<typename T> std::vector<T> NormalizeVector(const std::vector<T>& v) {
  return DivideVectorByScalar(v, VectorSum(v));
}  // NormalizeVector

double sigma2fwhm(double sigma) {
  // Sigma to Full Width Half Maximum
  return sigma * sqrt(8 * log(2));
}

double fwhm2sigma(double fwhm) {
  // Full Width Half Maximum to Sigma
  return fwhm / sqrt(8 * log(2));
}

std::vector<double> GaussianKernel(std::vector<double> r_vals, double r0, double sigma) {
  // Gaussian Kernel
  int    i;
  double a = 1 / (2 * sigma);
  double b = -2 * a * a;
  std::vector<double> aux(r_vals.size(), 0);
  for (i = 0; i < r_vals.size(); i++) {
    aux[i] = a * exp(b * pow(r_vals[i] - r0, 2));
  }
  return aux;
}

std::vector<double> KernelSmoothing(std::vector<double> r,
  std::vector<double>                                   y,
  double                                                sigma,
  int                                                   _kernel_) {
  //  Kernel Smoothing
  int    i, j;
  int    n_ = r.size();
  double kernel_sum;
  std::vector<double> kernel_at_pos(n_, 0);
  std::vector<double> smoothed(n_, 0);
  for (i = 0; i < n_; i++) {
    if (_kernel_ == 1) {
      kernel_at_pos = GaussianKernel(r, r[i], sigma);
      kernel_at_pos = NormalizeVector(kernel_at_pos);
    }
    for (j = 0; j < n_; j++) {
      smoothed[i] += kernel_at_pos[j] * y[j];
    }
  }
  return smoothed;
}  //  Kernel Smoothing
