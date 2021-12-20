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


double sigma2fwhm(double sigma) {
  // Sigma to Full Width Half Maximum
  return sigma * sqrt(8 * log(2));
}

double fwhm2sigma(double fwhm) {
  // Full Width Half Maximum to Sigma
  return fwhm / sqrt(8 * log(2));
}

std::vector<double> GaussianKernel(double r0, std::vector<double> r_vals, double sigma) {
  // Gaussian Kernel
  int    i;
  double a = 1 / (2 * sigma);
  double b = -2 * a * a;
  std::vector<double> aux(r_vals.size(), 0);
  for (i = 0; i < r_vals.size(); i++) {
    aux[i] = a * exp(b * (r_vals[i] - r0) * (r_vals[i] - r0));
  }
  return aux;
}

std::vector<double> KernelSmoothing(std::vector<double> r,
  std::vector<double>                                   y,
  std::vector<double> (*func)(double,
  std::vector<double>,
  double)) {
  //  Kernel Smoothing
  double n_ = r.size();
  double s = 0.1;
  int    i, j;
  std::vector<double> kernel(n_, 0);
  std::vector<double> smoothed(n_, 0);
  for (i = 0; i < n_; i++) {
    kernel = func(r[i], r, s);
    for (j = 0; j < n_; j++) {
      smoothed[i] += kernel[j] + y[j];
    }
  }
  return smoothed;
}  //  Kernel Smoothing
