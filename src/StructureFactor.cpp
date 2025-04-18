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
#include "../include/StructureFactor.hpp"

#include <numeric>

#include "../include/Constants.hpp"

//---------------------------------------------------------------------------//
//---------------------------------- Sinc -----------------------------------//
//---------------------------------------------------------------------------//
double sinc(double x) {
  if (x == 0.0) {
    return 1.0; // sinc(0) = 1
  }
  return sin(x) / x;
}

//---------------------------------------------------------------------------//
//--------------------------------- Gaussian --------------------------------//
//---------------------------------------------------------------------------//
double gaussian(double x, double a, double b) {
  double aux = x / (4 * constants::pi);
  return a * exp(-1.0 * b * aux * aux);
}

//---------------------------------------------------------------------------//
//--------------------------- Alpha Weight Factor ---------------------------//
//---------------------------------------------------------------------------//
double alphaWeightFactor(double x, double a, double b) {
  double aux = x / (4 * constants::pi);
  return a * exp(-1.0 * b * aux * aux);
}

//---------------------------------------------------------------------------//
//--------------------------- Atomic Form Factor  ---------------------------//
//---------------------------------------------------------------------------//
double atomicFormFactor(double q, std::string_view element) {
  std::vector<double> param = atomicFormFactorParameters(element);
  double sum = 0.0;
  for (int i = 0; i < 4; ++i) {
    sum += gaussian(q, param[i * 2], param[i * 2 + 1]);
  }
  return sum + param[8];
}

//---------------------------------------------------------------------------//
//-------------------- Average Atomic Scattering Factor ---------------------//
//---------------------------------------------------------------------------//
double avrgAtomicScatteringFactor(std::vector<double> q_,
				     std::string_view element) {
  int n = q_.size();
  std::vector<double> aux(n);
  std::vector<double> param = atomicFormFactorParameters(element);
  double sum;
  for (int i = 0; i < n; ++i) {
    sum = 0.0;
    for (int i = 0; i < 4; ++i) {
      sum += gaussian(q_[i], param[i * 2], param[i * 2 + 1]);
    }
    aux[i] = sum + param[8];
  }
  return std::accumulate(aux.begin(), aux.end(), 0.0) / n;
}

//---------------------------------------------------------------------------//
//------------------------ Average Scattering Factor ------------------------//
//---------------------------------------------------------------------------//
double avrgScatteringFactor(std::vector<double> q_,
			      std::vector<std::string> elements) {
  int n_ = elements.size();
  double aux = 0.0;
  for (int i = 0; i < n_; ++i) {
    aux += avrgAtomicScatteringFactor(q_, elements[i]);
  }
  return aux / n_;
}
