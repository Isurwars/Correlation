#ifndef SRC_STRUCTURE_FACTOR_H_
#define SRC_STRUCTURE_FACTOR_H_
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
#include <cmath>
#include <string>
#include <vector>

//---------------------------------------------------------------------------//
//------------------------------- Methods -----------------------------------//
//---------------------------------------------------------------------------//

// Sinc Function, sin(x)/x
double sinc(double);
// Gaussian Function a*exp (-b * (x / 4Pi)**2)
double gaussian(double, double, double);
// Atomic Form Factor
double atomicFormFactor(double, std::string_view);
// Average Scattering Factor
double avrgAtomicScatteringFactor(std::vector<double>, std::string_view);
double avrgScatteringFactor(std::vector<double>, std::vector<std::string>);
#endif // SRC_STRUCTURE_FACTOR_H_
