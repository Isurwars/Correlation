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
#include <fstream>
#include <iomanip>
#include <iostream>

#include "WriteFiles.hpp"

//---------------------------------------------------------------------------//
//------------------------------- Templates ---------------------------------//
//---------------------------------------------------------------------------//

template <typename T>
void write2DVectorToCSV(const std::vector<std::vector<T>> &vec,
			const std::string &filename,
			const std::string &header = "") {
  std::ofstream file(filename);
  int rows = vec[0].size();
  int cols = vec.size();
  if (!file.is_open()) {
    std::cerr << "Error: Could not open the file for writing." << std::endl;
    return;
  }
  file << std::fixed << std::setprecision(5);
  file << header << std::endl;
  // Set precision for floating-point numbers
  for (int j = 0; j < rows; ++j) {
    for (int i = 0; i < cols; ++i) {
      file << std::setw(12) << std::right << vec[i][j];
      if (i < cols - 1) {
	file << ",";
      }
    }
    file << std::endl;
  }

  file.close();
  file.clear();
  if (!file) {
    std::cerr << "Error: Could not close the file properly." << std::endl;
  }
}

//---------------------------------------------------------------------------//
//-------------------------------- Methods ----------------------------------//
//---------------------------------------------------------------------------//

// Write CSV Histograms (Cell, out_file, Smoothing)
void WriteCSV(Cell cell, std::string filename, bool smoothing) {
  std::stringstream header;
  std::vector<std::vector<double>> temp_hist;
  std::vector<std::vector<int>> temp_hist_int;
  int i, j, k, n = cell.distancesSize();

  //-------------------------------------------------------------------------//
  //------------------------------- Writing J -------------------------------//
  //-------------------------------------------------------------------------//
  temp_hist = cell.J();
  header << std::right << std::setw(14) << "r (Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header << std::right << std::setw(14)
	     << cell.elements()[i] + "-" + cell.elements()[j] + " (1/Å),";
    }
  }
  header << std::right << std::setw(12) << "J(r)";
  write2DVectorToCSV(temp_hist, filename + "_J.csv", header.str());
  header.str("");
  header.clear();
  std::cout << "Writing J(r)" << std::endl;

  //-------------------------------------------------------------------------//
  //------------------------------- Writing g -------------------------------//
  //-------------------------------------------------------------------------//
  temp_hist = cell.g();
  header << std::right << std::setw(14) << "r (Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header << std::right << std::setw(14)
	     << cell.elements()[i] + "-" + cell.elements()[j] + " (1/Å),";
    }
  }
  header << std::right << std::setw(12) << "g(r)";
  write2DVectorToCSV(temp_hist, filename + "_g.csv", header.str());
  header.str("");
  header.clear();
  std::cout << "Writing g(r)" << std::endl;

  //-------------------------------------------------------------------------//
  //------------------------------- Writing G -------------------------------//
  //-------------------------------------------------------------------------//
  temp_hist = cell.G();
  header << std::right << std::setw(14) << "r (Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header << std::right << std::setw(14)
	     << cell.elements()[i] + "-" + cell.elements()[j] + " (1/Å),";
    }
  }
  header << std::right << std::setw(12) << "G(r)";
  write2DVectorToCSV(temp_hist, filename + "_G_.csv", header.str());
  header.str("");
  header.clear();
  std::cout << "Writing G(r)" << std::endl;

  //-------------------------------------------------------------------------//
  //------------------------------ Writing PAD ------------------------------//
  //-------------------------------------------------------------------------//
  temp_hist = cell.F();
  header << std::right << std::setw(14) << "theta (°),";
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      for (k = j; k < n; k++) {
	header << std::right << std::setw(13)
	       << cell.elements()[j] + "-" + cell.elements()[i] + "-" +
		      cell.elements()[k] + ",";
      }
    }
  }
  header << std::right << std::setw(12) << "F(theta)";
  write2DVectorToCSV(temp_hist, filename + "_PAD.csv", header.str());
  header.str("");
  header.clear();
  std::cout << "Writing F(theta)" << std::endl;

  //-------------------------------------------------------------------------//
  //------------------------------- Writing Z -------------------------------//
  //-------------------------------------------------------------------------//
  temp_hist_int = cell.Z();
  header << std::right << std::setw(13) << "Number (#),";
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      header << std::right << std::setw(13)
	     << cell.elements()[i] + " by " + cell.elements()[j] + ",";
    }
  }
  for (i = 0; i < n; i++) {
    header << std::right << std::setw(12) << cell.elements()[i] + "  by any";
    if (i < n - 1) {
      header << ",";
    }
  }
  write2DVectorToCSV(temp_hist_int, filename + "_Z.csv", header.str());
  header.str("");
  header.clear();
  std::cout << "Writing Coordination Number (Z)" << std::endl;

  //-------------------------------------------------------------------------//
  //------------------------------ Writing SQ -------------------------------//
  //-------------------------------------------------------------------------//
  temp_hist = cell.S();
  header << std::right << std::setw(13) << "q (1/Å),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header << std::right << std::setw(12)
	     << cell.elements()[i] + "-" + cell.elements()[j] + " (Å),";
    }
  }
  header << std::right << std::setw(12) << "S(Q)";
  write2DVectorToCSV(temp_hist, filename + "_S.csv", header.str());
  header.str("");
  header.clear();
  std::cout << "Writing S(Q)" << std::endl;

  //-------------------------------------------------------------------------//
  //------------------------------ Writing XRD ------------------------------//
  //-------------------------------------------------------------------------//
  temp_hist = cell.X();
  header << std::right << std::setw(13) << "2-theta (°),";
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      header << std::right << std::setw(12)
	     << cell.elements()[i] + "-" + cell.elements()[j] + " (%),";
    }
  }
  header << std::right << std::setw(12) << "XRD";
  write2DVectorToCSV(temp_hist, filename + "_XRD.csv", header.str());
  header.str("");
  header.clear();
  std::cout << "Writing XRD(2theta)" << std::endl;

  if (smoothing) {
    //-----------------------------------------------------------------------//
    //-------------------------- Writing smoothed J -------------------------//
    //-----------------------------------------------------------------------//
    temp_hist = cell.J_smoothed();
    header << std::right << std::setw(14) << "r (Å),";
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
	header << std::right << std::setw(14)
	       << cell.elements()[i] + "-" + cell.elements()[j] + " (1/Å), ";
      }
    }
    header << std::right << std::setw(12) << "J(r)";
    write2DVectorToCSV(temp_hist, filename + "_J_smoothed.csv", header.str());
    header.str("");
    header.clear();
    std::cout << "Writing J(r) smoothed" << std::endl;

    //-----------------------------------------------------------------------//
    //-------------------------- Writing smoothed g -------------------------//
    //-----------------------------------------------------------------------//
    temp_hist = cell.g_smoothed();
    header << std::right << std::setw(14) << "r (Å),";
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
	header << std::right << std::setw(14)
	       << cell.elements()[i] + "-" + cell.elements()[j] + " (1/Å), ";
      }
    }
    header << std::right << std::setw(12) << "g(r)";
    write2DVectorToCSV(temp_hist, filename + "_g_smoothed.csv", header.str());
    header.str("");
    header.clear();
    std::cout << "Writing g(r) smoothed" << std::endl;

    //-----------------------------------------------------------------------//
    //-------------------------- Writing smoothed G -------------------------//
    //-----------------------------------------------------------------------//
    temp_hist = cell.G_smoothed();
    header << std::right << std::setw(14) << "r (Å),";
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
	header << std::right << std::setw(14)
	       << cell.elements()[i] + "-" + cell.elements()[j] + " (1/Å), ";
      }
    }
    header << std::right << std::setw(12) << "G(r)";
    write2DVectorToCSV(temp_hist, filename + "_G_smoothed_.csv", header.str());
    header.str("");
    header.clear();
    std::cout << "Writing G(r) smoothed" << std::endl;

    //-----------------------------------------------------------------------//
    //------------------------ Writing smoothed PAD -------------------------//
    //-----------------------------------------------------------------------//
    temp_hist = cell.F_smoothed();
    header << std::right << std::setw(14) << "theta (°),";
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	for (k = j; k < n; k++) {
	  header << std::right << std::setw(13)
		 << cell.elements()[j] + "-" + cell.elements()[i] + "-" +
			cell.elements()[k] + ",";
	}
      }
    }
    header << std::right << std::setw(12) << "F(theta)";
    write2DVectorToCSV(temp_hist, filename + "_PAD_smoothed.csv", header.str());
    header.str("");
    header.clear();
    std::cout << "Writing F(theta) smoothed" << std::endl;

    //-----------------------------------------------------------------------//
    //------------------------ Writing smoothed S(q) ------------------------//
    //-----------------------------------------------------------------------//
    temp_hist = cell.S_smoothed();
    header << std::right << std::setw(13) << "q (1/Å),";
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
	header << std::right << std::setw(12)
	       << cell.elements()[i] + "-" + cell.elements()[j] + " (Å),";
      }
    }
    header << std::right << std::setw(12) << "S(Q)";
    write2DVectorToCSV(temp_hist, filename + "_S_smoothed.csv", header.str());
    header.str("");
    header.clear();
    std::cout << "Writing S(Q) smoothed" << std::endl;

    //-----------------------------------------------------------------------//
    //------------------------ Writing smoothed XRD -------------------------//
    //-----------------------------------------------------------------------//
    temp_hist = cell.X();
    header << std::right << std::setw(13) << "2-theta (°),";
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
	header << std::right << std::setw(12)
	       << cell.elements()[i] + "-" + cell.elements()[j] + " (%),";
      }
    }
    header << std::right << std::setw(12) << "XRD";
    write2DVectorToCSV(temp_hist, filename + "_XRD_smoothed.csv", header.str());
    header.str("");
    header.clear();
    std::cout << "Writing XRD(2 theta) smoothed" << std::endl;
  }
};
