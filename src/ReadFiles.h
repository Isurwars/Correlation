#ifndef SRC_READFILES_H_
#define SRC_READFILES_H_

/*
 * Copyright [2021] <@isurwars>
 */
#include <string>

#include "Cell.h"

Cell read_CAR(std::string);
Cell read_CELL(std::string);
Cell read_ONETEP_DAT(std::string);
Cell read_CIF(std::string);

#endif  // SRC_READFILES_H_
