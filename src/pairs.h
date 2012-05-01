//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _PAIRS_H
#define _PAIRS_H

#include <string>
#include <vector>
#include <iostream>

#include "ftagparams.h"

typedef std::pair<std::string, std::string> PairAlignment;
typedef std::vector<PairAlignment> Pairs;
   
std::istream& operator >> (std::istream& is, Pairs& pairs);
std::ostream& operator << (std::ostream& os, const Pairs& pairs);
std::ostream& writeFasta(std::ostream& os, const Pairs& pairs, size_t ncols);

// should use typenames / need to clean up.
FTAGParams estParams(const std::vector<std::pair<std::string, std::string> >& 
                     alignments, bool sym);

#endif
