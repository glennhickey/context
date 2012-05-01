//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _AXTREADER_H
#define _AXTREADER_H

#include <vector>
#include "axtpair.h"

class AxtReader
{

public:
   
   void read(const std::string& filename);

   const std::vector<AxtPair>& getAlignments() const {
      return _alignments;
   }

   std::vector<std::pair<std::string, std::string> > 
      sample(size_t n, size_t minLen, size_t maxLen, size_t minCutDist,
         bool onlyGapped) const;

protected:

   std::vector<AxtPair> _alignments;
};

#endif
