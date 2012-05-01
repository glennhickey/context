//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _PAIRALIGNMENT_H
#define _PAIRALIGNMENT_H

#include <deque>
#include <ostream>

struct PairAlignment
{
   std::deque<char> _seqA;
   std::deque<char> _seqB;
};

inline std::ostream& operator<<(std::ostream& os, const PairAlignment& pa)
{
   size_t i;
   for (i = 0; i < pa._seqA.size(); ++i)
   {
      os << pa._seqA[i];
   }
   os << std::endl;
   for (i = 0; i < pa._seqB.size(); ++i)
   {
      os << pa._seqB[i];
   }
   return os;
}
#endif
