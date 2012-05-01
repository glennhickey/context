//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _AXTPAIR_H
#define _AXTPAIR_H

#include <vector>
#include <string>
#include <istream>

#include "dnasubmodel.h"

class AxtPair
{
   friend std::istream& operator>>(std::istream& is, AxtPair& ap);

public:

   struct GapInfo 
   {
      size_t _distLeft;  // distance to nearest gap to the left
      size_t _lenLeft;   // length of that gap
      size_t _distRight; // distance to nearest gap to the right
      size_t _lenRight;  // length of that gap
   };


   AxtPair();
   AxtPair(const std::string& seqA, const std::string& seqB);
   
   const std::pair<std::string, std::string>& getAlignment() const {
      return _alignment;
   }

   size_t getIdx() const {
      return _idx;
   }

   const std::string& getNameA() const {
      return _nameA;
   }

   const std::string& getNameB() const {
      return _nameB;
   }

   const size_t getStartA() const { 
      return _startA;
   }
   
   const size_t getStartB() const {
      return _startB;
   }
   
   const size_t getEndA() const { 
      return _endA;
   }
   
   const size_t getEndB() const {
      return _endB;
   }

   const size_t getLength() const {
      assert(_alignment.first.length() == _alignment.second.length());
      return _alignment.first.length();
   }

   
   bool before(const AxtPair& ap) const;
   bool after(const AxtPair& ap) const;

   AxtPair operator+(const AxtPair& ap) const;

   bool hasN() const;
   bool hasGap() const;
   
   std::vector<std::pair<std::string, std::string> > 
      chop(size_t minLen, size_t maxLen, size_t minCutDist = 5) const;

   int diff(const AxtPair& ap) const;

   void annotate(const std::pair<std::string, std::string> input, 
                 std::vector<GapInfo>& gapInfo) const;

protected:

   static const size_t ReadBufferSize;
   static char* ReadBuffer;
   
   size_t _idx;
   std::string _nameA;
   std::string _nameB;
   size_t _startA;
   size_t _startB;
   size_t _endA;
   size_t _endB;
   
   std::pair<std::string, std::string> _alignment;

};

std::istream& operator>>(std::istream& is, AxtPair& ap);

#endif
