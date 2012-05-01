//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cctype>
#include "axtpair.h"

using namespace std;

AxtPair::AxtPair() : _idx(0), _startA(0), _startB(0), _endA(0), _endB(0)
{

}

AxtPair::AxtPair(const string& seqA, const string& seqB)
   : _idx(0), _startA(0), _startB(0), _endA(0), _endB(0),
     _alignment(pair<string, string>(seqA, seqB))
{

}

istream& operator>>(istream& is, AxtPair& ap)
{
   if (!is) 
   {
      assert(false);
      throw string("File Read Error 1");
   }

   string buffer;
   char lineBuf[512];

   while (is.peek() == '#')
   {
      is.getline(lineBuf, 512);
   }

   if (!is) 
   {
      assert(false);
      throw string("File Read Error 2");
   }
  
   is >> ap._idx
      >> ap._nameA >> ap._startA
      >> ap._endA >> ap._nameB >> ap._startB >> ap._endB;
   
   if (!is) 
   {
      assert(false);
      throw string("File Read Error 3");
   }

   is >> buffer >> buffer;

   if (!is) 
   {
      assert(false);
      throw string("File Read Error 4");
   }

   is >> ap._alignment.first >> ap._alignment.second;

   if (!is) 
   {
      assert(false);
      throw string("File Read Error 5");
   }
   
   // EAT UP TRAILING WHITESPACE
   while(is && !isalnum(is.peek()))
   {
      is.get();
   }
   
   return is;
}

bool AxtPair::before(const AxtPair& ap) const
{
   return 
      ap._nameA == _nameA && 
      ap._nameB == _nameB &&
      ap._idx == _idx + 1 &&
      ap._startA == _endA + 1 &&
      ap._startB == _endB + 1;
}

bool AxtPair::after(const AxtPair& ap) const
{
   return ap.before(*this);
}

AxtPair AxtPair::operator+(const AxtPair& ap) const
{
   AxtPair merged;
   if (after(ap) == true)
   {
      return ap + *this;
   }
   else if (before(ap) == true)
   {
      merged = *this; // default assignment operator
      merged._alignment.first += ap._alignment.first;
      merged._alignment.second += ap._alignment.second;
      merged._endA = ap._endA;
      merged._endB = ap._endB;
   }
   else
   {
      assert(false);
   }
   return merged;
}

bool AxtPair::hasN() const
{
   return _alignment.first.find_first_of("Nn") != string::npos ||
      _alignment.second.find_first_of("Nn") != string::npos;
}

bool AxtPair::hasGap() const
{
   return _alignment.first.find_first_of("-") != string::npos ||
      _alignment.second.find_first_of("-") != string::npos;
}

vector<pair<string, string> > 
AxtPair::chop(size_t minLen, size_t maxLen, size_t minCutDist) const
{
   size_t i;

   // STRIP THE N's
   pair<string, string> nStripped;
   for (i = 0; i < getLength(); ++i)
   {
      if (_alignment.first[i] != 'n' && _alignment.second[i] != 'n' &&
          _alignment.first[i] != 'N' && _alignment.second[i] != 'N')
      {
         nStripped.first += _alignment.first[i];
         nStripped.second += _alignment.second[i];
      }
   }
   assert(nStripped.first.length() == nStripped.second.length());
   size_t n = nStripped.first.length();

   // COMPUTE DISTANCES TO NEAREST GAP
   vector<GapInfo> gapInfo;
   annotate(nStripped, gapInfo);

   // CHOP UP SEQUENCE TO BITS BETWEEN MIN AND MAX LENGTHS, SUCH THAT
   // NO CUT HAPPENS WITHIN minCutDist OF A GAP
   int curStartPos = -1;
   pair<string, string> seqPair;
   vector<pair<string, string> > seqPairs;

   for (i = 0; i < n; ++i)
   {
      // start new sequence
      if (curStartPos == -1 && 
          gapInfo[i]._distRight >= minCutDist + gapInfo[i]._lenRight)
      {
         curStartPos = i;
         seqPair.first = nStripped.first[i];
         seqPair.second = nStripped.second[i];
      }
      // end sequence
      else if ((i - curStartPos == maxLen || 
               (i == n - 1 && i - curStartPos >= minLen)) &&
               gapInfo[i]._distLeft >= minCutDist + gapInfo[i]._lenLeft)
      {
         seqPair.first += nStripped.first[i];
         seqPair.second += nStripped.second[i];
         seqPairs.push_back(seqPair);
         curStartPos = -1;
      }
      // give up on current sequence
      else if (i - curStartPos > maxLen)
      {
         curStartPos = -1;
      }
      // continute current sequence
      else if (curStartPos >= 0)
      {
         seqPair.first += nStripped.first[i];
         seqPair.second += nStripped.second[i];
      }
   }
     
   return seqPairs;
}

int AxtPair::diff(const AxtPair& ap) const
{
   vector<pair<size_t, size_t> > o1, o2;
   pair<size_t, size_t> count1(0, 0), count2(0, 0);
   for (size_t i = 0; i < getLength(); ++i)
   {
      if (_alignment.first[i] != '-' && _alignment.second[i] != '-')
      {
         o1.push_back(count1);
      }
      if (_alignment.first[i] != '-')
      {
         ++count1.first;
      }
      if (_alignment.second[i] != '-')
      {
         ++count1.second;
      }
   }

   for (size_t i = 0; i < ap.getLength(); ++i)
   {
      if (ap._alignment.first[i] != '-' && ap._alignment.second[i] != '-')
      {
         o2.push_back(count2);
      }
      if (ap._alignment.first[i] != '-')
      {
         ++count2.first;
      }
      if (ap._alignment.second[i] != '-')
      {
         ++count2.second;
      }
   }

   assert(count1 == count2);
   
   vector<pair<size_t, size_t> > intersection(max(o1.size(), o2.size()));
   vector<pair<size_t, size_t> >::iterator last = 
      set_intersection(o1.begin(), o1.end(), o2.begin(), o2.end(), 
                       intersection.begin());

   return o1.size() - (last - intersection.begin());
}

void AxtPair::annotate(const pair<string, string> input, 
                       vector<GapInfo>& gapInfo) const
{
   assert(input.first.length() == input.second.length());
   size_t n = input.first.length();
   gapInfo.resize(n);
   
   size_t curLen = 0;
   size_t curDist = n;
   bool inGap = false;

   // forward
   
   for (size_t i = 0; i < n; ++i)
   {
      if (input.first[i] == '-' || input.second[i] == '-')
      {
         if (inGap == false)
         {
            inGap = true;
            curLen = 0;
         }
         ++curLen;
         curDist = 0;
      }
      else
      {
         inGap = false;
         ++curDist;
      }
      gapInfo[i]._distLeft = curDist;
      gapInfo[i]._lenLeft = curLen;
   }

   // backward (same code different direction)

   curLen = 0;
   curDist = n;
   inGap = false;
   
   for (int i = n-1; i>= 0; --i)
   {
      if (input.first[i] == '-' || input.second[i] == '-')
      {
         if (inGap == false)
         {
            inGap = true;
            curLen = 0;
         }
         ++curLen;
         curDist = 0;
      }
      else
      {
         inGap = false;
         ++curDist;
      }
      gapInfo[i]._distRight = curDist;
      gapInfo[i]._lenRight = curLen;
   }
}
