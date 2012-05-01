//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <fstream>
#include <iostream>
#include <cstdlib>

#include "axtreader.h"

using namespace std;

void AxtReader::read(const string& filename)
{
   ifstream inFile(filename.c_str());
   
   while (!inFile.eof())
   {
      AxtPair ap;
      inFile >> ap;
      _alignments.push_back(ap);
   }
}

// quite inefficeint as strings get copied several times thourhgout 
// this method and chop, but should not be a big enough deal to warrant
// improving.
vector<pair<string, string> > 
AxtReader::sample(size_t n, size_t minLen, size_t maxLen, 
                  size_t minCutDist, bool onlyGapped) const
{
   vector<pair<string, string> > chopped;
   for (size_t i = 0; i < _alignments.size(); ++i)
   {
      vector<pair<string, string> > temp = 
         _alignments[i].chop(minLen, maxLen, minCutDist);
      
      for (size_t j = 0; j < temp.size(); ++j)
      {
         if (onlyGapped == false ||
             temp[j].first.find_first_of("-") != string::npos ||
             temp[j].second.find_first_of("-") != string::npos)
         {
            chopped.push_back(temp[j]);
         }
      }
   }
   
   double p = (double)n / (double)chopped.size();

   vector<pair<string, string> > sampled;
   for (size_t i = 0; i < chopped.size(); ++i)
   {
      if (drand48() <= p)
      {
         sampled.push_back(chopped[i]);
      }
      if (i == chopped.size() - 1 && sampled.size() < n)
      {
         i = 0;
      }
      else if (sampled.size() == n)
      {
         break;
      }
   }
   return sampled;
}
