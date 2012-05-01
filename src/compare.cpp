//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <fstream>
#include <iostream>
#include <algorithm>
#include <locale>
#include "axtpair.h"
#include "ftagargs.h"
#include "pairs.h"
#include "ftagparams.h"

using namespace std;

int main(int argc, char** argv)
{
   FTAGArgs args;
   args.getArgs(argc, argv);
   srand48(args._seed);

   ifstream ifile(args._inFile.c_str());
   if (!ifile)
   {
      cerr << "Specifiy input file with infile=<path>" << endl;
      exit(1);
   }

   ifstream ifile2(args._inFile2.c_str());
   if (!ifile2)
   {
      cerr << "Specifiy input file2 with infile2=<path>" << endl;
      exit(1);
   }

   ofstream ofile(args._outFile.c_str());

   Pairs setOne;
   Pairs setTwo;
   
   ifile >> setOne;
   ifile2 >> setTwo;

   if (setOne.size() != setTwo.size() || setOne.size() == 0)
   {
      cerr << "Files contain diff. num of alignments (" << 
         setOne.size() << " vs " << setTwo.size() << endl;
      exit(1);
   }

   size_t diffTotal = 0.;
   size_t len1Total = 0;
   size_t len2Total = 0;
   for (size_t i = 0; i < setOne.size(); ++i)
   {
      AxtPair ap1(setOne[i].first, setOne[i].second);
      AxtPair ap2(setTwo[i].first, setTwo[i].second);

      int diff = ap1.diff(ap2);
      diffTotal += abs(diff);
      len1Total += ap1.getLength();
      len2Total += ap2.getLength();
      if (ofile)
      {
         ofile << i + 1 << ") diff: " << diff 
               << " len1: " << ap1.getLength()
               << " len2: " << ap2.getLength() << endl;
      }
   }

   if (ofile)
   {
      ofile << "total) diff: " << diffTotal
            << " len1: " << len1Total
            << " len2: " << len2Total << endl;
      
   }

   cout << "total) diff: " << diffTotal
        << " len1: " << len1Total
        << " len2: " << len2Total << endl;

   return 0;
}
