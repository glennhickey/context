//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <fstream>
#include <iostream>
#include <algorithm>
#include <locale>
#include "axtreader.h"
#include "ftagargs.h"
#include "pairs.h"
#include "ftagparams.h"
#include "dnasubmodel.h"

using namespace std;


int main(int argc, char** argv)
{
   FTAGArgs args;
   args.getArgs(argc, argv);
   srand48(args._seed);
   cout << "SEED=" << args._seed << endl;

   ofstream ofile(args._outFile.c_str());
   if (!ofile)
   {
      cerr << "Specifiy output file with ofile=<path>" << endl;
      exit(1);
   }

   ofstream pfile(args._fpFileOut.c_str());

   AxtReader ar;

   try
   {
      ar.read(args._inFile);
      
      const vector<AxtPair>& wholeThing = ar.getAlignments();
      Pairs allAlignments;
      for (size_t i = 0; i < wholeThing.size(); ++i)
      {
         allAlignments.push_back(wholeThing[i].getAlignment());
      }
      Pairs alignments = ar.sample(args._numPairs, args._maxLength / 2, 
                                   args._maxLength, 0, false);

      FTAGParams parAll = estParams(allAlignments, args._symmetric);
      FTAGParams parInput = estParams(alignments, args._symmetric);

      cout << "From Whole File: " << parAll << endl;
      cout << "From input: " << parInput << endl;

      ofile << alignments;

      if (pfile.is_open())
      {
         pfile << parInput;
      }
   }
   catch(string message)
   {
      cerr << message << endl;
      return 1;
   }

   return 0;
}


