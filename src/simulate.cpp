//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <string>
#include <cstdlib>
#include <ostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include "ftaggen.h"
#include "ftagargs.h"
#include "pairs.h"

using namespace std;
using namespace FTAGStates;
using namespace FTAGTableStructs;

int main(int argc, char** argv)
{
   FTAGArgs args;
   args.getArgs(argc, argv);
   srand48(args._seed);
   cout << "SEED=" << args._seed << endl;

   FTAGParams params;
   ifstream pfile(args._fpFile.c_str());
   if (!pfile)
   {
      cerr << "Specify params with fpfile=<path>" << endl;
      exit(1);
   }
   ofstream outParFile(args._fpFileOut.c_str());
   if (args._fpFileOut.length() > 0 && ! outParFile.is_open())
   {
      cerr << "error opening fpoout" << endl;
      exit(1);
   }

   pfile >> params;
   params.setSymmetric(args._symmetric);
   params.setSingleF84(args._singleF84);
   params.setDoubleF84(args._doubleF84);
   
   ofstream ofile(args._outFile.c_str());
   if (!ofile)
   {
      cerr << "Specifiy output file with ofile=<path>" << endl;
      exit(1);
   }

   ContextModel em;
   em.setSingleAndDouble(params);
   TransitionModel tm;
   tm.setSimplified(params);

   FTAGGen sGen;
   sGen.setEmissionModel(em);
   sGen.setTransitionModel(tm);

   Pairs pairs;
   for (int i = 0; i < (int)args._numPairs; ++i)
   {
      PairAlignment pa;
      deque<Trace> trace;

      sGen.genAlignment(pa.first, pa.second, trace, args._keepGaps, args._maxLength);
      if (drand48() < args._flipProb)
      {
         reverse(pa.first.begin(), pa.first.end());
         reverse(pa.second.begin(), pa.second.end());
      }
      if (pa.first.length() >= args._winSize && 
          pa.second.length() >= args._winSize &&
          pa.first.length() >= args._minLength &&
          pa.second.length() >= args._minLength)
      {
         pairs.push_back(pa);
      }
      else
      {
         --i;
      }
   }

   if (args._simFasta)
      writeFasta(ofile, pairs, 80);
   else
      ofile << pairs;

   if (outParFile.is_open())
   {
      FTAGParams estPar = estParams(pairs, args._symmetric);
      outParFile << estPar;
   }
   
   return 0;
}
