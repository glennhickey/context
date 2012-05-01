//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>

#include "ftagtrain.h"
#include "ftagargs.h"
#include "pairs.h"

using namespace std;
using namespace FTAGTableStructs;
using namespace FTAGStates;


int main(int argc, char** argv)
{
   try{
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
   pfile >> params;
   params.setSymmetric(args._symmetric);
   params.setSingleF84(args._singleF84);
   params.setDoubleF84(args._doubleF84);

   ifstream ifile(args._inFile.c_str());
   if (!pfile)
   {
      cerr << "Specify input sequences with infile=<path>" << endl;
      exit(1);
   }

   ofstream ofile(args._outFile.c_str());
   if (!ofile)
   {
      cerr << "Specify output file with outfile=<path>" << endl;
      exit(1);
   }

   Pairs inputPairs, outputPairs;
   ifile >> inputPairs;

   cout << "input file contains " << inputPairs.size() << " pairs" << endl;
   
   ContextModel em;
   TransitionModel tm;

   em.setSingleAndDouble(params);
   tm.setSimplified(params);

   double totalFll = 0.;
   
   for (size_t i = 0; i < inputPairs.size(); ++i)
   {
      FTAGModel ftag;
      ftag.setSequences(inputPairs[i].first, inputPairs[i].second,
                        args._winSize, args._winSize);
      
      ftag.setTransitionModel(tm);
      ftag.setEmissionModel(em);
      double fll = log(ftag.forward());
      
      totalFll += fll;
      
      ofile << i + 1 << " " << fll << endl;
   }
   ofile << "total " << totalFll << endl;
   cout << "Total Forward Probability = " << totalFll << endl;
   }
   catch (string message)
   {
      cerr << message << endl;
   }
   return 0;
}
