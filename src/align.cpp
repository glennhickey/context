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

static void writePosterior(size_t num, const FTAGModel& ftag, 
									const PairAlignment& alignment,
                           const deque<Trace>& trace, ostream& postfile);

int main(int argc, char** argv)
{
   try
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
      pfile >> params;
      params.setSymmetric(args._symmetric);
      params.setSingleF84(args._singleF84);
      params.setDoubleF84(args._doubleF84);

      ifstream ifile(args._inFile.c_str());
      if (!ifile)
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

      ofstream postfile;
      if (args._outFile2.empty() == false)
      {
         postfile.open(args._outFile2.c_str());
      }

      Pairs inputPairs, outputPairs;
      ifile >> inputPairs;

      cout << "input file contains " << inputPairs.size() << " pairs" << endl;
   
      ContextModel em;
      TransitionModel tm;
      em.setSingleAndDouble(params);
      tm.setSimplified(params);

      double totalVll = 0.;
   
      for (size_t i = 0; i < inputPairs.size(); ++i)
      {
         FTAGModel ftag;
         ftag.setSequences(inputPairs[i].first, inputPairs[i].second,
                           args._winSize, args._winSize);

         ftag.setTransitionModel(tm);
         ftag.setEmissionModel(em);
         double vll = log(ftag.viterbi());

         deque<Trace> trace;
         string a, b;
         ftag.viterbiTrace(a, b, trace);

         outputPairs.push_back(PairAlignment(a, b));

         if (postfile.is_open())
         {
            ftag.forward();
            ftag.backward();
            writePosterior(i, ftag, outputPairs.back(), trace, postfile);
         }

         totalVll += vll;
      }

      ofile << outputPairs;
      cout << "Total Viterbi Probability = " << totalVll << endl;
   
      return 0;
   }
   catch(const string& message)
   {
      cout << "ERROR: " << message << endl;
      return -1;
   }
}

void writePosterior(size_t num, const FTAGModel& ftag, 
                    const PairAlignment& alignment,
                    const deque<Trace>& trace, ostream& postfile)
{
   size_t pos = 0;
   vector<pair<double, Trace> > posts(alignment.first.length(), 
                                      pair<double, Trace>(0., Trace()));
   double lprob = 0;
   size_t i, j;
   State s;

   postfile << num << ")" << endl;

   for (i = 0; i < trace.size(); ++i)
   {
      s = trace[i]._s;
      if (s == S_Mx || s == S_Dx || s == S_Ix)
      {
         lprob = log(ftag.posterior(trace[i]));
         posts[pos].first = lprob;
         posts[pos].second = trace[i];
         ++pos;
      }
      else if (s == S_MDx || s == S_MIx)
      {
         size_t len = s == S_MDx ? trace[i]._x3 - trace[i]._x1 :
            trace[i]._x4 - trace[i]._x2;
         assert ((s == S_MDx && trace[i]._x2 == trace[i]._x3) ||
                 (s == S_MIx && trace[i]._x3 == trace[i]._x4));

         for (j = 0; j < len; ++j)
         {
            lprob = log(ftag.posterior(trace[i+j]));
            posts[pos + j].first = lprob;
            posts[pos + j + len].first = lprob;
            posts[pos + j].second = trace[i+j];
            posts[pos + j + len].second = trace[i+j];
         } 
         pos = pos + 2 * len;
         i = i + len; // also cover close state
      }
   }
   
   for (size_t i = 0; i < posts.size(); ++i)
   {
      postfile << i << " " << posts[i].first 
					<< " " << posts[i].second << endl;
   }
   
   postfile << endl;
}
