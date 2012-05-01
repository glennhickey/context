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

#include "contextmodel.h"
#include "ftagmodel.h"
#include "ftagargs.h"

using namespace std;
using namespace FTAGTableStructs;
using namespace FTAGStates;
pair<string, string> readSeq(const string& fname)
{
   pair<string, string> pa;
   ifstream ifile(fname.c_str());
   if (!ifile)
   {
      cerr << fname.c_str() << " not found" <<endl;
      exit(1);
   }
   ifile >> pa.first >> pa.second;
   return pa;
}

int main(int argc, char** argv)
{
   FTAGArgs args;
   args.getArgs(argc, argv);

   pair<string, string> seqs = readSeq(args._inFile);

   FTAGParams params;
   ifstream pfile(args._fpFile.c_str());
   if (!pfile)
   {
      cout << "Specify params with fpfile=<path>" << endl;
   }
   pfile >> params;

   ContextModel em;
   em.setDoubleJC(params);
   
   TransitionModel tm;
   tm.setSimplified(params);

   cout << em.getTable() << endl << tm << endl;

   FTAGModel fm;
   fm.setEmissionModel(em);
   fm.setTransitionModel(tm);
   fm.setSequences(seqs.first, seqs.second, args._winSize, args._winSize);
   cout << "fp(A,B) = " << fm.forward() << endl;
   cout << "bp(A,B) = " << fm.backward() << endl;
   cout << "v(A,B) = " << fm.viterbi() << endl;

   string a,b;
   deque<Trace> trace;
   fm.viterbiTrace(a, b, trace);
   vector<double> sdist;
   fm.computeStateDistribution(sdist);
   cout << a << endl << b << endl;

//    copy(trace.begin(), trace.end(), ostream_iterator<Trace>(cout, "\n"));
   for (deque<Trace>::iterator i = trace.begin(); i != trace.end(); ++i)
      cout << *i << endl;
   for (size_t s = 0; s < S_Max; ++s)
   {
      cout << (State)s << "=" << sdist[s] << " ";
   }
      

   return 0;
}
