//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <string>
#include <cstdlib>
#include <ostream>
#include <iterator>
#include <fstream>

#include "ftaggen.h"
#include "ftagargs.h"

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
      cout << "Specify params with fpfile=<path>" << endl;
   }
   pfile >> params;

   ContextModel em;
   em.setDoubleJC(params);
   TransitionModel tm;
   tm.setSimplified(params);

   FTAGGen sGen;
   sGen.setEmissionModel(em);
   sGen.setTransitionModel(tm);
      
   string a,b;
   deque<Trace> trace;
   sGen.genAlignment(a, b, trace, true, 10000);
   
   if (drand48() < args._flipProb)
   {
      reverse(a.begin(), a.end());
      reverse(b.begin(), b.end());
   }
   
   cout << a << endl << b << endl;
//   copy(trace.begin(), trace.end(), ostream_iterator<Trace>(cout, "\n"));
   for (deque<Trace>::iterator i = trace.begin(); i != trace.end(); ++i)
      cout << *i << endl;
      
   cout << trace.size() << " state transitions used" << endl;
   
   return 0;
}
