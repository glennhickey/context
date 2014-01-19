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

#include "ftagscore.h"
#include "ftagargs.h"
#include "pairs.h"
#include "axtreader.h"

using namespace std;
using namespace FTAGTableStructs;
using namespace FTAGStates;


int main(int argc, char** argv)
{
   try
   {
      FTAGArgs args;
      args.getArgs(argc, argv);

      // Load the model parameters
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

      // Create the model
      ContextModel em;
      TransitionModel tm;
      em.setSingleAndDouble(params);
      tm.setSimplified(params);
      FTAGModel ftag;
      ftag.setTransitionModel(tm);
      ftag.setEmissionModel(em);

      // Create the scorer
      FTAGScore tagScore(&ftag);

      // Load the AXT alignments into memory (note that this is 
      // pretty wasteful since we only need to parse one at a time, but
      // am working with the existing classes...)
      AxtReader ar;
      ar.read(args._inFile);
      const vector<AxtPair>& wholeThing = ar.getAlignments();
      
      // Score each alignment
      double totalScore = 0.0;
      for (size_t i = 0; i < wholeThing.size(); ++i)
      {
        const pair<string, string>& alignment = wholeThing[i].getAlignment();
        double score = tagScore.logViterbiScore(alignment);
        totalScore += score;
      }

      cout << totalScore << endl;
   
      return 0;
   }
   catch(const string& message)
   {
      cout << "ERROR: " << message << endl;
      return -1;
   }
}
