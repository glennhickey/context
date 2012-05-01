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
#include "ftaggen.h"
#include "ftagargs.h"
#include "axtpair.h"

using namespace std;
using namespace FTAGTableStructs;
using namespace FTAGStates;

static void makeEmGraph(const FTAGTrain::EMTrace& trace, 
                        double prOrig, const FTAGParams& parOrig,
                        size_t numSites,
                        ostream& os);

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
   
   vector<pair<string, string> > inputPairs;
   vector<pair<size_t, size_t> > inputWins;
   vector<deque<Trace> > traces;
   pair<string, string> seqPair;
   pair<size_t, size_t> winPair;
   deque<Trace> trace;
    
   size_t extras = 0;
   for (size_t i = 0; i < args._numPairs + extras; ++i)
   {
      sGen.genAlignment(seqPair.first, seqPair.second, trace, true, 
                        args._maxLength);

      // TODO Test on empty.small sequences
      // (but for now we skip)
      if (seqPair.first.length() >= args._minLength &&
          seqPair.second.length() >= args._minLength)
      {
         if (drand48() < args._flipProb)
         {
            reverse(seqPair.first.begin(), seqPair.first.end());
            reverse(seqPair.second.begin(), seqPair.second.end());
         }

         inputPairs.push_back(seqPair);

         winPair.first = min(seqPair.first.length(), args._winSize);
         winPair.second = min(seqPair.second.length(), args._winSize);
      
         inputWins.push_back(winPair);

         cout << i << "A) " << seqPair.first << endl
              << i << "B) " << seqPair.second << endl;
      }
      else
      {
         ++extras;
      }
   }
   
   FTAGParams parSeed;
   parSeed = params; // copy fixed flags
//   parSeed.setClamp(FTAGParams::Par_Gamma, 0, 0.05);
   if (args._context == false)
   {
      parSeed.setRMD(0.);
      parSeed.setRMI(0.);
      parSeed.setPCD(1.);
      parSeed.setPCI(1.);
      parSeed.setRMDFixed(true);
      parSeed.setRMIFixed(true);
      parSeed.setPCDFixed(true);
      parSeed.setPCIFixed(true);
   }
   parSeed.setSymmetric(args._symmetric);
   parSeed.randomizeNonFixed();


   ContextModel emSeed;
   TransitionModel tmSeed;
   emSeed.setDoubleJC(parSeed);
   tmSeed.setSimplified(parSeed);

   FTAGParams offset;
   offset.setAll(args._offset);
   AGDOptimizer optimizer;
   optimizer.setFixed(params);
   optimizer.setRelative(args._relative);
   optimizer.setSavedInit(parSeed);
   optimizer.setLoopParams(args._numOptTrials, args._maxOptIt, 
                           args._optThreshold, offset); 

   FTAGTrain* trainer = new FTAGTrain();
   EmissionEstimator ee;
   TransitionEstimator te;
   ee.setBias(args._emBias);
   te.setBias(args._tmBias);
   trainer->initialize(optimizer, ee, te);
   trainer->setSequences(inputPairs, inputWins);

   ofstream logFile;
   FTAGTrain::EMTrace emTrace;
   if (!args._outFile.empty())
   {
      logFile.open(args._outFile.c_str());
   }

   cout << "Origina param\n" << params << endl;
   cout << "Seed param\n" << parSeed << endl;
   trainer->emLoop(args._maxEMIt, args._emThreshold, args._emConvRepeats,
                   args._emCRThreshold, &emTrace);
      
   FTAGParams parOut = trainer->getParams();
   ContextModel emOut = trainer->getEmissionModel();
   TransitionModel tmOut = trainer->getTransitionModel(); 

   emOut.setDoubleJC(parOut);
   tmOut.setSimplified(parOut);

   cout << "original params\n" << params << endl << endl;
   cout << "seed params\n" << parSeed << endl << endl;
   cout << "new params\n" << parOut << endl << endl;

   cout << "ORIGINAL TRNASITITION MODEL\n" << tm << endl << endl;
   cout << "SEED TRNASITITION MODEL\n" << tmSeed << endl << endl;
   cout << "NEW TRANSITION MODEL\n" << tmOut << endl << endl;

   cout << "ORIGINAL EMISSION TABLE\n" << em.getTable() << endl << endl;
   cout << "SEED EMISSION TABLE\n" << emSeed.getTable() << endl << endl;
   cout << "NEW EMISSION TABLE\n" << emOut.getTable() << endl << endl;

   double tprobOriginal = 0;
   double tvitOriginal = 0;
   double tprobSeed = 0;
   double tprobNew = 0;
   double tvitNew = 0;
   double tprobAD = 0.;
   size_t tdiff = 0;
   size_t tlength = 0;

   delete trainer;
   trainer = NULL;
  
   for (size_t i = 0; i < inputPairs.size(); ++i)
   {
      FTAGModel ftag;
      ftag.setSequences(inputPairs[i].first, inputPairs[i].second,
                        inputWins[i].first, inputWins[i].second);
      ftag.setTransitionModel(tm);
      ftag.setEmissionModel(em);
      double fw = log(ftag.forward());
      double vw = log(ftag.viterbi());
      ftag.setTransitionModel(tmSeed);
      ftag.setEmissionModel(emSeed);
      double fws = log(ftag.forward());
      ftag.setTransitionModel(tmOut);
      ftag.setEmissionModel(emOut);
      double fwo = log(ftag.forward());
      double vwo = log(ftag.viterbi());

      deque<Trace> trace;
      string a, b;
      ftag.viterbiTrace(a, b, trace);
      cout << i << "A) " << a << endl << i << "B) " << b << endl;
      AxtPair apEstimated(a, b);
      AxtPair apOriginal(inputPairs[i].first, inputPairs[i].second);
      size_t adiff = apOriginal.diff(apEstimated);
      
      cout << i << ") Pr[Original] = " << fw
           << " Vit[Orig] = " << vw << " Pr[New] = " << fwo 
           << " Vit[New] = " << vwo << " PrDiff = " << adiff << endl;
      tprobOriginal += fw;
      tvitOriginal += vw;
      tprobSeed += fws;
      tprobNew += fwo;
      tvitNew += vwo;
      tprobAD += (fwo - fw);
      tdiff += adiff;
      tlength += inputPairs[i].first.length();
   }


   FTAGParams parDelta = params - parOut;
   parDelta = parDelta.abs();
   double parMSE = parOut.mse(params);

   cout << "tot) Pr[Original] = " << tprobOriginal
        << " Pr[Seed] = " << tprobSeed << " Pr[New] = " << tprobNew 
        << " AveDelta = " << tprobAD << " tdiff = " << tdiff << endl;
   cout << "   Vit[Oroginal] = " << tvitOriginal
        << " vit[New] = " << tvitNew << endl;

   cout << "Pr Input= " << tprobOriginal << "\nPar_Input= " << params << endl
        << "Par_Delta= " << parDelta << "\nPar_MSE= " << parMSE << endl
        << "Pr = " << tprobNew  << " tdiff = " << tdiff
        << "\nParamEst =\n" << parOut << endl;

   if (logFile)
      makeEmGraph(emTrace, tprobOriginal, params, tlength, logFile);

   return 0;
}

static void makeEmGraph(const FTAGTrain::EMTrace& trace,
                        double prOrig, const FTAGParams& parOrig,
                        size_t numSites,
                        ostream& os)
{
   // iteration prob delta error squared
   os << "#it   numSites  prOriginal  prEst  paramMSE \n";
   for (size_t i = 0; i < trace.size(); ++i)
   {
      double mse = trace[i].second.mse(parOrig);

      os << i << "  " << numSites << "  "
         << prOrig << "  " 
         << trace[i].first << "  "
         << mse;
      os << "\n";
   }
}
