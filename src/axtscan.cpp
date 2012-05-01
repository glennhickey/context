//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <fstream>
#include <iostream>
#include <algorithm>
#include <locale>
#include "axtreader.h"
#include "ftagtrain.h"
#include "ftagargs.h"

using namespace std;
using namespace FTAGTableStructs;
using namespace FTAGStates;

vector<pair<string, string> > filter(const vector<AxtPair>& aps,
                                     size_t minLen, size_t maxLen)
{
   vector<pair<string, string> > alignments;
   for (size_t i = 0; i < aps.size(); ++i)
   {
      if (aps[i].getLength() <= maxLen && aps[i].getLength() >= minLen &&
          aps[i].hasN() == false)
      {
         alignments.push_back(aps[i].getAlignment());
      }
   }
   return alignments;
}

static FTAGParams estParams(const vector<pair<string, string> >& alignments);

static void makeEmGraph(const FTAGTrain::EMTrace& trace,
                        double prOrig, const FTAGParams& parOrig,
                        size_t numSites,
                        ostream& os);

static void makeLenDistSimpleSymmetric(const deque<Trace>& vTrace, 
                                       vector<size_t>& dist,
                                       vector<size_t>& dist_nc);

int main(int argc, char** argv)
{
   FTAGArgs args;
   args.getArgs(argc, argv);
   
   srand48(args._seed);
   cout << "SEED=" << args._seed << endl;

   ofstream logFile;
   FTAGTrain::EMTrace emTrace;
   if (!args._outFile.empty())
   {
      logFile.open(args._outFile.c_str());
   }

   ofstream distFile;
   vector<size_t> dist, dist_nc;
   if (!args._distFile.empty())
   {
      distFile.open(args._distFile.c_str());
      dist = vector<size_t>(args._maxLength, 0);
      dist_nc = vector<size_t>(args._maxLength, 0);
   }

   AxtReader ar;
   FTAGParams parSeed;
   
   try
   {
      ar.read(args._inFile);
      ifstream pfile(args._fpFile.c_str());
      if (!pfile)
      {
         cout << "Specify params with fpfile=<path>" << endl;
      }
      pfile >> parSeed;
   }
   catch(string message)
   {
      cout << message << endl;
      return 1;
   }
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


//   vector<pair<string, string> > alignments = filter(ar.getAlignments(),
//                                                     minLength, maxLength);
   const std::vector<AxtPair>& wholeThing = ar.getAlignments();
   vector<pair<string, string> > allAlignments;
   for (size_t i = 0; i < wholeThing.size(); ++i)
   {
      allAlignments.push_back(wholeThing[i].getAlignment());
   }
   vector<pair<string, string> > alignments 
      = ar.sample(args._numPairs, args._maxLength / 2, args._maxLength, 0, 
                  false);

   for (size_t i = 0; i < alignments.size(); ++i)
   {
      cout << i << "a) " << alignments[i].first << endl
           << i << "b) " << alignments[i].second << endl;
   }
   FTAGParams parAll = estParams(allAlignments);
   FTAGParams parInput = estParams(alignments);
   ContextModel emInput;
   emInput.setDoubleJC(parInput);
   TransitionModel tmInput;
   tmInput.setSimplified(parInput);

   parSeed.setSymmetric(args._symmetric);
   parSeed.setUniGap(args._unigap);
   parSeed.randomizeNonFixed();

   cout << "From Whole File: " << parAll << endl;
   cout << "From input: " << parInput << endl;
   cout << "Random Seed: " << parSeed << endl;

   FTAGParams offset;
   offset.setAll(args._offset);

   AGDOptimizer optimizer;
   optimizer.setRelative(args._relative);
   optimizer.setSavedInit(parSeed);
   optimizer.setLoopParams(args._numOptTrials, args._maxOptIt, 
                           args._optThreshold, offset);

   EmissionEstimator ee;
   TransitionEstimator te;
   ee.setBias(args._emBias);
   te.setBias(args._tmBias);

   FTAGTrain* trainer = new FTAGTrain();
   trainer->initialize(optimizer, ee, te);
   trainer->setSequences(alignments, args._winSize);
   trainer->emLoop(args._maxEMIt, args._emThreshold, args._emConvRepeats,
                   args._emCRThreshold, &emTrace);

   FTAGParams params = trainer->getParams();

   TransitionModel tm;
   ContextModel em;
   tm.setSimplified(params);
   em.setDoubleJC(params);

   delete trainer;
   trainer = NULL;
   double pr = 0.;
   double prInput = 0.;
   double tdiff = 0;
   size_t tlength = 0;

   for (size_t i = 0; i < alignments.size(); ++i)
   {
      FTAGModel ftag;
      ftag.setSequences(alignments[i].first, alignments[i].second,
                        args._winSize, args._winSize);

      ftag.setTransitionModel(tm);
      ftag.setEmissionModel(em);
      pr += log(ftag.forward());
      string a,b;
      deque<Trace> trace;
      ftag.viterbi();
      ftag.viterbiTrace(a, b, trace);
      cout << i << "a) " << a << endl 
           << i << "b) " << b << endl;
      AxtPair apEstimated(a, b);
      AxtPair apOriginal(alignments[i].first, alignments[i].second);
      size_t adiff = apOriginal.diff(apEstimated);

      ftag.setTransitionModel(tmInput);
      ftag.setEmissionModel(emInput);
      prInput += log(ftag.forward());

      tdiff += adiff;
      tlength += alignments[i].first.length();

      makeLenDistSimpleSymmetric(trace, dist, dist_nc);
   }

   FTAGParams parDelta = parInput - params;
   parDelta = parDelta.abs();
   double parMSE = params.mse(parInput);

   cout << "Pr Input= " << prInput << "\nPar_Input= " << parInput << endl
        << "Par_Delta= " << parDelta << "\nPar_MSE= " << parMSE << endl
        << "Num Alignments=" << alignments.size() 
        << " Length=[" << args._maxLength/2 << ", " << args._maxLength << "]"
        << " Pr = " << pr  << " tdiff = " << tdiff
        << "\nParamEst =\n" << params << endl;

   if (logFile)
      makeEmGraph(emTrace, prInput, params, tlength, logFile);

   if (distFile)
   {
      distFile << "len  freq  tot  freqNC  totNC\n";
      for (size_t i = 1; i < dist.size(); ++i)
      {
         distFile << i << "  " << dist[i] << "  " << dist[0] << "  "
                  << dist_nc[i] << "  " << dist_nc[0] << "\n";
      }
   }
   return 0;
}

static FTAGParams estParams(const vector<pair<string, string> >& alignments)
{
   size_t mc = 0;
   size_t ic = 0;
   size_t dc = 0;
   size_t is = 0;
   size_t ds = 0;
   size_t ta = 0;
   size_t tb = 0;
   size_t ts = 0;
   bool ind;
   bool ini;
   size_t tlen = 0;

   for (size_t i = 0; i < alignments.size(); ++i)
   {
      ind = false;
      ini = false;
      tlen += alignments[i].first.length();

      for (size_t j = 0; j < alignments[i].first.length(); ++j)
      {
         if (alignments[i].second[j] == '-')
         {
            ++tb;
            if (ind  == false)
            {
               ++dc;
               ind = true;
               ini = false;
            }
            ++ds;
         }
         else if (alignments[i].first[j] == '-')
         {
            ++ta;
            if (ini == false)
            {
               ++ic;
               ini = true;
               ind = false;
            }
            ++is;
         }
         else
         {
            ++ta;
            ++tb;
            ++ts;
            ind = false;
            ini = false;
            if (toupper(alignments[i].first[j]) != 
                toupper(alignments[i].second[j]))
            {
         if (toupper(alignments[i].first[j]) != 'N' &&
             toupper(alignments[i].second[j]) != 'N')

               ++mc;
            }
         }
      }
   }

   double pe = (double)alignments.size() / (double)tlen;
   double mexp = (double)mc / (double)ts;
   double rd = (double)dc / (double)ta;
   double ri = (double)ic / (double)tb;
   double kd = dc == 0 ? 0 : 1. - (double)dc / (double)ds;
   double ki = ic == 0 ? 0 : 1. - (double)ic / (double)is;

   double mu = -log(1. - (4./3.) * mexp);

   FTAGParams params;
   params.setT(1.0);
   params.setTFixed(true);
   params.setMu(mu);
   params.setRD(rd);
   params.setRI(ri);
   params.setKD(kd);
   params.setKI(ki);
   params.setRMD(0);
   params.setRMDFixed(true);
   params.setRMI(0);
   params.setRMIFixed(true);
   params.setPCD(1);
   params.setPCDFixed(true);
   params.setPCI(0);
   params.setPCIFixed(true);
   params.setGamma(0);
   params.setGammaFixed(true);
   params.setPE(pe);

   return params;
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

static void makeLenDistSimpleSymmetric(const deque<Trace>& vTrace, 
                                       vector<size_t>& dist,
                                       vector<size_t>& dist_nc)
{
   if (dist.size() == 0)
      return;

   for (size_t i = 0; i < vTrace.size(); ++i)
   {
      if (vTrace[i]._s == S_Dx)
      {
         size_t s = i;
         while (i+1 < vTrace.size() && vTrace[i+1]._s == S_Dx) 
            ++i;
         ++dist_nc[i - s + 1];
         ++dist_nc[0]; // store total count in 0
      }
      else if (vTrace[i]._s == S_Ix)
      {
         size_t s = i;
         while (i+1 < vTrace.size() && vTrace[i+1]._s == S_Ix) 
            ++i;
         ++dist_nc[i - s + 1];
         ++dist_nc[0];
      }
      else if (vTrace[i]._s == S_MDx)
      {
         assert(vTrace[i].isIJKPState() == true);
         size_t s = i;
         do {
            ++i;
         } while (vTrace[i]._s != S_MDx);
         assert(vTrace[i]._s == S_MDx && vTrace[i].isIPState());
         ++dist[i - s];
         ++dist[0];
      }
      else if (vTrace[i]._s == S_MIx)
      {
         assert(vTrace[i].isIJKPState() == true);
         size_t s = i;
         do {
            ++i;
         } while (vTrace[i]._s != S_MIx);
         assert(vTrace[i]._s == S_MIx && vTrace[i].isIPState());
         ++dist[i - s];
         ++dist[0];
      }
   }
}
