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
   FTAGArgs args;
   args.getArgs(argc, argv);
   srand48(args._seed);
   cout << "SEED=" << args._seed << endl;

   ofstream pfile(args._fpFileOut.c_str());
   if (!pfile)
   {
      cerr << "Specify output params with fpout=<path>" << endl;
      exit(1);
   }

   ifstream ifile(args._inFile.c_str());
   if (!ifile)
   {
      cerr << "Specify input alignments with infile=<path>" << endl;
      exit(1);
   }
 
   ifstream seedfile;
   if (args._fpFile.length() > 0)
   {
      seedfile.open(args._fpFile.c_str());
      if (!seedfile)
      {
         cerr << "Specify params with fpfile=<path>" << endl;
         exit(1);
      }
   }
 
   Pairs inputPairs;
   ifile >> inputPairs;

   cout << "input file contains " << inputPairs.size() << " pairs" << endl;
   
   FTAGParams parSeed;
   if (seedfile.is_open())
   {
      seedfile >> parSeed;
   }
   parSeed.setT(1.);
   parSeed.setTFixed(true);
   
   if (args._context == false)
   {
      parSeed.setRMD(0.);
      parSeed.setRMI(0.);
      parSeed.setPCD(1.);
      parSeed.setPCI(1.);
      parSeed.setGamma(1.);
      parSeed.setRMDFixed(true);
      parSeed.setRMIFixed(true);
      parSeed.setPCDFixed(true);
      parSeed.setPCIFixed(true);
      parSeed.setGammaFixed(true);
      parSeed.setQFlatFixed(true);
   }

   if (args._fixA >= 0)
   {
      parSeed.setA(args._fixA);
      parSeed.setAFixed(true);
   }
   if (args._fixB >= 0)
   {
      parSeed.setB(args._fixB);
      parSeed.setBFixed(true);
   }
   if (args._fixC >= 0)
   {
      parSeed.setC(args._fixC);
      parSeed.setCFixed(true);
   }
   if (args._fixD >= 0)
   {
      parSeed.setD(args._fixD);
      parSeed.setDFixed(true);
   }

   parSeed.setSymmetric(args._symmetric);
   parSeed.setSingleF84(args._singleF84);
   parSeed.setDoubleF84(args._doubleF84);
   parSeed.setPLinkedGC(args._plgc);
   parSeed.setPLinkedGA(args._plga);
   parSeed.setQLinkedGC(args._qlgc);
   parSeed.setQLinkedGA(args._qlga);
   parSeed.setPFlatFixed(args._pflat);
   parSeed.setQFlatFixed(args._qflat);
    
   if (!seedfile.is_open())
   {
      parSeed.randomizeNonFixed();
   }
/*
   parSeed.setPhase(FTAGParams::Par_P0, 0, 10);
   parSeed.setPhase(FTAGParams::Par_P1, 0, 10);
   parSeed.setPhase(FTAGParams::Par_P2, 0, 10);
   parSeed.setPhase(FTAGParams::Par_P3, 0, 10);
   parSeed.setPhase(FTAGParams::Par_RD, 0, 10);
   parSeed.setPhase(FTAGParams::Par_RI, 0, 10);
   parSeed.setPhase(FTAGParams::Par_KD, 0, 10);
   parSeed.setPhase(FTAGParams::Par_KI, 0, 10);
   parSeed.setPhase(FTAGParams::Par_PE, 0, 10);
   parSeed.setPhase(FTAGParams::Par_A, 0, 10);
   parSeed.setPhase(FTAGParams::Par_B, 0, 10);
*/
/*
   parSeed.setPhase(FTAGParams::Par_Q0, 3);
   parSeed.setPhase(FTAGParams::Par_Q1, 3);
   parSeed.setPhase(FTAGParams::Par_Q2, 3);
   parSeed.setPhase(FTAGParams::Par_Q3, 3);
   parSeed.setPhase(FTAGParams::Par_C, 2);
   parSeed.setPhase(FTAGParams::Par_D, 1);
*/
/*
   parSeed.setPhase(FTAGParams::Par_RMD, 1, 20);
   parSeed.setPhase(FTAGParams::Par_RMI, 1, 20);
   parSeed.setPhase(FTAGParams::Par_PCD, 1, 20);
   parSeed.setPhase(FTAGParams::Par_PCI, 1, 20);
*/

   FTAGParams offset;
   offset.setAll(args._offset);
   offset.setDirect(FTAGParams::Par_Q0, args._offset/50.);
   offset.setDirect(FTAGParams::Par_Q1, args._offset/50.);
   offset.setDirect(FTAGParams::Par_Q2, args._offset/50.);   
   offset.setDirect(FTAGParams::Par_Q3, args._offset/50.);
                               
   AGDOptimizer optimizer;
   optimizer.setFixed(parSeed);
   optimizer.setRelative(args._relative);
   optimizer.setRandomOrder(args._randOrder);
   optimizer.setSavedInit(parSeed);
   optimizer.setLoopParams(args._numOptTrials, args._maxOptIt, 
                           args._optThreshold, offset); 

   FTAGTrain* trainer = new FTAGTrain();
   EmissionEstimator ee;
   TransitionEstimator te;
   ee.setBias(args._emBias);
   te.setBias(args._tmBias);
   trainer->initialize(optimizer, ee, te);
   trainer->setSequences(inputPairs, args._winSize);
   FTAGTrain::EMTrace emTrace;
   
   cout << "Seed param\n" << parSeed << endl;
   trainer->emLoop(args._maxEMIt, args._emThreshold, args._emConvRepeats,
                   args._emCRThreshold, &emTrace);

   // with a neg. emthreshold, the returned parameters may not be
   // the best, so we scan through each iteration
   double prMin = emTrace[0].first;
   size_t itMin = 0;
   for (size_t it = 1; it < emTrace.size(); ++it)
   {
      if (emTrace[it].first > prMin)
      {
         prMin = emTrace[it].first;
         itMin = it;
      }
   }
   
   FTAGParams parOut = emTrace[itMin].second;

   cout << "estimated param\n" << parOut << endl;
   pfile << parOut;

   return 0;
}
