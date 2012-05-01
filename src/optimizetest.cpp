//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <fstream>
#include "agdoptimizer.h"
#include "ftagargs.h"
using namespace std;

int main(int argc, char** argv)
{
   FTAGArgs args;
   args.getArgs(argc, argv);
   
   srand48(args._seed);
   cout << "SEED=" << args._seed << endl;

   FTAGParams params;
   try
   {
      ifstream pfile(args._fpFile.c_str());
      if (!pfile)
      {
         cout << "Specify params with fpfile=<path>" << endl;
      }
      pfile >> params;
   }
   catch(string message)
   {
      cout << message << endl;
      return 1;
   }
   
   ContextModel em;
   em.setDoubleJC(params);
   TransitionModel tm;
   tm.setSimplified(params);

   AGDOptimizer optimizer;
   optimizer.setFixed(params);
   optimizer.setModels(em, tm);

   FTAGParams offset(params);
   offset.setAll(args._offset);
   FTAGParams parSeed(params);
   parSeed.setSymmetric(args._symmetric);
   parSeed.randomizeNonFixed();

   optimizer.setSavedInit(parSeed);
   optimizer.setRelative(args._relative);
   optimizer.setRandomOrder(args._randOrder);
   optimizer.setLoopParams(args._numOptTrials, args._maxOptIt, 
                           args._optThreshold, offset); 

   double mse = optimizer.optimize();

   double mymse = params.mse(optimizer.getParams());

   TransitionModel tmEst;
   tmEst.setSimplified(optimizer.getParams());
   cout << "ORIGINAL\n" << tm << "\nESTIMATED\n" << tmEst << endl << endl;

   cout << "ORIGINAL\n" << params
        << "\nNEW\n" << optimizer.getParams()
        << "\nDELTA\n" << params - optimizer.getParams() 
        << endl;
   cout << "MSE=" << mse << endl;


   return 0;
}
