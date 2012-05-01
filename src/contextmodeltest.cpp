//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "contextmodel.h"
#include "ftagargs.h"
#include "ftagparams.h"

using namespace std;

static 
double sumOneTest(ContextModel& model, 
                double (ContextModel::*prSingle)(DNA) const)
{
   double sum = 0.;
   for (size_t i = 0; i < DNA_Max; ++i) 
   {
      sum += (model.*prSingle)((DNA)i);
   }
   return sum;
}

static
double sumOneTest(ContextModel& model, 
                  double (ContextModel::*prDouble)(DNA, 
                                                   DNA) const)
{
   double sum = 0.;
   for (size_t i = 0; i < DNA_Max; ++i) 
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         sum += (model.*prDouble)((DNA)i, (DNA)j);
      }
   }
   return sum;
}

static
double sumOneTest(ContextModel& model, 
                  double (ContextModel::*prTriple)(DNA, 
                                                   DNA,
                                                   DNA) const)
{
   double sum = 0.;
   for (size_t i = 0; i < DNA_Max; ++i) 
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            sum += (model.*prTriple)((DNA)i, 
                                     (DNA)j,
                                     (DNA)k);         
         }
      }
   }
   return sum;
}

static
double sumOneTestMM(ContextModel& model)
{
   double sum = 0.;
   for (size_t i = 0; i < DNA_Max; ++i) 
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            for (size_t l = 0; l < DNA_Max; ++l)
            {
               sum += 
                  model.prFullM((DNA)i, (DNA)j) * 
                  model.prFullM((DNA)k, (DNA)l);
            }
         }
      }
   }
   return sum;
}



int main(int argc, char** argv)
{
   try{
   FTAGArgs args;
   args.getArgs(argc, argv);

   FTAGParams params;
   ifstream pfile(args._fpFile.c_str());
   if (!pfile)
   {
      cout << "Specify params with fpfile=<path>" << endl;
   }
   pfile >> params;
   params.setSymmetric(args._symmetric);
   params.setSingleF84(args._singleF84);
   params.setDoubleF84(args._doubleF84);

   ContextModel model;
   model.setSingleAndDouble(params);
   
   cout << "A " << model.prFullD(DNA_A) << endl;
   cout << "AA " 
        << model.prFullD(DNA_A) * 
      model.prFullD(DNA_A) 
        << endl;
   cout << "Match(AC) " 
        << model.prFullM(DNA_A, DNA_C)
        << "  Match(AA) " 
        << model.prFullM(DNA_A, DNA_A)
        << endl;
   cout << "Match(AG) " 
        << model.prFullM(DNA_A, DNA_G)
        << "  Match(GG) " 
        << model.prFullM(DNA_G, DNA_G)
        << endl;
   cout << "Match(AT) " 
        << model.prFullM(DNA_A, DNA_T)
        << "  Match(TA) " 
        << model.prFullM(DNA_T, DNA_A)
        << endl;

   cout << "AA -> A-: " << model.prFullMD(DNA_A, 
                                          DNA_A, 
                                          DNA_A) << endl;
   cout << "AC -> A-: " << model.prFullMD(DNA_A, 
                                          DNA_C, 
                                          DNA_A) << endl;

   cout << "A- -> AA: " << model.prFullMI(DNA_A, 
                                          DNA_A, 
                                          DNA_A) << endl;

   cout << "CA -> A-: " << model.prFullMD(DNA_C, 
                                          DNA_A, 
                                          DNA_A) << endl;
   cout << "C- -> AC: " << model.prFullMI(DNA_C, 
                                          DNA_A, 
                                          DNA_C) << endl;

   cout << "CC -> C-: " << model.prFullMD(DNA_C, 
                                          DNA_C, 
                                          DNA_C) << endl;
   cout << "TT -> C-: " << model.prFullMD(DNA_T, 
                                          DNA_T, 
                                          DNA_C) << endl;
   cout << "TG -> C-: " << model.prFullMD(DNA_T, 
                                          DNA_G, 
                                          DNA_C) << endl;


   cout << "Total M Prob= " << sumOneTest(model, &ContextModel::prFullM)
        << endl;
   cout << "Total D Prob= " << sumOneTest(model, &ContextModel::prFullD)
        << endl;

   cout << "Total MD Prob= " << sumOneTest(model, &ContextModel::prFullMD)
        << endl;
   cout << "Total MI Prob= " << sumOneTest(model, &ContextModel::prFullMI)
        << endl;
   cout << "Total MM Prob= " << sumOneTestMM(model)
        << endl;
   }
   catch(string message)
   {
      cerr << message << endl;
   }

   return 0;
}

/*
double oneterm(double b, double a, double bm, double am, double t)
{
   double mult = pow(b, bm);
   double expon = am * -a * t;
   double divis = am * -a;

   return divis ? (mult * exp(expon)) / divis : mult * t;
}

double DNASubModel::temptest(double b, double a, double k, double t) const
{
   double sum = 
      oneterm(b, a, 4, 1, t) +
      oneterm(b, a, 4, 2, t) * 2 +
      oneterm(b, a, 4, 3, t) +
      oneterm(b, a, 4, 2, t) +
      oneterm(b, a, 4, 3, t) * 2 +
      oneterm(b, a, 4, 4, t);
      
   sum *= k;

   _A[0] = 0.;
   _A[1] = b;
   _A[2] = -a;
   _B[0] = _C[0] = _D[0] = b;
   _B[1] = _C[1] = _D[1] = b;
   _B[2] = _C[2] = _D[2] = -a;
   _E[0] = k;
   _E[1] = 0;
   _E[2] = 0;
   
   double evret = intProdFiveExp(t);

   cout << "sum " << sum << "   eval " << evret << endl;
   return sum;
}

*/
