//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cstdlib>
#include <iostream>

#include "ftagtable.h"
#include "dnasubmodel.h"

using namespace std;
using namespace FTAGStates;

double ***** mkTab5(size_t lenA, size_t lenB, size_t win, 
                    size_t states)
{
   double ***** tab = new double ****[lenA];
   for (size_t i = 0; i < lenA; ++i)
   {
      tab[i] = new double***[lenB];
      for (size_t p = 0; p < lenB; ++p)
      {
         tab[i][p] = new double**[win];
         for (size_t w1 = 0; w1 < win; ++w1)
         {
            tab[i][p][w1] = new double*[win];
            for (size_t w2 = 0; w2 < win; ++w2)
            {
               tab[i][p][w1][w2] = new double[states];
            }
         }
      }
   }
   return tab;
}

void delTab5(size_t lenA, size_t lenB, size_t win, 
             size_t states, double ***** tab)
{
   for (size_t i = 0; i < lenA; ++i)
   {
      for (size_t p = 0; p < lenB; ++p)
      {
         for (size_t w1 = 0; w1 < win; ++w1)
         {
            for (size_t w2 = 0; w2 < win; ++w2)
            {
               delete [] tab[i][p][w1][w2];
            }
            
            delete [] tab[i][p][w1];
         }
         delete [] tab[i][p];
      }
      delete [] tab[i];
   }
   delete [] tab;
}

double *** mkTab3(size_t lenA, size_t lenB, size_t states)
{
   double *** tab = new double **[lenA];
   for (size_t i = 0; i < lenA; ++i)
   {
      tab[i] = new double*[lenB];
      for (size_t p = 0; p < lenB; ++p)
      {
         tab[i][p] = new double[states];
      }
   }
   return tab;
}

void delTab3(size_t lenA, size_t lenB,
             size_t states, double *** tab)
{
   for (size_t i = 0; i < lenA; ++i)
   {
      for (size_t p = 0; p < lenB; ++p)
      {
         delete [] tab[i][p];
      }
      delete [] tab[i];
   }
   delete [] tab;
}

void initRand(double*** mtab1, double***** mtab2, double***** mtab3,
              FTAGTable<double, FTAGTableStructs::Sum>& ftab, 
              size_t lenA, size_t lenB, size_t winA,
              size_t winB)
{
   double rv;
   for (size_t i = 0; i < lenA; ++i)
   {
      for (size_t p = 0; p < lenB; ++p)
      {
         // TAB1
         for (size_t s = 0; s < IPStatesSize; ++s)
         {
            rv = drand48();
            mtab1[i][p][IPStates[s]] = rv;
            ftab.set1(i, p, IPStates[s], rv);
         }

         // TAB2
         for (size_t j = 0; j < winA && i + j + 1 < lenA; ++j)
         {
            for (size_t k = j; k < winA && i + k + 1 < lenA; ++k)
            {
               for (size_t s = 0; s < IJKPStatesSize; ++s)
               {
                  rv = drand48();
                  mtab2[i][p][j][k][IJKPStates[s]] = rv;
                  ftab.set2(i, i+1+j, i+1+k, p, IJKPStates[s], rv);
               }
            }
         }
         
         //TAB3
         for (size_t q = 0; q < winB && p + q + 1 < lenB; ++q)
         {
            for (size_t r = q; r < winB && r + q + 1 < lenB; ++r)
            {
               for (size_t s = 0; s < IPQRStatesSize; ++s)
               {
                  rv = drand48();
                  mtab3[i][p][q][r][IPQRStates[s]] = rv;
                  ftab.set3(i, p, p+1+q, p+1+r, IPQRStates[s], rv);
               }
            }
         }
      }
   }
}

bool compare(double*** mtab1, double***** mtab2, double***** mtab3,
             FTAGTable<double, FTAGTableStructs::Sum>& ftab, 
             size_t lenA, size_t lenB, size_t winA,
             size_t winB)
{
   for (size_t i = 0; i < lenA; ++i)
   {
      for (size_t p = 0; p < lenB; ++p)
      {
         // TAB1
         for (size_t s = 0; s < IPStatesSize; ++s)
         {
            double mval = mtab1[i][p][IPStates[s]];
            double fval = ftab.get1(i, p, IPStates[s]);
            if (mval != fval)
            {
               cout << "m1[" << i << "][" << p << "][" << s << "] = "
                    << mval << " =/= " << fval << endl;
               return false;
            }
         }

         // TAB2
         for (size_t j = 0; j < winA && i + 1 + j < lenA; ++j)
         {
            for (size_t k = j; k < winA && i + 1 + k < lenA; ++k)
            {
               for (size_t s = 0; s < IJKPStatesSize; ++s)
               {
                  double mval = mtab2[i][p][j][k][IJKPStates[s]];
                  double fval = ftab.get2(i, i+1+j, i+1+k, p, IJKPStates[s]);
                  if (mval != fval)
                  {
                     cout << "m2[" << i << "][" << p << "][" << j << "]["
                          << k << "][" << s << "] = " << mval 
                          << " =/= " << fval << endl;
                     return false;
                  }
               }
            }
         }
         
         //TAB3
         for (size_t q = 0; q < winB && p + 1 + q < lenB; ++q)
         {
            for (size_t r = q; r < winB && p + 1 + r < lenB; ++r)
            {
               for (size_t s = 0; s < IPQRStatesSize; ++s)
               {
                  double mval = mtab3[i][p][q][r][IPQRStates[s]];
                  double fval = ftab.get3(i, p, p+1+q, p+1+r, IPQRStates[s]);
                  if (mval != fval)
                  {
                     cout << "m3[" << i << "][" << p << "][" << q << "]["
                          << r << "][" << s << "] = " << mval 
                          << " =/= " << fval << endl;
                     return false;
                  }

               }
            }
         }
      }
   }
   return true;
}


int main(int argc, char** argv)
{
   size_t lenA = 55;
   size_t lenB = 17;
   size_t winA = 22;
   size_t winB = 17;
   
   DNASubModel sm;
   sm.setJukesCantor(0.2, 1.2);
   FTAGTable<double, FTAGTableStructs::Sum> ftab;
   ftab.resize(lenA, lenB, winA, winB);
   
   double*** mtab1 = mkTab3(lenA, lenB, S_Max);
   double***** mtab2 = mkTab5(lenA, lenB, winA, S_Max);
   double***** mtab3 = mkTab5(lenA, lenB, winB, S_Max);
   
   initRand(mtab1, mtab2, mtab3, ftab, lenA, lenB, winA, winB);
   compare(mtab1, mtab2, mtab3, ftab, lenA, lenB, winA, winB);

   delTab3(lenA, lenB, S_Max, mtab1);
   delTab5(lenA, lenB, winA, S_Max, mtab2);
   delTab5(lenA, lenB, winB, S_Max, mtab3);

   return 0;
}
