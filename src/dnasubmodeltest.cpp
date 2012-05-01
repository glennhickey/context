//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include "dnasubmodel.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

pair<double, double> getJCMat(double mu, double t, size_t size)
{
   double e = 2.71828183;
   double f = 1. / (double)size;
   double c = 1. - f;
   double stay = f + c * pow(e, -t * mu);
   double tran = f - f * pow(e, -t * mu);
   return pair<double, double>(stay, tran);
}

int main(int argc, char** argv)
{
   double t = 1.3;
   double mu = 0.2;
   double e = 2.71828183;

   DNASubModel model;

   model.setJukesCantor(mu, t);
   
   cout << "Qt \n";
   cout << model.getProbMatrix();
   cout << "Q0 \n";
   cout << model.getRateMatrix();
   cout << endl;
   
   double stayProb = getJCMat(mu, t, 4).first;
   double transProb = getJCMat(mu, t, 4).second;

   cout << "should be < " << transProb << " " << stayProb << " >" << endl;
   
   return 0;
}
