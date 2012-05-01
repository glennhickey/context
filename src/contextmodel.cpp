//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <cmath>

#include "contextmodel.h"
#include "ftagparams.h"

using namespace std;

ContextModel::ContextModel() : 
   _mu(0), _gamma(0), _t(0), 
   _b(0), _d(0),
   _rd(0), _ri(0), _rmd(0), _rmi(0)
{
   _p[0] = _p[1] = _p[2] = _p[3] = 0;
   _q[0] = _q[1] = _q[2] = _q[3] = 0;
   _a[0] = 0;
   _c[0] = 0;
}

ContextModel::~ContextModel()
{
   
}

ContextModel::ContextModel(const ContextModel& model)
{
   *this = model;
}

ContextModel& ContextModel::operator=(const ContextModel& model)
{
   _subModel = model._subModel;

   _mu = model._mu;
   _gamma = model._gamma;
   _t = model._t;
   _b = model._b;
   _d = model._d;
   _a[0] = model._a[0];
   _c[0] = model._c[0];
   _p[0] = model._p[0]; _p[1] = model._p[1]; 
   _p[2] = model._p[2]; _p[3] = model._p[3];
   _q[0] = model._q[0]; _q[1] = model._q[1]; 
   _q[2] = model._q[2]; _q[3] = model._q[3];
   _rd = model._rd;
   _ri = model._ri;
   _rmd = model._rmd;
   _rmi = model._rmi;
   _tab = model._tab;
   return *this;
}

void ContextModel::setSingleAndDouble(const FTAGParams& params)
{
   if (params.isSingleF84() == true)
      setSingleF84(params);
   else setSingleJC(params);
   if (params.isDoubleF84() == true)
      setDoubleF84(params);
   else setDoubleJC(params);
   assert(getTable().verify(0.0005) == true);
}

void ContextModel::setSingleJC(const FTAGParams& params)
{
   _mu = params.getMu();
   _t = params.getT();
   _subModel.setJukesCantor(_mu, _t);
   buildSingleTables();
}

void ContextModel::setSingleF84(const FTAGParams& params)
{
   _p[0] = params.getP0();
   _p[1] = params.getP1();
   _p[2] = params.getP2();
   _p[3] = params.getP3();
   _mu = 0;
   _t = params.getT();
   _a[0] = params.getA();
   _b = params.getB();

   _subModel.setF84(_p[0], _p[1], _p[2], _p[3], _a[0], _b, _t);
   buildSingleTables();
}

void ContextModel::setDoubleJC(const FTAGParams& params)
{
   _mu = params.getMu();
   _gamma = params.getGamma();
   _rd = params.getRD();
   _ri = params.getRI();
   _rmd = params.getRMD();
   _rmi = params.getRMI();
   _t = params.getT();

   computeMXFaster();
}

void ContextModel::setDoubleF84(const FTAGParams& params)
{
   _mu = params.getMu();
   _gamma = params.getGamma();
   _rd = params.getRD();
   _ri = params.getRI();
   _rmd = params.getRMD();
   _rmi = params.getRMI();
   _t = params.getT();
   if (_t != 1.)
   {
      throw("DblF84 only works for t=1");
   }
   _p[0] = params.getP0();
   _p[1] = params.getP1();
   _p[2] = params.getP2();
   _p[3] = params.getP3();
   _a[0] = params.getA();
   _b = params.getB();

   _q[0] = params.getQ0();
   _q[1] = params.getQ1();
   _q[2] = params.getQ2();
   _q[3] = params.getQ3();
   _c[0] = params.getC();
   _d = params.getD();

   // prevent division by zero in generated code. 
   double th = 0.00000000001;
   if (_a[0] < th)
      _a[0] = th;  
   if (_b < th)
      _b = th; 
   if (fabs(_a[0] - _b) < th)
      _a[0] += th;
   if (_d < th)
      _d = th;
   if (fabs(_c[0] - _d) < th)
      _c[0] += th;
   
   computeDblF84Table();
}

void ContextModel::buildSingleTables()
{
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      computeI((DNA)i);
      computeD((DNA)i);
   }

   for (size_t i = 0; i < DNA_Max; ++i)
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         computeM((DNA)i, (DNA)j);
      }
   }
}

// use jukes cantor and summetry properties to not recompute
// unnessarily (can be done elsewhere but the bottleneck was here)
void ContextModel::computeMXFaster()
{      
   computeMD((DNA)0, (DNA)0, (DNA)0);
   computeMD((DNA)0, (DNA)0, (DNA)1);
   computeMD((DNA)0, (DNA)1, (DNA)0);
   computeMD((DNA)0, (DNA)1, (DNA)1);
   computeMD((DNA)0, (DNA)1, (DNA)2);

   if (_rd != _ri || _rmd != _rmi)
   {
      computeMI((DNA)0, (DNA)0, (DNA)0);
      computeMI((DNA)0, (DNA)0, (DNA)1);
      computeMI((DNA)0, (DNA)1, (DNA)0);
      computeMI((DNA)0, (DNA)1, (DNA)1);
      computeMI((DNA)0, (DNA)1, (DNA)2);
   }
   else
   {
      _tab.setMI((DNA)0, (DNA)0, (DNA)0, _tab.getMD((DNA)0, (DNA)0, (DNA)0));
      _tab.setMI((DNA)0, (DNA)0, (DNA)1, _tab.getMD((DNA)0, (DNA)0, (DNA)1));
      // note: flipped
      _tab.setMI((DNA)0, (DNA)1, (DNA)0, _tab.getMD((DNA)0, (DNA)1, (DNA)1));
      _tab.setMI((DNA)0, (DNA)1, (DNA)1, _tab.getMD((DNA)0, (DNA)1, (DNA)0));
      //  .. because insert is upside down delete in for those cases
      _tab.setMI((DNA)0, (DNA)1, (DNA)2, _tab.getMD((DNA)0, (DNA)1, (DNA)2));
   }
 
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            if (i == j)
            {
               if (j == k)
               {
                  _tab.setMD((DNA)i, (DNA)j, (DNA)k, 
                             _tab.getMD((DNA)0, (DNA)0, (DNA)0));
                  _tab.setMI((DNA)i, (DNA)j, (DNA)k, 
                             _tab.getMI((DNA)0, (DNA)0, (DNA)0));
               }
               else
               {
                  _tab.setMD((DNA)i, (DNA)j, (DNA)k, 
                             _tab.getMD((DNA)0, (DNA)0, (DNA)1));
                  _tab.setMI((DNA)i, (DNA)j, (DNA)k, 
                             _tab.getMI((DNA)0, (DNA)0, (DNA)1));
               }
            }
            else if (i == k)
            {
               _tab.setMD((DNA)i, (DNA)j, (DNA)k, 
                          _tab.getMD((DNA)0, (DNA)1, (DNA)0));
               _tab.setMI((DNA)i, (DNA)j, (DNA)k, 
                          _tab.getMI((DNA)0, (DNA)1, (DNA)0));
            }
            else if (j == k)
            {
               _tab.setMD((DNA)i, (DNA)j, (DNA)k, 
                          _tab.getMD((DNA)0, (DNA)1, (DNA)1));
               _tab.setMI((DNA)i, (DNA)j, (DNA)k, 
                          _tab.getMI((DNA)0, (DNA)1, (DNA)1));
            }
            else
            {
               _tab.setMD((DNA)i, (DNA)j, (DNA)k, 
                          _tab.getMD((DNA)0, (DNA)1, (DNA)2));
               _tab.setMI((DNA)i, (DNA)j, (DNA)k, 
                          _tab.getMI((DNA)0, (DNA)1, (DNA)2));
            }
         }
      }
   }
}

void ContextModel::computeD(DNA i)
{
   _tab.setD(i, _subModel.pr(i));
}

void ContextModel::computeI(DNA p)
{
   _tab.setI(p, _subModel.pr(p));
}

void ContextModel::computeM(DNA i, DNA p)
{
   _tab.setM(i, p, _subModel.prFull(i, p));
}
   
// total probability of match delete (todo: table look-up!)
void ContextModel::computeMD(DNA i, DNA k, DNA p)
{
   double total = 0.;

   for (DNA x = (DNA)0; x < DNA_Max; x = (DNA)(int(x) + 1))
   {
      for (DNA y = (DNA)0; y < DNA_Max; y = (DNA)(int(y) + 1))
      {
         // Pr[D|t] = 1/4 rd exp(-rd * t)
         _A[0] = 0.;
         _A[1] = _rd * 0.25; //prFullD(k);
         _A[2] = -_rd;

         // Pr[x|i,t] = same: 1/4 + 3/4exp(-mu * t)  diff: 1/4 -1/4exp(-mu *t)
         _B[0] = 0.25;
         _B[1] = x == i ? 0.75 : -0.25;
         _B[2] = -_mu;

         // Pr[p|x,tmax-t] = 
         // same: 1/4 + 3/4(exp(-mu * tmax) * exp(mu * t))
         // diff: 1/4 - 1/4(exp(-mu * tmax) * exp(mu * t))
         _C[0] = 0.25;
         _C[1] = p == x ? 0.75 * exp(-_mu * _t) : -0.25 * exp(-_mu * _t);
         _C[2] = _mu;
         
         // Pr[y|k,t] = same: 1/4 + 3/4exp(mu * t)  diff: 1/4 - 1/4exp(mu *t)
         _D[0] = 0.25;
         _D[1] = y == k ? 0.75 : -0.25;
         _D[2] = -_mu;

         // Pr[y|x,t] = same: 1/4 + 3/4exp(-gamma)  diff: 1/4 - 1/4exp(-gamma)
         _E[0] = 0.25; 
         _E[1] = y == x ? 0.75 * exp(-_gamma) : -0.25 * exp(-_gamma);
         _E[2] = 0.;

         // integrate at time=t
         double term = intProdFiveExp(_t); 

         // integrate at time=0
         // add difference to running sum
         term -= intProdFiveExp(0.);
         total += term;
      }
   }

   // Condition on deletion happening within time _t
   if (_rd != 0. && _rmd != 0.)
   {
      total /= (1. - exp(-_rd * _t));
   }
   else 
   {
      total = 0.;
   }
   assert(!isnan(total));
   _tab.setMD(i, k, p, total);
}

// context: kp  op: -r
void ContextModel::computeMI(DNA i, DNA p, DNA r)
{
   double total = 0.;
   
   for (DNA x = (DNA)0; x < DNA_Max; x = (DNA)(int(x) + 1))
   {
      for (DNA y = (DNA)0; y < DNA_Max; y = (DNA)(int(y) + 1))
      {
         // Pr[I|t] = 1/4 rd exp(-rd * t)
         _A[0] = 0.;
         _A[1] = _ri * 0.25; //prFullI(y);
         _A[2] = -_ri;

         // Pr[x|i,t] = same: 1/4 + 3/4exp(-mu * t)  diff: 1/4 -1/4exp(-mu *t)
         _B[0] = 0.25;
         _B[1] = x == i ? 0.75 : -0.25;
         _B[2] = -_mu;

         // Pr[p|x,tmax-t] = 
         // same: 1/4 + 3/4(exp(-mu * tmax) * exp(mu * t))
         // diff: 1/4 - 1/4(exp(-mu * tmax) * exp(mu * t))
         _C[0] = 0.25;
         _C[1] = p == x ? 0.75 * exp(-_mu * _t) : -0.25 * exp(-_mu * _t);
         _C[2] = _mu;
         
         // Pr[r|y,t] = same: 1/4 + 3/4exp(mu * t)  diff: 1/4 - 1/4exp(mu *t)
         _D[0] = 0.25;
         _D[1] = y == r ? 0.75 * exp(-_mu * _t) : -0.25 * exp(-_mu * _t);
         _D[2] = _mu;

         // Pr[y|x,t] = same: 1/4 + 3/4exp(-gamma)  diff: 1/4 - 1/4exp(-gamma)
         _E[0] = 0.25; 
         _E[1] = y == x ? 0.75 * exp(-_gamma) : -0.25 * exp(-_gamma);
         _E[2] = 0.;

         // integrate at time=t
         double term = intProdFiveExp(_t); 

         // integrate at time=0
         // add difference to running sum
         term -= intProdFiveExp(0.);
         total += term;
      }
   }

   // Condition on insertion happening within time _t
   if (_ri != 0. && _rmi != 0.)
   {
      total /= (1. - exp(-_ri * _t));
   }
   else
   {
      total = 0.;
   }
   assert(!isnan(total));

   _tab.setMI(i, p, r, total);
}

// compute integral (evaluated at t) of product of five exponential terms:
// (A[0] + A[1]e^(A[2]x)) * (B[0] + B[1]e^(B[2]x)) ....
// not written for speed.. so these values will probably need to be 
// stored in a table or something. 
// Basic idea:
// the product will multiply out to a sum of 2^5 terms, of the form 
// mult * exp(expon * t)
// mult is the product of five values, which can either be "x" or "y",
// and the corresponding expon term is the product of all "y" terms.
// the sum is over all possible values of mult. The integral is computed
// by dividing by expon (by definition).  
double ContextModel::intProdFiveExp(double t) const
{
   // 0 term
   double sum = 0;
   double temp = 0;
   for (unsigned int mask = 0; mask < (1<<5); ++mask)
   {
      double mult = 1.;
      mult *= mask & 1 ? _A[1] : _A[0];
      mult *= mask & (1<<1) ? _B[1] : _B[0];
      mult *= mask & (1<<2) ? _C[1] : _C[0];
      mult *= mask & (1<<3) ? _D[1] : _D[0];
      mult *= mask & (1<<4) ? _E[1] : _E[0];
      
      double expon = 0.;
      expon += mask & 1 ? _A[2] : 0.;
      expon += mask & (1<<1) ? _B[2] : 0.;
      expon += mask & (1<<2) ? _C[2] : 0.;
      expon += mask & (1<<3) ? _D[2] : 0.;
      expon += mask & (1<<4) ? _E[2] : 0.;
      
      // constant term
      if (expon == 0.)
      {
         sum += mult * t;
         temp += mult;
      }
      else
      {
         sum += (mult * exp(expon * t)) / expon;
      }
      assert(!isnan(sum) && !isinf(sum));
   }
   return sum;
}
