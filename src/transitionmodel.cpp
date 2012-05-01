//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cmath>

#include "transitionmodel.h"
#include "ftagstates.h"

using namespace std;
using namespace FTAGStates;

TransitionModel::TransitionModel() : CtMarkov(S_Max)
{

}

TransitionModel::~TransitionModel()
{

}

TransitionModel::TransitionModel(const TransitionModel& tm) 
   : CtMarkov(tm)
{
   
}

TransitionModel& TransitionModel::operator=(const TransitionModel& tm)
{
   this->CtMarkov::operator=(tm);
   return *this;
}

// Set transitions allowing only match contexts (hence simplified).
void TransitionModel::setSimplified(const FTAGParams& tp)
{
   setAllZero();
   double t = tp.getT(); 
   double rd = tp.getRD();
   double ri = tp.getRI();
   double rmd = tp.getRMD();
   double rmi = tp.getRMI();
   double pe = tp.getPE();
   double pcd = tp.getPCD();
   double pci = tp.getPCI();
   double kd = tp.getKD();
   double ki = tp.getKI();
   
   // OUTSIDE
   double rtot = rd + ri + rmd + rmi;
   double rlt;

   // Start
   _Qt.set(S_Start, S_End, pe);
   if (rtot == 0.)
   {
      _Qt.set(S_Start, S_Mx, (1.-pe));
   }
   else
   {
      _Qt.set(S_Start, S_Mx, (1.-pe) * exp(-t * rtot));
      _Qt.set(S_Start, S_Dx, (1.-pe) * (1.-exp(-t * rtot)) * (rd / rtot));
      _Qt.set(S_Start, S_Ix, (1.-pe) * (1.-exp(-t * rtot)) * (ri / rtot));
      _Qt.set(S_Start, S_MDx, (1.-pe) * (1.-exp(-t * rtot)) * (rmd / rtot));
      _Qt.set(S_Start, S_MIx, (1.-pe) * (1.-exp(-t * rtot)) * (rmi / rtot));
   }
   // Match
   _Qt.set(S_Mx, S_End, pe);
   if (rtot == 0.)
   {
      _Qt.set(S_Mx, S_Mx, (1.-pe));
   }
   else
   {
      _Qt.set(S_Mx, S_Mx, (1.-pe) * exp(-t * rtot));
      _Qt.set(S_Mx, S_Dx, (1.-pe) * (1.-exp(-t * rtot)) * (rd / rtot));
      _Qt.set(S_Mx, S_Ix, (1.-pe) * (1.-exp(-t * rtot)) * (ri / rtot));
      _Qt.set(S_Mx, S_MDx, (1.-pe) * (1.-exp(-t * rtot)) * (rmd / rtot));
      _Qt.set(S_Mx, S_MIx, (1.-pe) * (1.-exp(-t * rtot)) * (rmi / rtot));
   }

   // Delete
   rlt = rtot - rd;
   if (rd == 0.)
   {
      _Qt.set(S_Dx, S_End, 1.);
   }
   else
   {
      _Qt.set(S_Dx, S_End, pe);
      _Qt.set(S_Dx, S_Mx, (1.-pe) * (1.-kd) * exp(-t * rlt));
      _Qt.set(S_Dx, S_Dx, (1.-pe) * kd);
      _Qt.set(S_Dx, S_Ix, (1.-pe) * (1.-kd) * (1.-exp(-t * rlt)) * (ri / rlt));
      _Qt.set(S_Dx, S_MDx, (1.-pe) * (1.-kd) * (1.-exp(-t * rlt)) * (rmd/ rlt));
      _Qt.set(S_Dx, S_MIx, (1.-pe) * (1.-kd) * (1.-exp(-t * rlt)) * (rmi/ rlt));
   }

   // Insert
   rlt = rtot - ri;
   if (ri == 0.)
   {
      _Qt.set(S_Ix, S_End, 1.);
   }
   else
   {
      _Qt.set(S_Ix, S_End, pe);
      _Qt.set(S_Ix, S_Mx, (1.-pe) * (1.-ki) * exp(-t * rlt));
      _Qt.set(S_Ix, S_Dx, (1.-pe) * (1.-ki) * (1.-exp(-t * rlt)) * (rd / rlt));
      _Qt.set(S_Ix, S_Ix, (1.-pe) * ki);
      _Qt.set(S_Ix, S_MDx, (1.-pe) * (1.-ki) * (1.-exp(-t * rlt)) * (rmd/ rlt));
      _Qt.set(S_Ix, S_MIx, (1.-pe) * (1.-ki) * (1.-exp(-t * rlt)) * (rmi/ rlt));
   }

   // Match-Delete
   if (rmd == 0.)
   {
      _Qt.set(S_MDx, S_End, 1.);
   }
   else
   {
      _Qt.set(S_MDx, S_End, pe);
      _Qt.set(S_MDx, S_Mx, (1.-pe) * exp(-t * rtot));
      _Qt.set(S_MDx, S_Dx, (1.-pe) * (1.-exp(-t * rtot)) * (rd / rtot));
      _Qt.set(S_MDx, S_Ix, (1.-pe) * (1.-exp(-t * rtot)) * (ri / rtot));
      _Qt.set(S_MDx, S_MDx, (1.-pe) * (1.-exp(-t * rtot)) * (rmd / rtot));
      _Qt.set(S_MDx, S_MIx, (1.-pe) * (1.-exp(-t * rtot)) * (rmi / rtot));
   }

   // Match-Insert
   if (rmi == 0.)
   {
      _Qt.set(S_MIx, S_End, 1.);
   }
   else
   {
      _Qt.set(S_MIx, S_End, pe);
      _Qt.set(S_MIx, S_Mx, (1.-pe) * exp(-t * rtot));
      _Qt.set(S_MIx, S_Dx, (1.-pe) * (1.-exp(-t * rtot)) * (rd / rtot));
      _Qt.set(S_MIx, S_Ix, (1.-pe) * (1.-exp(-t * rtot)) * (ri / rtot));
      _Qt.set(S_MIx, S_MDx, (1.-pe) * (1.-exp(-t * rtot)) * (rmd / rtot));
      _Qt.set(S_MIx, S_MIx, (1.-pe) * (1.-exp(-t * rtot)) * (rmi / rtot));
   }

   // INSIDE 

   // Match-Delete Start
   _Qt.set(S_MDx, S_Close, pcd);
   _Qt.set(S_MDx, S_MDyM, 1.-pcd);
   
   // Match-Delete Continue
   _Qt.set(S_MDyM, S_Close, pcd);
   _Qt.set(S_MDyM, S_MDyM, 1.-pcd);

   // Match-Insert Start
   _Qt.set(S_MIx, S_Close, pci);
   _Qt.set(S_MIx, S_MIzM, 1.-pci);
   
   // Match-Insert Continue
   _Qt.set(S_MIzM, S_Close, pci);
   _Qt.set(S_MIzM, S_MIzM, 1.-pci);
   
   assert(verify() == true);
}

void TransitionModel::scale(const vector<double>& svec)
{
   for (size_t i = 0; i < svec.size(); ++i)
   {
      for (size_t j = 0; j < svec.size(); ++j)
      {
         _Qt.set(i, j, svec[i] * _Qt.get(i, j));
      }
   }
}                                    

void TransitionModel::setMatchOnly(double p)
{
   setAllZero();
   
   _Qt.set(S_Start, S_Mx, p);
   _Qt.set(S_Start, S_End, 1-p);
   _Qt.set(S_Mx, S_Mx, p);
   _Qt.set(S_Mx, S_End, 1-p);
   
}

bool TransitionModel::verify(double threshold) const
{
   for (size_t i = 0; i < S_Max; ++i)
   {
      double totalInside = 0.;
      double totalOutside = 0.;

      for (size_t j = 0; j < InsideToStatesSize; ++j)
      {
         assert(!isnan(_Qt.get((State)i, InsideToStates[j])));
         totalInside += _Qt.get((State)i, InsideToStates[j]);
      }
      for (size_t j = 0; j < OutsideToStatesSize; ++j)
      {
         assert(!isnan(_Qt.get((State)i, OutsideToStates[j])));
         totalOutside += _Qt.get((State)i, OutsideToStates[j]);
      }
      
      if (fabs(1. - totalInside) > threshold && 
          fabs(0. - totalInside) > threshold)
      {
         return false;
      }
      if (fabs(1. - totalOutside) > threshold && 
          fabs(0. - totalOutside) > threshold)
      {
         return false;
      }
   
      double total = totalInside + totalOutside;
      if (fabs(0. - total) > threshold &&
          fabs(1. - total) > 2. * threshold &&
          fabs(2. - total) > 2. * threshold)
      {
         return false;
      }
   }
   return true;
}

ostream& operator<<(ostream& os, const TransitionModel& tm)
{
   for (size_t i = 0; i < S_Max; ++i)
   {
      if ((State)i != S_End && (State)i != S_Close)
      {
         for (size_t j = 0; j < S_Max; ++j)
         {
            if ((State)j != S_Start)
            {
               os << "Pr[" << (State)i << "->" << (State)j << "] = "
                  << tm.prTrans((State)i, (State)j) << endl;
            }
         }
      }
   }
   return os;
}
