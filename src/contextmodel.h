//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _CONTEXT_MODEL_H
#define _CONTEXT_MODEL_H

#include "dnasubmodel.h"
#include "contextmodeltable.h"

class FTAGParams;

class ContextModel
{
public:
   ContextModel();
   virtual ~ContextModel();

   ContextModel(const ContextModel& model);
   ContextModel& operator=(const ContextModel& model);
   
   void setSingleAndDouble(const FTAGParams& params);
   void setSingleJC(const FTAGParams& params);
   void setSingleF84(const FTAGParams& params);
   void setDoubleJC(const FTAGParams& params);
   void setDoubleF84(const FTAGParams& params);
   void setTable(const ContextModelTable& tab) {_tab = tab;}
   const ContextModelTable& getTable() const {return _tab;}

   // single events
   double prFullD(DNA i) const {return _tab.getD(i);}
   double prFullI(DNA p) const {return _tab.getI(p);}
   double prFullM(DNA i, DNA p) const {return _tab.getM(i, p);}

   double prFullMD(DNA i, DNA k, DNA p) const {return _tab.getMD(i, k, p);}
   
   double prFullMI(DNA i, DNA p, DNA r) const {return _tab.getMI(i, p, r);}

protected:
   
   void buildSingleTables();
   void computeMXFaster();
   void computeD(DNA p);   
   void computeI(DNA i);
   void computeM(DNA i, DNA p);
   void computeMD(DNA i, DNA k, DNA p);
   void computeMI(DNA i, DNA p, DNA r);

   double intProdFiveExp(double lambda) const;

protected:

   DNASubModel _subModel;

   // JC Substitution
   double _mu;
   
   // JC "Cross"
   double _gamma;

   double _t;

   // F84 Substitution
   long double _p[4];
   long double _a[1];
   long double _b;

   // F84 "Cross"
   long double _q[4];
   long double _c[1];
   long double _d;

   double _rd;
   double _ri;
   double _rmd;
   double _rmi;

   ContextModelTable _tab;

   mutable double _A[3];
   mutable double _B[3];
   mutable double _C[3];
   mutable double _D[3];
   mutable double _E[3];

   #include "dblf84/dblf84_dec.h"
};


#endif

// dual sub model here
