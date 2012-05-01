//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _DNASUBMODEL_H
#define _DNASUBMODEL_H

#include <ostream>
#include "ctmarkov.h"

enum DNA {DNA_A = 0, DNA_C, DNA_G, DNA_T, DNA_Max};

class DNASubModel : public CtMarkov
{
public:

   double prTrans(DNA from, DNA to) const { return CtMarkov::prTrans(from, to); }
   double pr(DNA state) const { return CtMarkov::pr(state); }
   double prFull(DNA from, DNA to) const { return CtMarkov::prFull(from, to); }
   
   DNASubModel();
   virtual ~DNASubModel();

   DNASubModel(const DNASubModel& model);
   DNASubModel& operator=(const DNASubModel& model);

   static DNA charToDNA(char c);
   static char DNAToChar(DNA dna);
   static char gap() { return '-'; }

   void setJukesCantor(double mu, double t);
   void setF84(double p0, double p1, double p2, double p3, double a, double b,
               double t);

protected:

};

inline std::ostream& operator<<(std::ostream& os, DNA dna)
{
   os << DNASubModel::DNAToChar(dna);
   return os;
}

#endif
