//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include "pairs.h"
#include "dnasubmodel.h"

using namespace std;

istream& operator >> (istream& is, Pairs& pairs)
{
   if (!is) 
   {
      throw string("File Read Error 1");
   }

   string buffer;
   char linebuffer[512];

   while (is)
   {
       // EAT UP LINES WE DON'T CARE ABOUT
      while(is && !isalnum(is.peek()) && is.peek() != '-')
      {
         is.getline(linebuffer, 512);
      }
      if (is)
      {
         PairAlignment pa;
         is >> buffer >> pa.first >> buffer >> pa.second;
         pairs.push_back(pa);
      }
   }
   
   return is;
}

ostream& operator << (ostream& os, const Pairs& pairs)
{
   for (size_t i = 0; i < pairs.size(); ++i)
   {
      os << i + 1 << "a) " << pairs[i].first 
         << "\n"
         << i + 1 << "b) " << pairs[i].second
         << "\n";
   }
   return os;
}

ostream& writeFasta(ostream& os, const Pairs& pairs, size_t ncols)
{
   for (size_t i = 0; i < pairs.size(); ++i)
   {
      os << ">a" << i + 1 << "\n";
      for (size_t j = 0; j < pairs[i].first.length(); ++j)
      {
         os << pairs[i].first[j];
         if (j && j % ncols == 0)
            os << "\n";
      }
      if (pairs[i].first.length() % ncols)
         os << "\n";

      os << ">b" << i + 1 << "\n";
      for (size_t j = 0; j < pairs[i].second.length(); ++j)
      {
         os << pairs[i].second[j];
         if (j && j % ncols == 0)
            os << "\n";
      }
      if (pairs[i].second.length() % ncols)
         os << "\n";
   }
   
   return os;
}

// should i go in a class? 
FTAGParams estParams(const vector<pair<string, string> >& alignments,
                     bool sym)
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
   size_t base[5] = {0, 0, 0, 0, 0};

   for (size_t i = 0; i < alignments.size(); ++i)
   {
      ind = false;
      ini = false;
      tlen += alignments[i].first.length();

      for (size_t j = 0; j < alignments[i].first.length(); ++j)
      {
         ++base[(int)DNASubModel::charToDNA(alignments[i].first[j])];
         ++base[(int)DNASubModel::charToDNA(alignments[i].second[j])];
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
   double baseTotal = base[0] + base[1] + base[2] + base[3];


   if (sym == true)
   {
      rd = ri = (rd + ri) / 2.;
      kd = ki = (kd + ki) / 2.;
   }

   FTAGParams params;
   params.setT(1.0);
   params.setTFixed(true);
   params.setMu(mu);
   params.setRD(rd / 2.);
   params.setRI(ri / 2.);
   params.setKD(kd);
   params.setKI(ki);
   params.setRMD(rd / 2.);
   params.setRMI(ri / 2.);
   params.setPCD(1.-kd);
   params.setPCI(1.-ki);
   params.setGamma(mu);
   params.setPE(pe);
   params.setP0(base[0] / baseTotal);
   params.setP1(base[1] / baseTotal);
   params.setP2(base[2] / baseTotal);
   params.setP3(base[3] / baseTotal);
   params.setA(0.0001);
   params.setB(mu - 0.0001);
   params.setC(0.0001);
   params.setD(mu - 0.0001);
   params.setQ0(base[0] / baseTotal);
   params.setQ1(base[1] / baseTotal);
   params.setQ2(base[2] / baseTotal);
   params.setQ3(base[3] / baseTotal);
   return params;
}

