//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGPARAMS_H
#define _FTAGPARAMS_H

#include <cstdlib>
#include <iostream>

class FTAGParams
{
public:

   enum Param { Par_T, Par_Mu, Par_A, Par_B, Par_P0, Par_P1, Par_P2, Par_P3, 
                Par_Gamma, Par_C, Par_D, Par_Q0, Par_Q1, Par_Q2, Par_Q3,
                Par_RD, Par_RI,
                Par_RMD, Par_RMI, Par_PE, Par_PCD, Par_PCI,
                Par_KD, Par_KI, Par_Max };

   FTAGParams();
   FTAGParams(const FTAGParams& params);
   ~FTAGParams();
   FTAGParams& operator=(const FTAGParams& params);
   FTAGParams operator-(const FTAGParams& params) const;
   FTAGParams operator/(const FTAGParams& params) const;
   FTAGParams operator*(const FTAGParams& params) const;
   bool operator==(const FTAGParams& params) const;
   bool operator!=(const FTAGParams& params) const{
      return !(*this == params);
   }
   FTAGParams abs() const;
   double mse(const FTAGParams& params) const;
   void set(Param param, double val);
   void setDirect(Param param, double val) {
      _params[param] = val;
   }
   void setAll(double val);
   void setFixed(Param param, bool val);
   void setClamp(Param param, double lower, double upper);
   void randomizeNonFixed();

   void setSymmetric(bool sym) {
      _isSymmetric = sym;
   }
   void setUniGap(bool ug) {
      _isUniGap = ug;
   }
   void setSingleF84(bool sf);
   void setDoubleF84(bool df);
   void setPLinkedGC(bool pgc)
   {
      _isPLinkedGC = pgc;
      if (pgc) _isPLinkedGA = false;
   }
   void setQLinkedGC(bool qgc)
   {
      _isQLinkedGC = qgc;
      if (qgc) _isQLinkedGA = false;
   }
   void setPLinkedGA(bool pga)
   {
      _isPLinkedGA = pga;
      if (pga) _isPLinkedGC = false;
   }
   void setQLinkedGA(bool qga)
   {
      _isQLinkedGA = qga;
      if (qga) _isQLinkedGC = false;
   }
   void setPhase(Param param, size_t imin, size_t imax = 99999999) {
      _phases[0][(size_t)param] = imin;
      _phases[1][(size_t)param] = imax;
   }
   bool isSymmetric() const { 
      return _isSymmetric; 
   }
   bool isUniGap() const {
      return _isUniGap;
   }
   bool isSingleF84() const {
      return _singleF84;
   }
   bool isDoubleF84() const {
      return _doubleF84;
   }
   double get(Param param) const{
      return _params[param];
   }
   bool isFixed(Param param) const {
      return _fixed[param];
   }
   bool isDependent(Param param) const;

   std::pair<double, double> getClamp(Param param) const {
      return std::pair<double, double>(_clamps[0][param], _clamps[1][param]);
   }

   double getT() const { return _params[Par_T]; }
   double getMu() const { return _params[Par_Mu]; }
   double getA() const { return _params[Par_A]; }
   double getB() const { return _params[Par_B]; }
   double getP0() const { return _params[Par_P0]; }
   double getP1() const { return _params[Par_P1]; }
   double getP2() const { return _params[Par_P2]; }
   double getP3() const { return _params[Par_P3]; }
   double getGamma() const { return _params[Par_Gamma]; }
   double getC() const { return _params[Par_C]; }
   double getD() const { return _params[Par_D]; }
   double getQ0() const { return _params[Par_Q0]; }
   double getQ1() const { return _params[Par_Q1]; }
   double getQ2() const { return _params[Par_Q2]; }
   double getQ3() const { return _params[Par_Q3]; }
   double getRD() const { return _params[Par_RD]; }
   double getRI() const { return _params[Par_RI]; }
   double getRMD() const { return _params[Par_RMD]; }
   double getRMI() const { return _params[Par_RMI]; }
   double getPE() const { return _params[Par_PE]; }
   double getPCD() const { return _params[Par_PCD]; }
   double getPCI() const { return _params[Par_PCI]; }
   double getKD() const { return _params[Par_KD]; }
   double getKI() const { return _params[Par_KI]; }
   std::pair<size_t, size_t> getPhase(Param param) const { 
      return std::pair<size_t, size_t>(_phases[0][param], _phases[1][param]);
   }

   void setT(double t) {
      _params[Par_T] = clamp(Par_T, t); 
   }
   void setMu(double mu) {
      _params[Par_Mu] = clamp(Par_Mu, mu); 
   }
   void setA(double a);
   void setB(double b);
   void setP0(double p0);
   void setP1(double p1);
   void setP2(double p2);
   void setP3(double p3);
   void setPFlatFixed(bool ff) {
      if (ff) {
         setDirect(Par_P0, 0.25); setP0Fixed(ff);
         setDirect(Par_P1, 0.25); setP1Fixed(ff);
         setDirect(Par_P2, 0.25); setP2Fixed(ff);
         setDirect(Par_P3, 0.25); setP3Fixed(ff);
      }
   }
   void setGamma(double gamma) { 
      _params[Par_Gamma] = clamp(Par_Gamma, gamma); 
   }
   void setC(double c);
   void setD(double d); 
   void setQ0(double q0);
   void setQ1(double q1);
   void setQ2(double q2);
   void setQ3(double q3);

   void setQFlatFixed(bool ff) {
      if (ff) {
         setDirect(Par_Q0, 0.25); setQ0Fixed(ff);
         setDirect(Par_Q1, 0.25); setQ1Fixed(ff);
         setDirect(Par_Q2, 0.25); setQ2Fixed(ff);
         setDirect(Par_Q3, 0.25); setQ3Fixed(ff);
      }
   }
   void normalizeP();
   void normalizeQ();
   void randomSwapP();
   void randomSwapQ();

   void setRD(double rd) { 
      _params[Par_RD] = clamp(Par_RD, rd); 
      if (_isSymmetric) _params[Par_RI] = _params[Par_RD];
   }
   void setRI(double ri) { 
      if (!_isSymmetric) _params[Par_RI] = clamp(Par_RI, ri); 
   }
   void setRMD(double rmd) {
      _params[Par_RMD] = clamp(Par_RMD, rmd); 
      if (_isSymmetric) _params[Par_RMI] = _params[Par_RMD];
   }
   void setRMI(double rmi) { 
      if (!_isSymmetric) _params[Par_RMI] = clamp(Par_RMI, rmi); 
   }
   void setPE(double pe) {
      _params[Par_PE] = clamp(Par_PE, pe); 
   }
   void setPCD(double pcd) {
      if (!_isUniGap) {
         _params[Par_PCD] = clamp(Par_PCD, pcd); 
         if (_isSymmetric) _params[Par_PCI] = _params[Par_PCD];
      }
   }
   void setPCI(double pci) { 
      if (!_isUniGap && !_isSymmetric) _params[Par_PCI] = clamp(Par_PCI, pci); 
   }
   void setKD(double kd) {
      _params[Par_KD] = clamp(Par_KD, kd); 
      if (_isSymmetric) _params[Par_KI] = _params[Par_KD];
      if (_isUniGap) _params[Par_PCD] = 1. - _params[Par_KD];
      if (_isUniGap && _isSymmetric) _params[Par_PCI] = 1. - _params[Par_KD];
   }
   void setKI(double ki) { 
      if (!_isSymmetric) _params[Par_KI] = clamp(Par_KI, ki); 
      else if (_isUniGap) _params[Par_PCI] = 1. - _params[Par_KI];
   }

   bool isTFixed() const { return _fixed[Par_T]; }
   bool isMuFixed() const { return _fixed[Par_Mu]; }
   bool isAFixed() const { return _fixed[Par_A]; }
   bool isBFixed() const { return _fixed[Par_B]; }
   bool isP0Fixed() const { return _fixed[Par_P0]; }
   bool isP1Fixed() const { return _fixed[Par_P1]; }
   bool isP2Fixed() const { return _fixed[Par_P2]; }
   bool isP3Fixed() const { return _fixed[Par_P3]; }
   bool isGammaFixed() const { return _fixed[Par_Gamma]; }
   bool isCFixed() const { return _fixed[Par_C]; }
   bool isDFixed() const { return _fixed[Par_D]; }
   bool isQ0Fixed() const { return _fixed[Par_Q0]; }
   bool isQ1Fixed() const { return _fixed[Par_Q1]; }
   bool isQ2Fixed() const { return _fixed[Par_Q2]; }
   bool isQ3Fixed() const { return _fixed[Par_Q3]; }
   bool isRDFixed() const { return _fixed[Par_RD]; }
   bool isRIFixed() const { return _fixed[Par_RI]; }
   bool isRMDFixed() const { return _fixed[Par_RMD]; }
   bool isRMIFixed() const { return _fixed[Par_RMI]; }
   bool isPEFixed() const { return _fixed[Par_PE]; }
   bool isPCDFixed() const { return _fixed[Par_PCD]; }
   bool isPCIFixed() const { return _fixed[Par_PCI]; }
   bool isKDFixed() const { return _fixed[Par_KD]; }
   bool isKIFixed() const { return _fixed[Par_KI]; }

   void setTFixed(bool t) {
      _fixed[Par_T] = t; 
   }
   void setMuFixed(bool mu) { 
      _fixed[Par_Mu] = mu; 
   }
   void setAFixed(bool a) {
      _fixed[Par_A] = a;
   }
   void setBFixed(bool b) {
      _fixed[Par_B] = b;
   }
   void setP0Fixed(bool p0) { 
      _fixed[Par_P0] = p0; 
   }
   void setP1Fixed(bool p1) { 
      _fixed[Par_P1] = p1; 
   }
   void setP2Fixed(bool p2) { 
      _fixed[Par_P2] = p2; 
   }
   void setP3Fixed(bool p3) { 
      _fixed[Par_P3] = p3; 
   }
   void setGammaFixed(bool gamma) { 
      _fixed[Par_Gamma] = gamma; 
   }
   void setCFixed(bool c) {
      _fixed[Par_C] = c;
   }
   void setDFixed(bool d) {
      _fixed[Par_D] = d;
   }
   void setQ0Fixed(bool q0) { 
      _fixed[Par_Q0] = q0; 
   }
   void setQ1Fixed(bool q1) { 
      _fixed[Par_Q1] = q1; 
   }
   void setQ2Fixed(bool q2) { 
      _fixed[Par_Q2] = q2; 
   }
   void setQ3Fixed(bool q3) { 
      _fixed[Par_Q3] = q3; 
   }
   void setRDFixed(bool rd) { 
      _fixed[Par_RD] = rd; 
      if (_isSymmetric) _fixed[Par_RI] = rd;
   }
   void setRIFixed(bool ri) { 
      if (!_isSymmetric) _fixed[Par_RI] = ri; 
   }
   void setRMDFixed(bool rmd) { 
      _fixed[Par_RMD] = rmd; 
      if (_isSymmetric) _fixed[Par_RMI] = rmd;
   }
   void setRMIFixed(bool rmi) { 
      if (!_isSymmetric) _fixed[Par_RMI] = rmi; 
   }
   void setPEFixed(bool pe) { 
      _fixed[Par_PE] = pe; 
   }
   void setPCDFixed(bool pcd) { 
      _fixed[Par_PCD] = pcd; 
      if (_isSymmetric) _fixed[Par_PCI] = pcd;
   }
   void setPCIFixed(bool pci) { 
      if (!_isSymmetric) _fixed[Par_PCI] = pci; 
   }
   void setKDFixed(bool kd) {
      _fixed[Par_KD] = kd; 
      if (_isSymmetric) _fixed[Par_KI] = kd;
   }
   void setKIFixed(bool ki) { 
      if (!_isSymmetric) _fixed[Par_KI] = ki; 
   }

   std::string asRow() const;

protected:

   double clamp(Param par, double val) const {
      if (val < _clamps[0][par]) return _clamps[0][par];
      else if (val > _clamps[1][par]) return _clamps[1][par];
      else return val;
   }

   double _params[Par_Max];
   bool _fixed[Par_Max];
   double _clamps[2][Par_Max];
   size_t _phases[2][Par_Max];
   bool _isSymmetric;
   bool _isUniGap;
   bool _singleF84;
   bool _doubleF84;
   bool _isPLinkedGC;
   bool _isQLinkedGC;
   bool _isPLinkedGA;
   bool _isQLinkedGA;
};

std::ostream& operator<<(std::ostream& os, const FTAGParams& params);
std::istream& operator>>(std::istream& is, FTAGParams& params);
#endif
