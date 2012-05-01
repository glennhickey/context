//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGARGS_H
#define _FTAGARGS_H

#include <ostream>
#include <string>

struct FTAGArgs
{
   FTAGArgs();
   size_t getArgs(int argc, char** argv);

   bool _relative;
   static const bool RelativeDefault;
   bool _symmetric;
   static const bool SymmetricDefault;
   bool _unigap;
   static const bool UniGapDefault;
   bool _scale;
   static const bool ScaleDefault;
   bool _context;
   static const bool ContextDefault;
   bool _pflat;
   static const bool PFlatDefault;
   bool _qflat;
   static const bool QFlatDefault;
   bool _plgc;
   static const bool PLinkedGCDefault;
   bool _plga;
   static const bool PLinkedGADefault;
   bool _qlgc;
   static const bool QLinkedGCDefault;
   bool _qlga;
   static const bool QLinkedGADefault;
   size_t _numOptTrials;
   static const size_t NumOptTrialsDefault;
   size_t _maxOptIt;
   static const size_t MaxOptItDefault;
   double _optThreshold;
   static const double OptThresholdDefault;
   double _offset;
   static const double OffsetDefault;
   bool _randOrder;
   static const bool RandOrderDefault;
   bool _keepGaps;
   static const bool KeepGapsDefault;
   bool _simFasta;
   static const bool SimFastaDefault;
   bool _singleF84;
   static const bool SingleF84Default;
   bool _doubleF84;
   static const bool DoubleF84Default;
   
   double _fixA;
   static const double FixADefault;
   double _fixB;
   static const double FixBDefault;
   double _fixC;
   static const double FixCDefault;
   double _fixD;
   static const double FixDDefault;

   double _tmBias;
   static const double TMBiasDefault;
   double _emBias;
   static const double EMBiasDefault;

   size_t _maxEMIt;
   static const size_t MaxEMItDefault;
   double _emThreshold;
   static const double EMThresholdDefault;
   size_t _emConvRepeats;
   static const size_t EMConvRepeatsDefault;
   double _emCRThreshold;
   static const double EMCRThresholdDefault;

   size_t _seed;
   std::string _fpFile;
   std::string _inFile;
   std::string _inFile2;
   std::string _outFile;
   std::string _outFile2;
   std::string _distFile;
   std::string _fpFileOut;

   size_t _numPairs;
   static const size_t NumPairsDefault;
   size_t _maxLength;
   static const size_t MaxLengthDefault;
   size_t _minLength;
   static const size_t MinLengthDefault;
   size_t _winSize;
   static const size_t WinSizeDefault;
   double _flipProb;
   static const double FlipProbDefault;

   static char* OSSep;
};

std::ostream& operator<<(std::ostream& os, const FTAGArgs& params);

#endif
