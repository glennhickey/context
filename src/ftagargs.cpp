//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>

#include "ftagargs.h"

using namespace std;

const bool FTAGArgs::RelativeDefault = true;
const bool FTAGArgs::SymmetricDefault = true;
const bool FTAGArgs::UniGapDefault = false;
const bool FTAGArgs::ScaleDefault = true;
const bool FTAGArgs::ContextDefault = true;
const bool FTAGArgs::PFlatDefault = false;
const bool FTAGArgs::QFlatDefault = false;
const bool FTAGArgs::PLinkedGCDefault = false;
const bool FTAGArgs::PLinkedGADefault = false;
const bool FTAGArgs::QLinkedGCDefault = false;
const bool FTAGArgs::QLinkedGADefault = false;
const size_t FTAGArgs::NumOptTrialsDefault = 1;
const size_t FTAGArgs::MaxOptItDefault = 1000;
const double FTAGArgs::OptThresholdDefault = 1e-12;
const double FTAGArgs::OffsetDefault = 0.05;
const bool FTAGArgs::RandOrderDefault = false;
const bool FTAGArgs::KeepGapsDefault = true;
const bool FTAGArgs::SimFastaDefault = false;
const bool FTAGArgs::SingleF84Default = true;
const bool FTAGArgs::DoubleF84Default = true;

const double FTAGArgs::FixADefault = -1;
const double FTAGArgs::FixBDefault = -1;
const double FTAGArgs::FixCDefault = -1;
const double FTAGArgs::FixDDefault = -1;

const double FTAGArgs::TMBiasDefault = 0.;
const double FTAGArgs::EMBiasDefault = 0.;

const size_t FTAGArgs::MaxEMItDefault = 15;
const double FTAGArgs::EMThresholdDefault = -100;
const size_t FTAGArgs::EMConvRepeatsDefault = 3;
const double FTAGArgs::EMCRThresholdDefault = 1e-6;

const size_t FTAGArgs::NumPairsDefault = 25;
const size_t FTAGArgs::MaxLengthDefault = 25;
const size_t FTAGArgs::MinLengthDefault = 15;
const size_t FTAGArgs::WinSizeDefault = 7;
const double FTAGArgs::FlipProbDefault = 0.;

char* FTAGArgs::OSSep = " ";

FTAGArgs::FTAGArgs() :
   _relative(RelativeDefault),
   _symmetric(SymmetricDefault),
   _unigap(UniGapDefault),
   _scale(ScaleDefault),
   _context(ContextDefault),
   _pflat(PFlatDefault),
   _qflat(QFlatDefault),
   _plgc(PLinkedGCDefault),
   _plga(PLinkedGADefault),
   _qlgc(QLinkedGCDefault),
   _qlga(QLinkedGADefault),
   _numOptTrials(NumOptTrialsDefault),
   _maxOptIt(MaxOptItDefault),
   _optThreshold(OptThresholdDefault),
   _offset(OffsetDefault),
   _randOrder(RandOrderDefault),
   _keepGaps(KeepGapsDefault),
   _simFasta(SimFastaDefault),
   _singleF84(SingleF84Default),
   _doubleF84(DoubleF84Default),
   _fixA(FixADefault),
   _fixB(FixBDefault),
   _fixC(FixCDefault),
   _fixD(FixDDefault),
   _tmBias(TMBiasDefault),
   _emBias(EMBiasDefault),
   _maxEMIt(MaxEMItDefault),
   _emThreshold(EMThresholdDefault),
   _emConvRepeats(EMConvRepeatsDefault),
   _emCRThreshold(EMCRThresholdDefault),
   _numPairs(NumPairsDefault),
   _maxLength(MaxLengthDefault),
   _winSize(WinSizeDefault),
   _flipProb(FlipProbDefault)
{
   _seed = time(0);
}

size_t FTAGArgs::getArgs(int argc, char** argv)
{
   size_t count = 0;
   for (int i = 1; i < argc; ++i)
   {
      if (strncmp(argv[i], "relative=", strlen("relative=")) == 0)
      {
         if (strcmp(argv[i] + strlen("relative="), "true") == 0)
         {
            _relative = true;
         }
         else if (strcmp(argv[i] + strlen("relative="), "false") == 0)
         {
            _relative = false;
         }
         else throw string("bad value for relative=");
         ++count;
      }
      else if (strncmp(argv[i], "sym=", strlen("sym=")) == 0)
      {
         if (strcmp(argv[i] + strlen("sym="), "true") == 0)
         {
            _symmetric = true;
         }
         else if (strcmp(argv[i] + strlen("sym="), "false") == 0)
         {
            _symmetric = false;
         }
         else throw string("bad value for sym=");
         ++count;
      }
      else if (strncmp(argv[i], "unigap=", strlen("unigap=")) == 0)
      {
         if (strcmp(argv[i] + strlen("unigap="), "true") == 0)
         {
            _unigap = true;
         }
         else if (strcmp(argv[i] + strlen("unigap="), "false") == 0)
         {
            _unigap = false;
         }
         else throw string("bad value for unigap=");
         ++count;
      }
      else if (strncmp(argv[i], "scale=", strlen("scale=")) == 0)
      {
         if (strcmp(argv[i] + strlen("scale="), "true") == 0)
         {
            _scale = true;
         }
         else if (strcmp(argv[i] + strlen("scale="), "false") == 0)
         {
            _scale = false;
         }
         else throw string("bad value for scale=");
         ++count;
      }
      else if (strncmp(argv[i], "context=", strlen("context=")) == 0)
      {
         if (strcmp(argv[i] + strlen("context="), "true") == 0)
         {
            _context = true;
         }
         else if (strcmp(argv[i] + strlen("context="), "false") == 0)
         {
            _context = false;
         }
         else throw string("bad value for context=");
         ++count;
      }
      else if (strncmp(argv[i], "pflat=", strlen("pflat=")) == 0)
      {
         if (strcmp(argv[i] + strlen("pflat="), "true") == 0)
         {
            _pflat = true;
         }
         else if (strcmp(argv[i] + strlen("pflat="), "false") == 0)
         {
            _pflat = false;
         }
         else throw string("bad value for pflat=");
         ++count;
      }
      else if (strncmp(argv[i], "qflat=", strlen("qflat=")) == 0)
      {
         if (strcmp(argv[i] + strlen("qflat="), "true") == 0)
         {
            _qflat = true;
         }
         else if (strcmp(argv[i] + strlen("qflat="), "false") == 0)
         {
            _qflat = false;
         }
         else throw string("bad value for qflat=");
         ++count;
      }
      else if (strncmp(argv[i], "plgc=", strlen("plgc=")) == 0)
      {
         if (strcmp(argv[i] + strlen("plgc="), "true") == 0)
         {
            _plgc = true;
         }
         else if (strcmp(argv[i] + strlen("plgc="), "false") == 0)
         {
            _plgc = false;
         }
         else throw string("bad value for plgc=");
         ++count;
      }
      else if (strncmp(argv[i], "plga=", strlen("plga=")) == 0)
      {
         if (strcmp(argv[i] + strlen("plga="), "true") == 0)
         {
            _plga = true;
         }
         else if (strcmp(argv[i] + strlen("plga="), "false") == 0)
         {
            _plga = false;
         }
         else throw string("bad value for plga=");
         ++count;
      }
      else if (strncmp(argv[i], "qlgc=", strlen("qlgc=")) == 0)
      {
         if (strcmp(argv[i] + strlen("qlgc="), "true") == 0)
         {
            _qlgc = true;
         }
         else if (strcmp(argv[i] + strlen("qlgc="), "false") == 0)
         {
            _qlgc = false;
         }
         else throw string("bad value for qlgc=");
         ++count;
      }
      else if (strncmp(argv[i], "qlga=", strlen("qlga=")) == 0)
      {
         if (strcmp(argv[i] + strlen("qlga="), "true") == 0)
         {
            _qlga = true;
         }
         else if (strcmp(argv[i] + strlen("qlga="), "false") == 0)
         {
            _qlga = false;
         }
         else throw string("bad value for qlga=");
         ++count;
      }
      else if (strncmp(argv[i], "optruns=", strlen("optruns=")) == 0)
      {
         _numOptTrials = atoi(argv[i] + strlen("optruns="));
         ++count;
      }
      else if (strncmp(argv[i], "optits=", strlen("optits=")) == 0)
      {
         _maxOptIt = atoi(argv[i] + strlen("optits="));
         ++count;
      }
      else if (strncmp(argv[i], "optthreshold=", strlen("optthreshold=")) == 0)
      {
         _optThreshold = atof(argv[i] + strlen("optthreshold="));
         ++count;
      }
      else if (strncmp(argv[i], "offset=", strlen("offset=")) == 0)
      {
         _offset = atof(argv[i] + strlen("offset="));
         ++count;
      }
      else if (strncmp(argv[i], "rorder=", strlen("rorder=")) == 0)
      {
         if (strcmp(argv[i] + strlen("rorder="), "true") == 0)
         {
            _randOrder = true;
         }
         else if (strcmp(argv[i] + strlen("rorder="), "false") == 0)
         {
            _randOrder = false;
         }
         else throw string("bad value for rorder=");
         ++count;
      }
      else if (strncmp(argv[i], "keepgaps=", strlen("keepgaps=")) == 0)
      {
         if (strcmp(argv[i] + strlen("keepgaps="), "true") == 0)
         {
            _keepGaps = true;
         }
         else if (strcmp(argv[i] + strlen("keepgaps="), "false") == 0)
         {
            _keepGaps = false;
         }
         else throw string("bad value for keepgaps=");
         ++count;
      }
      else if (strncmp(argv[i], "simfasta=", strlen("simfasta=")) == 0)
      {
         if (strcmp(argv[i] + strlen("simfasta="), "true") == 0)
         {
            _simFasta = true;
         }
         else if (strcmp(argv[i] + strlen("simfasta="), "false") == 0)
         {
            _simFasta = false;
         }
         else throw string("bad value for simfasta=");
         ++count;
      }
      else if (strncmp(argv[i], "sf84=", strlen("sf84=")) == 0)
      {
         if (strcmp(argv[i] + strlen("sf84="), "true") == 0)
         {
            _singleF84 = true;
         }
         else if (strcmp(argv[i] + strlen("sf84="), "false") == 0)
         {
            _singleF84 = false;
         }
         else throw string("bad value for sf84=");
         ++count;
      }
      else if (strncmp(argv[i], "df84=", strlen("df84=")) == 0)
      {
         if (strcmp(argv[i] + strlen("df84="), "true") == 0)
         {
            _doubleF84 = true;
         }
         else if (strcmp(argv[i] + strlen("df84="), "false") == 0)
         {
            _doubleF84 = false;
         }
         else throw string("bad value for df84=");
         ++count;
      }
      else if (strncmp(argv[i], "fixa=", strlen("fixa=")) == 0)
      {
         _fixA = atof(argv[i] + strlen("fixa="));
         ++count;
      }
      else if (strncmp(argv[i], "fixb=", strlen("fixb=")) == 0)
      {
         _fixB = atof(argv[i] + strlen("fixb="));
         ++count;
      }
      else if (strncmp(argv[i], "fixc=", strlen("fixc=")) == 0)
      {
         _fixC = atof(argv[i] + strlen("fixc="));
         ++count;
      }
      else if (strncmp(argv[i], "fixd=", strlen("fixd=")) == 0)
      {
         _fixD = atof(argv[i] + strlen("fixd="));
         ++count;
      }
      else if (strncmp(argv[i], "tmbias=", strlen("tmbias=")) == 0)
      {
         _tmBias = atof(argv[i] + strlen("tmbias="));
         ++count;
      }
      else if (strncmp(argv[i], "embias=", strlen("embias=")) == 0)
      {
         _emBias = atof(argv[i] + strlen("embias="));
         ++count;
      }
      else if (strncmp(argv[i], "emits=", strlen("emits=")) == 0)
      {
         _maxEMIt = atoi(argv[i] + strlen("emits="));
         ++count;
      }
      else if (strncmp(argv[i], "emthreshold=", strlen("emthreshold=")) == 0)
      {
         _emThreshold = atof(argv[i] + strlen("emthreshold="));
         ++count;
      }
      else if (strncmp(argv[i], "ecr=", strlen("ecr=")) == 0)
      {
         _emConvRepeats = atoi(argv[i] + strlen("ecr="));
         ++count;
      }
      else if (strncmp(argv[i], "ecrthreshold=", strlen("ecrthreshold=")) == 0)
      {
         _emCRThreshold = atof(argv[i] + strlen("ecrthreshold="));
         ++count;
      }
      else if (strncmp(argv[i], "seed=", strlen("seed=")) == 0)
      {
         _seed = atoi(argv[i] + strlen("seed="));
         ++count;
      }
      else if (strncmp(argv[i], "fpfile=", strlen("fpfile=")) == 0)
      {
         _fpFile = argv[i] + strlen("fpfile=");
         ++count;
      }
      else if (strncmp(argv[i], "fpout=", strlen("fpout=")) == 0)
      {
         _fpFileOut = argv[i] + strlen("fpout=");
         ++count;
      }
      else if (strncmp(argv[i], "infile=", strlen("infile=")) == 0)
      {
         _inFile = argv[i] + strlen("infile=");
         ++count;
      }
      else if (strncmp(argv[i], "infile2=", strlen("infile2=")) == 0)
      {
         _inFile2 = argv[i] + strlen("infile2=");
         ++count;
      }
      else if (strncmp(argv[i], "outfile=", strlen("outfile=")) == 0)
      {
         _outFile = argv[i] + strlen("outfile=");
         ++count;
      }
      else if (strncmp(argv[i], "outfile2=", strlen("outfile2=")) == 0)
      {
         _outFile2 = argv[i] + strlen("outfile2=");
         ++count;
      }
      else if (strncmp(argv[i], "n=", strlen("n=")) == 0)
      {
         _numPairs = atoi(argv[i] + strlen("n="));
         ++count;
      }
      else if (strncmp(argv[i], "maxlen=", strlen("maxlen=")) == 0)
      {
         _maxLength = atoi(argv[i] + strlen("maxlen="));
         ++count;
      }
      else if (strncmp(argv[i], "minlen=", strlen("minlen=")) == 0)
      {
         _minLength = atoi(argv[i] + strlen("minlen="));
         ++count;
      }
      else if (strncmp(argv[i], "win=", strlen("win=")) == 0)
      {
         _winSize = atoi(argv[i] + strlen("win="));
         ++count;
      }
      else if (strncmp(argv[i], "fprob=", strlen("fprob=")) == 0)
      {
         _flipProb = atof(argv[i] + strlen("fprob="));
         ++count;
      }
      else if (strncmp(argv[i], "distfile=", strlen("distfile=")) == 0)
      {
         _distFile = argv[i] + strlen("distfile=");
         ++count;
      }
      else
      {
         cout << "warning: unrecognized parameter: " << argv[i] << endl;
      }
   }
   return count;
}
                               
ostream& operator<<(ostream& os, const FTAGArgs& params)
{
   os << "relative=" << (int)params._relative << FTAGArgs::OSSep
      << "sym=" << (int)params._symmetric << FTAGArgs::OSSep
      << "unigap=" << (int)params._unigap << FTAGArgs::OSSep
      << "scale=" << (int)params._scale << FTAGArgs::OSSep
      << "context=" << (int)params._context << FTAGArgs::OSSep
      << "pflat=" << (int)params._pflat << FTAGArgs::OSSep
      << "qflat=" << (int)params._qflat << FTAGArgs::OSSep
      << "plgc=" << (int)params._plgc << FTAGArgs::OSSep
      << "plga=" << (int)params._plga << FTAGArgs::OSSep
      << "qlgc=" << (int)params._qlgc << FTAGArgs::OSSep
      << "qlga=" << (int)params._qlga << FTAGArgs::OSSep
      << "optruns=" << params._numOptTrials << FTAGArgs::OSSep
      << "optits=" << params._maxOptIt << FTAGArgs::OSSep
      << "optthrehold=" << params._optThreshold << FTAGArgs::OSSep
      << "offset=" << params._offset << FTAGArgs::OSSep
      << "rorder=" << params._randOrder << FTAGArgs::OSSep
      << "keepgaps=" << params._keepGaps << FTAGArgs::OSSep
      << "simfasta=" << params._simFasta << FTAGArgs::OSSep
      << "sf84=" << params._singleF84 << FTAGArgs::OSSep
      << "df84=" << params._doubleF84 << FTAGArgs::OSSep
      << "fixa=" << params._fixA << FTAGArgs::OSSep
      << "fixb=" << params._fixB << FTAGArgs::OSSep
      << "fixc=" << params._fixC << FTAGArgs::OSSep
      << "fixd=" << params._fixD << FTAGArgs::OSSep
      << "tmbias=" << params._tmBias << FTAGArgs::OSSep
      << "embias=" << params._emBias << FTAGArgs::OSSep
      << "emits=" << params._maxEMIt << FTAGArgs::OSSep
      << "emthreshold=" << params._emThreshold << FTAGArgs::OSSep
      << "ecr=" << params._emConvRepeats << FTAGArgs::OSSep
      << "ecrthreshold=" << params._emCRThreshold << FTAGArgs::OSSep
      << "seed=" << params._seed << FTAGArgs::OSSep
      << "fpfile=" << params._fpFile << FTAGArgs::OSSep
      << "fpout=" << params._fpFileOut << FTAGArgs::OSSep
      << "infile=" << params._inFile << FTAGArgs::OSSep
      << "infile2=" << params._inFile2 << FTAGArgs::OSSep
      << "outfile=" << params._outFile << FTAGArgs::OSSep
      << "outfile2=" << params._outFile2 << FTAGArgs::OSSep
      << "n=" << params._numPairs << FTAGArgs::OSSep
      << "maxlen=" << params._maxLength << FTAGArgs::OSSep
      << "minlen=" << params._minLength << FTAGArgs::OSSep
      << "win=" << params._winSize << FTAGArgs::OSSep
      << "fprob=" << params._flipProb << FTAGArgs::OSSep
      << "distfile=" << params._distFile;
   return os; 
}


