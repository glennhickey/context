//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <limits>
#include "ftagscore.h"

using namespace std;
using namespace FTAGStates;

static const double NEGINF = -numeric_limits<double>::max();


FTAGScore::FTAGScore(const FTAGModel* model) 
{
  if (model != NULL)
  {
    setModel(*model);
  }
  else
  {
    _model = NULL;
  }
}

FTAGScore::~FTAGScore()
{

}

void FTAGScore::setModel(const FTAGModel& model) 
{
  _model = &model;
  _transModel = &_model->getTransitionModel();
  _emModel = &_model->getEmissionModel();
}

double FTAGScore::logViterbiScore(const pair<string, string>& alignment)
{
  assert(_model != NULL && _transModel != NULL && _emModel != NULL);
  _alignment = &alignment;
  assert(_alignment->first.length() == _alignment->second.length());
  _N = _alignment->first.length();

  // init the DP table
  init2(_table, 0., 0.);
  init2(_lastTable, 0., 0.);
  init2(_lastState, S_Start, S_Start);

  for (size_t pos = 0; pos < _N; ++pos)
  {
    _state[0] = getColID(pos);
    _state[1] = getContext(_state[0]);
    if (_state[0] == S_Mx)
    {
      scoreMatch(pos);
    }
    else
    {
      scoreIndel(pos);
      // skip ahead
      pos += _indelSize - 1; 
    }
    copy2(_lastState, _state);
    copy2(_lastTable, _table);
  }

  _table[0] = _lastTable[0] + log(_transModel->prTrans(_lastState[0], S_End));
  _table[1] = _lastTable[1] + log(_transModel->prTrans(_lastState[1], S_End));
  return std::max(_table[0], _table[1]);
}

void FTAGScore::scoreMatch(size_t pos)
{
  DNA c1 = DNASubModel::charToDNA(_alignment->first[pos]);
  DNA c2 = DNASubModel::charToDNA(_alignment->second[pos]);
  // emission probability
  double emProb = log(_emModel->prFullM(c1, c2));
  
  // DP step
  double tPr0 = _lastTable[0] + log(_transModel->prTrans(_lastState[0], S_Mx));
  double tPr1 = tPr0;
  if (_lastState[0] != S_Mx)
  { 
    tPr1 = _lastTable[1] + log(_transModel->prTrans(_lastState[1], S_Mx));
  }
  _table[0] = emProb + max(tPr0, tPr1);
  _table[1] = _table[0];
}

void FTAGScore::scoreIndel(size_t pos)
{
  // compute score of vanilla indel
  double insScore = 0.0;
  _indelSize = 0;
  for (size_t i = pos; i < _N && getColID(i) == _state[0]; ++i, ++_indelSize)
  {
    // compute emission log prob
    if (_state[0] == S_Ix)
    {
      DNA c = DNASubModel::charToDNA(_alignment->second[i]);
      insScore += log(_emModel->prFullI(c));
    }
    else
    {
      assert(_state[0] == S_Dx);
      DNA c = DNASubModel::charToDNA(_alignment->first[i]);
      insScore += log(_emModel->prFullD(c));
    } 
    // add transition log prob
    if (i > pos)
    {
      insScore += log(_transModel->prTrans(_state[0], _state[0]));
    }
  }
  // DP step.
  double tPr0 = _lastTable[0] + log(_transModel->prTrans(_lastState[0], 
                                                         _state[0]));
  double tPr1 = _lastTable[1] + log(_transModel->prTrans(_lastState[1],
                                                         _state[0]));
  _table[0] = insScore + max(tPr0, tPr1);


  // compute score of contexst indel;
  _table[1] = NEGINF;
  if (pos >= _indelSize)
  {
    State inside = _state[0] == S_MIx ? S_MIzM : S_MDyM;
    insScore = 0.0;
    for (size_t i = pos; i < _N && getColID(i) == _state[0]; ++i)
    {
      DNA c1 = DNASubModel::charToDNA(_alignment->first[i - _indelSize]);
      DNA c2 = DNASubModel::charToDNA(_alignment->second[i - _indelSize]);
      // compute emission log prob
      if (_state[0] == S_Ix)
      {
        DNA c3 = DNASubModel::charToDNA(_alignment->second[i]);
        insScore += log(_emModel->prFullMI(c1, c2, c3));
      }
      else
      {
        DNA c3 = DNASubModel::charToDNA(_alignment->first[i]);
        insScore += log(_emModel->prFullMD(c1, c2, c3));        
      }
      // compute transition log prob
      if (i == pos + 1)
      {
        insScore += log(_transModel->prTrans(_state[1], inside));
      }
      else if (i > pos + 1)
      {
        insScore += log(_transModel->prTrans(inside, inside));
      }
    }
    // close the context indel
    State eState = _indelSize > 1 ? inside : _state[1];
    insScore += log(_transModel->prTrans(eState, S_Close));

    // DP step
    double tPr0 = _lastTable[0] + log(_transModel->prTrans(_lastState[0], 
                                                           _state[1]));
    double tPr1 = _lastTable[1] + log(_transModel->prTrans(_lastState[1], 
                                                           _state[1]));
    _table[1] = insScore + max(tPr0, tPr1);
  }
}

State FTAGScore::getColID(size_t pos)
{
  bool topGap = _alignment->first[pos] == '-';
  bool botGap = _alignment->second[pos] == '-';
  assert(!topGap || !botGap);
  if (topGap == false and botGap == false)
  {
    return S_Mx;
  }
  else if (topGap == true)
  {
    return S_Ix;
  }
  else
  {
    return S_Dx;
  }
}

