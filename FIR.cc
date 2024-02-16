#include "FIR.h"
#include <math.h>


FIR::FIR()
{
  // default constructor for collections
}

FIR::FIR(const unsigned taps, const bool verbose) : 
  taps(taps), verbose(verbose), delayLine(taps), weights(taps)
{
}


void FIR::init(const unsigned taps, const bool verbose) 
{
    this->taps = taps;
    this->verbose = verbose;
    this->delayLine.resize(taps);
    this->weights.resize(taps);
}


void FIR::setWeights(float* newWeights)
{
  for (int tap = 0; tap < taps; tap++) {
    weights[tap] = newWeights[tap];
  }
}


fcomplex FIR::processNextSample(const fcomplex sample)
{
  fcomplex sum = sample * weights[0];
  delayLine[0] = sample;

  for (int tap = taps; -- tap > 0;) {
    sum += weights[tap] * delayLine[tap];
    delayLine[tap] = delayLine[tap - 1];
  }

  return sum;
}
