#ifndef FIR_H
#define FIR_H

#include <complex>
#include <vector>
#include "boost/multi_array.hpp"

typedef std::complex<float> fcomplex;

class FIR {
  public:

    // We need a default constructor, since we create vectors of FIR filters. In that case, we have to call init.
    FIR();
    void init(const unsigned taps, const bool verbose);

    FIR(const unsigned taps, const bool verbose);
    void setWeights(float* weights); // Weights must be an array with size [taps], where the taps are in reverse order.

    fcomplex processNextSample(const fcomplex sample);

private:
    unsigned taps;
    bool verbose;

    std::vector<fcomplex> delayLine; // size is nr taps

    std::vector<float> weights; // size is nr taps  
};

#endif
