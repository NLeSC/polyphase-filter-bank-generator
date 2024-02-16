#ifndef PPF_H
#define PPF_H

#include "FilterBank.h"
#include "FIR.h"

#include <boost/multi_array.hpp>
#include <fftw3.h>


class PPF
{
  public:
    PPF(unsigned nrChannels, unsigned nrTaps, unsigned nrSamplesPerIntegration, bool verbose);
    ~PPF();

   // in = 1D array of size nrSamplesPerIntegration * nrChannels
   // out = [nrChannels][nrSamplesPerIntegration]
   void filter(const std::vector<fcomplex> &in, boost::multi_array<fcomplex, 2> &out); 

  private:
    void init_fft(), destroy_fft();
 
    const unsigned nrSamplesPerIntegration;
    const unsigned nrChannels;
    const unsigned nrTaps;
    FilterBank     filterBank;
    std::vector<FIR> FIRs; // [nrChannels]

    boost::multi_array<fcomplex, 2> FFTinData; //[nrTaps - 1 + nrSamplesPerIntegration][itsNrChannels]

    fftwf_plan FFTWPlan;
};

#endif // PPF_H
