#include "FilterBank.h"
#include "PPF.h"

#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;


int main (int argc, char *argv[])
{
  unsigned channels = 4;
  unsigned taps = 256;
  unsigned nrSamplesPerIntegration = 64;

  WindowType window = KAISER;
  
  if(argc == 1) {
    ; // use default settings
  } else if(argc != 4) {
    cerr << "Usage: " << argv[0] << " [nrChannels] [nrTaps] [windowType]" << endl;
    cerr << "       where windowType is one of HAMMING, BLACKMAN, GAUSSIAN, KAISER" << endl;
    return -1;
  } else {
    channels = atoi(argv[1]);
    taps = atoi(argv[2]);
    window = FilterBank::getWindowTypeFromString(argv[3]);

    if(window == ERROR) {
      cerr << "WindowType should be one of HAMMING, BLACKMAN, GAUSSIAN, KAISER" << endl;
      return -1;
    }
  }

  FilterBank fb = FilterBank(false, channels, taps, window);
  //fb.printWeights();
  fb.reverseTaps();

  FIR fir(taps, false);
  fir.setWeights(fb.getWeights(0));
  fir.processNextSample(fcomplex(1,0));

  for (int i = 1; i < taps; i++)
  {
    fcomplex sample = fir.processNextSample(fcomplex(0,0));
    cerr  << sample.real() << endl;
  }
  

#if 0
  // Do some filtering
  // in = 1D array of size nrSamplesPerIntegration * nrChannels
  // out = [nrChannels][nrSamplesPerIntegration]

  vector<fcomplex> in(nrSamplesPerIntegration * channels);
  in[0] = fcomplex(1,0.0f);
/*
  for(int i=0; i< nrSamplesPerIntegration * channels; i++) {
    in[i] = fcomplex(i,0.0f);
    cerr << "i = " << in[i] << endl;
  }
*/
  boost::multi_array<fcomplex, 2> out(boost::extents[channels][nrSamplesPerIntegration]);

  PPF ppf(channels, taps, nrSamplesPerIntegration, false);

  ppf.filter(in, out);
  for(int i=0; i< nrSamplesPerIntegration; i++) {
    cerr << "out = " << out[0][i] << endl;
  }
#endif
#if 0
  ppf.getImpulseResponse(0, response);

  for(int i=0; i<taps; i++) {
    cerr << response[i].real() << endl;
  }
#endif

#if 0
  vector<complex<float> > response(taps);
  PPF ppf(channels, taps);
  ppf.getFrequencyResponse(0, response);

  for(int i=0; i<taps; i+=2) {
    cerr << response[i].real() << endl;
  }
#endif
}
