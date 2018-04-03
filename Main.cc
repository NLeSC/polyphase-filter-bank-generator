#define USE_PPF 0

#include "FilterBank.h"

#if USE_PPF
#include "PPF.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;


int main (int argc, char *argv[])
{
  unsigned channels = 4;
  unsigned taps = 16;
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
  fb.printWeights();
  fb.reverseTaps();
  
#if USE_PPF
  vector<complex<float> > response(taps);
  PPF ppf(channels, taps);
  ppf.getImpulseResponse(0, response);

  for(int i=0; i<taps; i++) {
    cerr << response[i].real() << endl;
  }
#endif

#if USE_PPF
  vector<complex<float> > response(taps);
  PPF ppf(channels, taps);
  ppf.getFrequencyResponse(0, response);

  for(int i=0; i<taps; i+=2) {
    cerr << response[i].real() << endl;
  }
#endif
}
