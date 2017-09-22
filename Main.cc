#include "FilterBank.h"
#include <stdlib.h>
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

  FilterBank fb = FilterBank(true, taps, channels, window);
  fb.printWeights();
}
