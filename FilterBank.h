#ifndef FILTER_BANK_H
#define FILTER_BANK_H

enum WindowType { HAMMING, BLACKMAN, GAUSSIAN, KAISER, ERROR };

class FilterBank
{
  public:

  // This constructor designs a new filter with the specified parameters, and initializes the weights array.
  FilterBank(bool verbose, unsigned channels, unsigned taps, WindowType windowType);

  unsigned getNrTaps();
  unsigned getNrChannels();
  
  float *getWeights(unsigned channel); // returns taps weights
  float* getWeights(); // returns all weights: [nrChannels][taps]

  // Reverse the order of the taps. This can be useful, because FIR implementatons might run over the taps backwards for efficiency.
  // This is the case for the LOFAR telescope, for example.
  void reverseTaps();


  // In CEP, the first subband is from -98 KHz to 98 KHz, rather than from 0 to 195 KHz.
  // To avoid that the FFT outputs the channels in the wrong order (from 128 to
  // 255 followed by channels 0 to 127), we multiply each second FFT input by -1.
  // This is efficiently achieved by negating the FIR filter constants of all
  // uneven FIR filters.
  void negateWeights();
  
  // Print the weights array in the natural order, in a format that can be read by gnuplot.
  void printWeights();

  static WindowType getWindowTypeFromString(char* s);

private:
  // Hamming window function
  void hamming(unsigned n, double d[]);

  // Blackman window function
  void blackman(unsigned n, double d[]);

  // Gaussian window function
  void gaussian(int n, double a, double d[]);

  // Kaiser window function
  void kaiser(int n, double beta, double d[]);

  // helper functions
  double besselI0(double x);
  void interpolate(const double x[], const double y[], unsigned xlen, unsigned n, double result[]);
  void generate_fir_filter(unsigned n, double w, const double window[], double result[]);
  void generate_filter();

  // Returns the first power of two higher than n.
  unsigned nextPowerOfTwo(unsigned n);

  // The window used for generating the filter, default is KAISER.
  WindowType itsWindowType;

  const unsigned itsNrTaps;
  const unsigned itsNrChannels;
  const bool itsVerbose;

  float* weights; // [nrChannels][taps]
};

#endif // FILTER_BANK_H
