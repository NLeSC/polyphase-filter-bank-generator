#include "FilterBank.h"
#include <cstring>
#include <vector>

#include <math.h>
#include <iostream>
#include <cstring>
#include <fftw3.h>

#define fftw_real(x)     ((x)[0])
#define fftw_imag(x)     ((x)[1])

using namespace std;

// For documentation on this class, see the header file.

FilterBank::FilterBank(bool verbose, unsigned channels, unsigned taps, WindowType windowType)
:
  itsWindowType(windowType), itsNrChannels(channels), itsNrTaps(taps), itsVerbose(verbose)
{
  generate_filter();
}


unsigned FilterBank::getNrTaps()
{
  return itsNrTaps;
}


unsigned FilterBank::getNrChannels()
{
  return itsNrChannels;
}


float *FilterBank::getWeights(unsigned channel)
{
  return &weights[channel*itsNrTaps];
}


float *FilterBank::getWeights()
{
  return weights;
}


WindowType FilterBank::getWindowTypeFromString(char* s)
{
  if(!strcmp(s, "HAMMING")) {
    return HAMMING;
  }

  if(!strcmp(s, "BLACKMAN")) {
    return BLACKMAN;
  }

  if(!strcmp(s, "GAUSSIAN")) {
    return GAUSSIAN;
  }

  if(!strcmp(s, "KAISER")) {
    return KAISER;
  }

  cerr << "Unknown window type" << endl;
  
  return ERROR;
}


// Returns the first power of two higher than n.
unsigned FilterBank::nextPowerOfTwo(unsigned n)
{
  unsigned p;

  for (p = 1; p < n; p <<= 1)
    ;

  return p;
}


// hamming window function
void FilterBank::hamming(unsigned n, double d[])
{
  if (n == 1) {
    d[0] = 1.0;
    return;
  }

  unsigned m = n - 1;

  for (unsigned i = 0; i < n; i++) {
    d[i] = 0.54 - 0.46 * cos((2.0 * M_PI * i) / m);
  }
}


// blackman window function
void FilterBank::blackman(unsigned n, double d[])
{
  if (n == 1) {
    d[0] = 1.0;
    return;
  }

  unsigned m = n - 1;

  for (unsigned i = 0; i < n; i++) {
    double k = i / m;
    d[i] = 0.42 - 0.5 * cos(2.0 * M_PI * k) + 0.08 * cos(4.0 * M_PI * k);
  }
}


// Guassian window function
void FilterBank::gaussian(int n, double a, double d[])
{
  int index = 0;

  for (int i = -(n - 1); i <= n - 1; i += 2) {
    d[index++] = exp(-0.5 * pow((a / n * i), 2));
  }
}


// Compute the modified Bessel function I_0(x) for any real x.
// This method was taken from the ROOT package, See http://root.cern.ch/root.
// It was released undet the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1
double FilterBank::besselI0(double x)
{
  // Parameters of the polynomial approximation
  const double p1 = 1.0, p2 = 3.5156229, p3 = 3.0899424, p4 = 1.2067492, p5 = 0.2659732, p6 = 3.60768e-2, p7 = 4.5813e-3;

  const double q1 = 0.39894228, q2 = 1.328592e-2, q3 = 2.25319e-3, q4 = -1.57565e-3, q5 = 9.16281e-3, q6 = -2.057706e-2, q7 = 2.635537e-2, q8 = -1.647633e-2,
      q9 = 3.92377e-3;

  const double k1 = 3.75;
  double ax = fabs(x);

  double y = 0, result = 0;

  if (ax < k1) {
    double xx = x / k1;
    y = xx * xx;
    result = p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7)))));
  } else {
    y = k1 / ax;
    result = (exp(ax) / sqrt(ax)) * (q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))))))));
  }

  return result;
}


// Kaiser window function
void FilterBank::kaiser(int n, double beta, double d[])
{
  if (n == 1) {
    d[0] = 1.0;
    return;
  }

  int m = n - 1;

  for (int i = 0; i < n; i++) {
    double k = 2.0 * beta / m * sqrt(i * (m - i));
    d[i] = besselI0(k) / besselI0(beta);
  }
}


// One-dimensional interpolation. Interpolate Y, defined at the points X, 
// at N evenly spaced points between 0 and 1. The sample points X must be strictly monotonic
void FilterBank::interpolate(const double x[], const double y[], unsigned xlen, unsigned n, double result[])
{
  unsigned nextX = 0;
  unsigned index = 0;

  for (double interpolatedX = 0; interpolatedX <= 1.0; interpolatedX += 1.0 / (n - 1), index++) {
    while (x[nextX] <= interpolatedX && nextX < xlen - 1)
      nextX++;

    if (nextX == 0) {
      cerr << "ERROR in FilterBank::interpolate";
      return;
    }

    double prevXVal = x[nextX - 1];
    double nextXVal = x[nextX];
    double prevYVal = y[nextX - 1];
    double nextYVal = y[nextX];

    double rc = (nextYVal - prevYVal) / (nextXVal - prevXVal);

    double newVal = prevYVal + (interpolatedX - prevXVal) * rc;
    result[index] = newVal;
  }
}


// Compute the filter, similar to Octave's fir2(n, f, m, grid_n, ramp_n, window);
// Window and result must be of size n+1.
// grid_n: length of ideal frequency response function
// ramp_n: transition width for jumps in filter response
// defaults to grid_n/20; a wider ramp gives wider transitions
// but has better stopband characteristics.
void FilterBank::generate_fir_filter(unsigned n, double w, const double window[], double result[])
{
  // make sure grid is big enough for the window
  // the grid must be at least (n+1)/2
  // for all filters where the order is a power of two minus 1, grid_n = n+1;
  unsigned grid_n = nextPowerOfTwo(n + 1);

  unsigned ramp_n = grid_n/20;

  // Apply ramps to discontinuities
  // this is a low pass filter
  // maybe we can omit the "w, 0" point?
  // I did observe a small difference
  double f[] = { 0.0, w - ramp_n / grid_n / 2.0, w, w + ramp_n / grid_n / 2.0, 1.0 };
  double m[] = { 1.0, 1.0, 0.0, 0.0, 0.0 };

  // grid is a 1-D array with grid_n+1 points. Values are 1 in filter passband, 0 otherwise
  vector<double> grid(grid_n + 1);

  // interpolate between grid points
  interpolate(f, m, 5 /* length of f and m arrays */, grid_n + 1, &grid[0]);

#if 0
  cerr << "interpolated = [";
  for(unsigned i=0; i<grid_n+1; i++) {
    cerr << grid[i];
    if(i != grid_n+1-1) cerr << ", ";
  }
  cerr << "];";
  cerr << cerr.str();
#endif

  // the grid we do an ifft on is:
  // grid appended with grid_n*2 zeros
  // appended with original grid values from indices grid_n..2, i.e., the values in reverse order
  // (note, arrays start at 1 in octave!)
  // the input for the ifft is of size 4*grid_n
  // input = [grid ; zeros(grid_n*2,1) ;grid(grid_n:-1:2)];

  fftwf_complex* cinput = (fftwf_complex*) fftwf_malloc(grid_n * 4 * sizeof(fftwf_complex));
  fftwf_complex* coutput = (fftwf_complex*) fftwf_malloc(grid_n * 4 * sizeof(fftwf_complex));

  if (cinput == NULL || coutput == NULL) {
    cerr << "cannot allocate buffers";
    return;
  }

  // wipe imaginary part
  for (unsigned i = 0; i < grid_n * 4; i++) {
    fftw_imag(cinput[i]) = 0.0;
  }

  // copy first part of grid
  for (unsigned i = 0; i < grid_n + 1; i++) {
    fftw_real(cinput[i]) = grid[i];
  }

  // append zeros
  for (unsigned i = grid_n + 1; i <= grid_n * 3; i++) {
    fftw_real(cinput[i]) = 0.0;
  }

  // now append the grid in reverse order
  for (unsigned i = grid_n - 1, index = 0; i >= 1; i --, index ++) {
    fftw_real(cinput[grid_n * 3 + 1 + index]) = grid[i];
  }

#if 0
  cerr << "ifft_in = [";
  for(unsigned i=0; i<grid_n*4; i++) {
    cerr << fftw_real(cinput[i]) << " " << fftw_imag(cinput[i]);
    if(i != grid_n*4-1) cerr << ", ";
  }
  cerr << "];";
#endif

  fftwf_plan plan = fftwf_plan_dft_1d(grid_n * 4, cinput, coutput, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_execute(plan);

#if 0
  for(unsigned i=0; i<grid_n*4; i++) {
    cerr << "ifft result [" << i << "] = " << fftw_real(coutput[i]) << " " << fftw_imag(coutput[i]);
  }
#endif

  //                        half                   end
  // 1 2       n+1          2(n+1)      3(n+1)     4(n+1)
  //                                    x x x x x x x x x     # last quarter
  //   x x x x x x                                            # first quarter

  // last_quarter  = b([end-n+1:2:end]); # the size is only 1/8, since we skip half of the elements
  // first_quarter = b(2:2:(n+1));       # the size is only 1/8, since we skip half of the elements

  unsigned index = 0;

  for (unsigned i = 4 * grid_n - n; i < 4 * grid_n; i += 2) {
    result[index] = fftw_real(coutput[i]);
    index++;
  }

  for (unsigned i = 1; i <= n; i += 2) {
    result[index] = fftw_real(coutput[i]);
    index++;
  }

  fftwf_destroy_plan(plan);
  fftwf_free(cinput);
  fftwf_free(coutput);

  // multiply with window
  for (unsigned i = 0; i <= n; i++) {
    result[i] *= window[i];
  }

  // normalize
  double factor = result[n / 2];
  for (unsigned i = 0; i <= n; i++) {
    result[i] /= factor;
  }

#if 0
  cerr << "result = [";
  for(unsigned i=0; i<=n; i++) {
    cerr << result[i];
    if(i != n) cerr << ", ";
  }
  cerr << "];";
  cerr << cerr.str();
#endif
}


// This method initializes the weights array.
void FilterBank::generate_filter()
{
  unsigned n = itsNrChannels * itsNrTaps;

  if (itsVerbose) {
    cerr << "generating FIR filter bank with " << itsNrChannels << " channels and " << itsNrTaps << " taps (" << n << " total), using a ";
  }

  vector<double> window(n);

  switch (itsWindowType) {
  case HAMMING: {
    // Use a n-point Hamming window.
    if (itsVerbose) {
      cerr << "Hamming window";
    }
    hamming(n, &window[0]);
    break;
  }
  case BLACKMAN: {
    // Use a n-point Blackman window.
    if (itsVerbose) {
      cerr << "Blackman window";
    }
    blackman(n, &window[0]);
    break;
  }
  case GAUSSIAN: {
    // Use a n-point Gaussian window.
    double alpha = 3.5;
    if (itsVerbose) {
      cerr << "Gaussian window with alpha = " << alpha;
    }
    gaussian(n, alpha, &window[0]);
    break;
  }
  case KAISER: {
    // Use a n-point Kaiser window.
    // The beta parameter is found in matlab / octave with
    // [n,Wn,bta,filtype]=kaiserord([fsin/channels 1.4*fsin/channels],[1 0],[10^(0.5/20) 10^(-91/20)],fsin);
    // where fsin is the sample freq
    double beta = 9.0695;
    if (itsVerbose) {
      cerr << "Kaiser window with beta = " << beta;
    }
    kaiser(n, beta, &window[0]);
    break;
  }
  default:
    cerr << "unknown window type";
    return;
  }

  if (itsVerbose) {
    cerr << endl;
  }
  
#if 0
  cerr << "window = [";
  for(unsigned i=0; i<n; i++) {
    cerr << window[i];
    if(i != n-1) cerr << ", ";
  }
  cerr << "];";
#endif

  vector<double> result(n);

  generate_fir_filter(n - 1, 1.0 / itsNrChannels, &window[0], &result[0]);

  weights = new float[itsNrChannels * itsNrTaps];

  unsigned index = 0;
  for (int tap = 0; tap < itsNrTaps; tap++) {
    for (unsigned channel = 0; channel < itsNrChannels; channel++) {
      // Correct total power.
      // we use the 256 channel case as a reference, so we
      // multiply by 256, and divide by the number of channels
      weights[channel*itsNrTaps + tap] = result[index] * 256.0 / itsNrChannels;
      index++;
    }
  }
}


void FilterBank::reverseTaps() {
  for (unsigned channel = 0; channel < itsNrChannels; channel++) {
    for (unsigned tap = 0; tap < itsNrTaps/2; tap++) {
      float tmp = weights[channel*itsNrTaps + itsNrTaps - tap - 1];
      weights[channel * itsNrTaps + itsNrTaps - tap - 1] = weights[channel * itsNrTaps + tap];
      weights[channel * itsNrTaps + tap] = tmp;
    }
  }
}


void FilterBank::negateWeights() {
  for (int tap = itsNrTaps - 1; tap >= 0; tap--) { // store the taps in reverse!
    for (unsigned channel = 0; channel < itsNrChannels; channel++) {
      // Negate all odd channels
      if (channel % 2 != 0) {
        weights[channel*itsNrTaps+tap] = -weights[channel*itsNrTaps+tap];
      }
    }
  }
}


// Print the weights array in the natural order, in a format that can be read by gnuplot.
void FilterBank::printWeights() {
  for (int tap = 0; tap < itsNrTaps; tap++) {
    for (unsigned channel = 0; channel < itsNrChannels; channel++) {
      cout << weights[channel*itsNrTaps+tap] << endl;
    }
  }
}
