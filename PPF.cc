#include "PPF.h"

#include <complex>
#include <cmath>
#include <stdexcept>


PPF::PPF(unsigned nrChannels, unsigned nrTaps, unsigned nrSamplesPerIntegration, bool verbose)
:
  nrSamplesPerIntegration(nrSamplesPerIntegration),
  nrChannels(nrChannels),
  nrTaps(nrTaps),
  filterBank(verbose, nrTaps, nrChannels, KAISER),
  FIRs(nrChannels),
  FFTinData(boost::extents[nrTaps - 1 + nrSamplesPerIntegration][nrChannels])
{
  init_fft();

  // In CEP, the first subband is from -98 KHz to 98 KHz, rather than from 0 to 195 KHz.
  // To avoid that the FFT outputs the channels in the wrong order (from 128 to
  // 255 followed by channels 0 to 127), we multiply each second FFT input by -1.
  // This is efficiently achieved by negating the FIR filter constants of all
  // uneven FIR filters.
  filterBank.negateWeights();

  // Init the FIR filters themselves with the weights of the filterbank.
  for (unsigned chan = 0; chan < nrChannels; chan ++) {
    FIRs[chan].init(nrTaps, verbose);
    FIRs[chan].setWeights(filterBank.getWeights(chan));
  }
}


PPF::~PPF()
{
  destroy_fft();
}


void PPF::init_fft()
{
  fftwf_complex *buf = static_cast<fftwf_complex *>(fftwf_malloc(2 * nrChannels * sizeof(fftwf_complex)));
  assert(buf);
     
  FFTWPlan = fftwf_plan_dft_1d(nrChannels, buf, buf + nrChannels, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_free(buf);
}


void PPF::destroy_fft()
{
  fftwf_destroy_plan(FFTWPlan);
}


void PPF::filter(const std::vector<fcomplex> &in, boost::multi_array<fcomplex, 2> &out)
{
  for (int chan = 0; chan < (int) nrChannels; chan ++) {
    for (unsigned time = 0; time < nrTaps - 1 + nrSamplesPerIntegration; time ++) {
	    fcomplex sample = in[chan*nrSamplesPerIntegration + time];
	    FFTinData[time][chan] = FIRs[chan].processNextSample(sample);
    }
  }

  std::vector<fcomplex> fftOutData(nrChannels);

  for (int time = 0; time < nrSamplesPerIntegration; time++) {
    fftwf_execute_dft(FFTWPlan,
    (fftwf_complex *) FFTinData[nrTaps - 1 + time].origin(),
    (fftwf_complex *) (void *) &fftOutData[0]);

  	for (unsigned chan = 0; chan < nrChannels; chan++) {
     out[chan][time] = fftOutData[chan];
    }
  }
}
