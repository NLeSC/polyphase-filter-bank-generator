POLYPHASE FILTER BANK GENERATOR
===============================

Implemented by Rob V. van Nieuwpoort, http://www.vannieuwpoort.com/,
at the Netherlands eScience center (https://www.esciencecenter.nl/).

This program is a stand-alone version of the polyphase filter bank
generator I designed and implemented for the LOFAR telescope (http://www.lofar.org/). This
code generates the filter weights for polyphase filter banks with
arbitrary numbers of channels, and with configurable windows.  Window
types currently supported are: HAMMING, BLACKMAN, GAUSSIAN, and
KAISER.  The original code is a part of the LOFAR real-time central
processor (the correlator). The code is completely generic, and can
be used for other telescopes, or even completely different signal
processing applications as well.

The paper below describes the LOFAR real-time central processing
pipeline. This production pipeline uses the filter bank generator to
generate the correct polyphase filter banks at run time, depending on
the telescope paramters.

John W. Romein, P. Chris Broekema, Jan David Mol, Rob V. van Nieuwpoort:
The LOFAR Correlator: Implementation and Performance Analysis,
ACM Symposium on Principles and Practice of Parallel Programming (PPoPP’10), Bangalore, India, pp. 169-178, January, 2010.
https://vannieuwpoort.com/wp-content/uploads/lofar.pdf

The paper describes the usage of the filter bank as follows.

The LOFAR subband data are processed by a Poly-Phase Filter bank
(PPF) that splits a frequency subband into a number of narrower
frequency channels. In this step, we trade time resolution for frequency
resolution: we split a subband into N separate channels, but
with an N-times lower sampling rate per channel. With the higher
frequency resolution, we can remove RFI artifacts with a higher accuracy
later in the pipeline. For LOFAR, typically a 195 KHz subband is split
into 256 channels of 763 Hz, but the filter supports any reasonable
power-of-two number of channels for different observation modes.
The PPF consists of two parts. First, the data are filtered using
Finite Impulse Response (FIR) filters. A FIR filter simply multiplies
a sample with a real weight factor, and also adds a number
of weighted samples from the past. Since we have to support different
numbers of channels, our software automatically designs a
filter bank with the desired properties and number of channels at
run time, generating the FIR filter weights on the fly. This again
demonstrates the flexibility of a software solution. For performance
reasons, the implementation of the filter is done in assembly. Next,
the filtered data are Fourier Transformed.

Please cite this paper if this code is useful to you.


dependencies:
-------------

This code needs FFTW3 to run. (On Debian / Ubuntu based systems, you can use "sudo apt install libfftw3-dev" to install it.
Gnuplot is used to show the output, but this is optional.

Compiling the code:
-------------------

just type "make". The code should be compiled, and you should now have an executable called "polyphase-filter-bank-generator".
You can run the code as follows: ./polyphase-filter-bank-generator [nrChannels] [nrTaps] [windowType]", 
where windowType is one of HAMMING, BLACKMAN, GAUSSIAN, KAISER.

Visualizing the output:
-----------------------

You can show the filter constants of the filter bank by running "make plot". 
This will generate a small filter bank with 32 channels, and 16 filter taps per channel, using all different window options. 
The filter constants are saved to a file called "[WINDOW]-example.data". Next, gnuplot is used to plot the data, saving the result to example.pdf.

.. image:: example.jpg?raw=true



RELATED WORK: polyphase filter bank implementations on CPUs and GPUs.
---------------------------------------------------------------------

You can use the filter bank constants generated with this generator program to create a filter bank.
I also worked on polyphase filter bank implementations for GPUs and multi-core processors.
The code runs on Intel CPUs (written in C), NVIDIA (with Cuda) and AMD GPUs (with OpenCL), and on the simulated MicroGrid architecture. 
The source code for the filter banks is available here:
https://vannieuwpoort.com/wp-content/uploads/2023/05/ppf.zip

For more information, see this paper:

Karel van der Veldt, Rob van Nieuwpoort, Ana Lucia Varbanescu and Chris Jesshope:
A Polyphase Filter For GPUs And Multi-Core Processors
First Workshop on High Performance Computing in Astronomy (AstroHPC 2012)
In conjunction with the 21-st International ACM Symposium on High-Performance Parallel and Distributed Computing (HPDC 2012) June 19, 2012, Delft, the Netherlands.
https://vannieuwpoort.com/wp-content/uploads/astro05-vanderveldt.pdf.

For more details on the implementation, you can also have a
look at Karel’s master thesis:
A Polyphase Filter For GPUs And Multi-Core Processors.
https://vannieuwpoort.com/wp-content/uploads/Karel-van-der-Veldt.pdf




For more information on Polyphase filter banks in general, please see the paper by Harris et al.:

F.J. Harris ; C. Dick ; M. Rice
Digital receivers and transmitters using polyphase filter banks for wireless communications
IEEE Transactions on Microwave Theory and Techniques ( Volume: 51, Issue: 4, Apr 2003 )
Page(s): 1395 - 1412
April 2003 
DOI: 10.1109/TMTT.2003.809176
