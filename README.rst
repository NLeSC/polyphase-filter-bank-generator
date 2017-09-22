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
processor (the correlator). The codes is completely generic, and can
be used for other telescopes, or even completely different signal
processing applications as well.

dependencies:

This code needs FFTW3 to run. (On Debian / Ubuntu based systems, you can use "sudo apt install libfftw3-dev" to install it.
Gnuplot is used to show the output, but this is optional.

Compiling the code:

just type "make". The code should be compiled, and you should now have an executable called "polyphase-filter-bank-generator".
You can run the code as follows: ./polyphase-filter-bank-generator [nrChannels] [nrTaps] [windowType]", 
where windowType is one of HAMMING, BLACKMAN, GAUSSIAN, KAISER.

Visualizing the output:

You can show the filter constants of the filter bank by running "make plot". 
This will generate a small filter bank with 32 channels, and 16 filter taps per channel, using a KAISER window. 
The filter constants are saved to a file called "example.data". Next, gnuplot is used to plot the data, saving the result to example.pdf.

.. image:: example.jpg?raw=true

RELATED WORK: polyphase filter bank implementations on CPUs and GPUs.
---------------------------------------------------------------------

You can use the filter bank constants generated with this generator program to create a filter bank.
I also worked on polyphase filter bank implementations for GPUs and multi-core processors.
The code runs on Intel CPUs (written in C), NVIDIA (with Cuda) and AMD GPUs (with OpenCL), and on the simulated MicroGrid architecture. 
The source code for the filter banks is available here:
http://rvannieuwpoort.synology.me/software/ppf.tgz

For more information, see this paper:

Karel van der Veldt, Rob van Nieuwpoort, Ana Lucia Varbanescu and Chris Jesshope:
A Polyphase Filter For GPUs And Multi-Core Processors
First Workshop on High Performance Computing in Astronomy (AstroHPC 2012)
In conjunction with the 21-st International ACM Symposium on High-Performance Parallel and Distributed Computing (HPDC 2012) June 19, 2012, Delft, the Netherlands.
http://rvannieuwpoort.synology.me/papers/astro05-vanderveldt.pdf.

For more details on the implementation, you can also have a
look at Karel’s master thesis:
A Polyphase Filter For GPUs And Multi-Core Processors.
http://rvannieuwpoort.synology.me/masters-theses/Karel-van-der-Veldt.pdf




For more information on Polyphase filter banks in general, please see the paper by Harris et al.:

F.J. Harris ; C. Dick ; M. Rice
Digital receivers and transmitters using polyphase filter banks for wireless communications
IEEE Transactions on Microwave Theory and Techniques ( Volume: 51, Issue: 4, Apr 2003 )
Page(s): 1395 - 1412
April 2003 
DOI: 10.1109/TMTT.2003.809176
