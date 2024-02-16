CC		= g++
CC_FLAGS	= -g

CC_SOURCES	= FilterBank.cc FIR.cc PPF.cc Main.cc

CC_OBJECTS	= $(CC_SOURCES:%.cc=%.o)

CC_LINK_FLAGS	= -lfftw3f

%.o:		%.cc
		$(CC) $(CC_FLAGS) -c $< -o $@


polyphase-filter-bank-generator:	$(CC_OBJECTS)
		$(CC) $^ -o $@ $(CC_LINK_FLAGS)

plot: polyphase-filter-bank-generator
	./polyphase-filter-bank-generator 32 16 HAMMING > HAMMING-example.data
	./polyphase-filter-bank-generator 32 16 BLACKMAN > BLACKMAN-example.data
	./polyphase-filter-bank-generator 32 16 GAUSSIAN > GAUSSIAN-example.data
	./polyphase-filter-bank-generator 32 16 KAISER > KAISER-example.data
	gnuplot example.gnuplot > example.pdf

clean:
	rm -f *.o polyphase-filter-bank-generator *~ HAMMING-example.data BLACKMAN-example.data GAUSSIAN-example.data KAISER-example.data example.pdf
