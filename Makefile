CC		= g++
CC_FLAGS	= -g

CC_SOURCES	= FilterBank.cc Main.cc

CC_OBJECTS	= $(CC_SOURCES:%.cc=%.o)

CC_LINK_FLAGS	= -lfftw3f

%.o:		%.cc
		$(CC) $(CC_FLAGS) -c $< -o $@


polyphase-filter-bank-generator:	$(CC_OBJECTS)
		$(CC) $^ -o $@ $(CC_LINK_FLAGS)

plot: polyphase-filter-bank-generator
	./polyphase-filter-bank-generator 32 16 KAISER > example.data
	gnuplot example.gnuplot > example.pdf

clean:
	rm -f *.o polyphase-filter-bank-generator *~ example.data example.pdf
