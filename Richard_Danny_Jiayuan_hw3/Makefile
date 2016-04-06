CC = gcc
UPCC = upcc

KMER_LENGTH 		= 51

# Add -std=gnu99 to CFLAGS if use gnu compiler
CFLAGS 	 = -O3 
DEFINE 	 = -DKMER_LENGTH=$(KMER_LENGTH)
LIBS	 =
UPCFLAGS = -O

TARGETS	= pgen

all: 	$(TARGETS)

pgen:	pgen.upc kmer.h
		$(UPCC) $(UPCFLAGS) -Wc,"-O3" -o $@ $< $(DEFINE) $(LIBS)

clean :
	rm -f *.o
	rm -rf $(TARGETS)
