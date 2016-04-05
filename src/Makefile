CC = gcc
UPCC = upcc

KMER_LENGTH 		= 19

# Add -std=gnu99 to CFLAGS if use gnu compiler
CFLAGS 	= -O3 
UPCFLAGS = -O3 -std=gnu99
DEFINE 	= -DKMER_LENGTH=$(KMER_LENGTH)
HEADERS	= contig_generation.h kmer_hash.h packingDNAseq.h
LIBS	=

TARGETS	= serial pgen

all: 	$(TARGETS)

serial: serial.c $(HEADERS)
		$(CC) $(CFLAGS) -o $@ $< $(DEFINE) $(LIBS) -DKMER_PACKED_LENGTH=$(shell echo $$((($(KMER_LENGTH)+3)/4)))

pgen:	richard-pgen.upc richard-kmer.h
		$(UPCC) $(UPCFLAGS) -Wc,"$(CFLAGS)" -o $@ $< $(DEFINE) $(LIBS)

clean :
	rm -f *.o
	rm -rf $(TARGETS)
