#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
//#include <upc_io.h>

//#include "packingDNAseq.h"
//#include "kmer_hash.h"

#define XSTR(s) STR(s)
#define STR(s) #s

#define KMER_PACKED_LENGTH (KMER_LENGTH/4+1)

shared int kmer_len;
shared int line_len;
shared int line_count;
int hashtable_size;

typedef struct {
  char kmer[KMER_PACKED_LENGTH];
  char l_ext;
  char r_ext;
} kmer;

typedef unsigned char ksym_t;

shared kmer *kmers;

char SymbolToBibit(char symbol){
  switch(symbol){
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
  }
  printf("Unknown character: %c\n",symbol);
  return -1;
}

void PrintCharBits(ksym_t a){
  for(int i=7;i>=0;i--)
    printf("%d",(int)((a&(1<<i))!=0));
  printf(" ");
}

void PrintPackedKmer(const ksym_t *kmer){
  for(int i=0;i<KMER_PACKED_LENGTH;i++)
    PrintCharBits(kmer[i]);
  printf("\n");
}

void PrintKmer(const kmer km){
  PrintPackedKmer(km.kmer);
  printf("\t%c %c\n",km.l_ext,km.r_ext);
}

void PackSequence(const unsigned char *unpacked, unsigned char *packed){
  char pack_offset = 0;
  for(int i=0;i<KMER_LENGTH;i++,unpacked++){
    if(pack_offset==0) //Make sure memory is initialized to 0
      *packed = 0;
    char bibit   = SymbolToBibit(*unpacked);
    *packed     |= (bibit<<pack_offset);
    pack_offset += 2;
    if(pack_offset==8){
      pack_offset = 0;
      packed++;
    }
  }
}

int CompareKmer(const ksym_t *seq1, const ksym_t *seq2){
  return memcmp(seq1, seq2, KMER_PACKED_LENGTH);
}

int64_t HashKmer(int64_t hashtable_size, const ksym_t *kmer){
  int64_t hashval;
  hashval = 5381;
  for(int i = 0; i < (KMER_PACKED_LENGTH); i++)
    hashval = kmer[i] + (hashval << 5) + hashval;

  return hashval % hashtable_size;
}

void AddKmer(const ksym_t *raw_kmer, ksym_t l_ext, ksym_t r_ext){
  kmer temp;
  //printf("Packing: %.19s\n",raw_kmer);  
  PackSequence(raw_kmer,temp.kmer);
  //PrintPackedKmer(kmer_packed);
  int64_t i = HashKmer(hashtable_size,temp.kmer);
  temp.l_ext = l_ext;
  temp.r_ext = r_ext;

  //Linear probing
  for(;kmers[i].l_ext!=0;i++){}

  //printf("Store at: %lld\n",i);

  kmers[i] = temp;
}


int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	/** Read input **/
	upc_barrier;
	//inputTime -= gettime();

  if(argc!=2){
    printf("Need to specify program file!\n");
    upc_global_exit(-4);
  }

  if(MYTHREAD==0){
    FILE *fin = fopen(argv[1],"r");
    if(!fin){
      printf("Failed to open kmers file!\n");
      upc_global_exit(-3);
    }
    char first_chars[200];
    if(fread(first_chars,sizeof(char),200,fin)!=200){
      printf("Failed to read first part of kmers file!\n");
      upc_global_exit(-2);
    }
    for(int i=0;i<200;i++){
      if(first_chars[i]=='\n'){
        line_len = i+1;
        break;
      } else if(first_chars[i]==0x09){
        kmer_len = i;
      } else if(i==200){
        printf("Control characters not found in kmers file!\n");
        upc_global_exit(-1);
      }
    }

    fseek(fin, sizeof(char), SEEK_END);
    long file_length = ftell(fin);
    //fseek(fin, 0L, SEEK_SET); //TODO
    line_count = file_length/line_len;

    printf("Line length: %d\n",line_len);
    printf("Kmer count:  %d\n",line_count);
    printf("Kmer length: %d\n",kmer_len);

    if(kmer_len!=KMER_LENGTH){
      printf("Detected kmer length (%d) doesn't match compiled length (%d)!",kmer_len,KMER_LENGTH);
      upc_global_exit(-5);
    }
    fclose(fin);
  }
  upc_barrier;

  int bucket_size = 3*line_count/THREADS;
  hashtable_size  = bucket_size*THREADS;
  kmers           = upc_all_alloc(THREADS, bucket_size*sizeof(kmer));

  //Declare all buckets unused
  upc_forall(int i=0;i<hashtable_size;i++;&kmers[i])
    kmers[i].l_ext = 0;

  if(MYTHREAD==0){
    FILE *fin = fopen(argv[1],"r");
    for(int i=0;i<line_count;i++){
      char kstr[KMER_LENGTH];
      char l_ext;
      char r_ext;
      fscanf(fin,"%" XSTR(KMER_LENGTH) "c %c%c ",kstr,&l_ext,&r_ext);
      AddKmer(kstr,l_ext,r_ext);
    }
    fclose(fin);
  }

  if(MYTHREAD==0){
    for(int i=0;i<hashtable_size;i++)
      if(kmers[i].l_ext!=0)
        PrintKmer(kmers[i]);
  }

	///////////////////////////////////////////
	// Your code for input file reading here //
	///////////////////////////////////////////
	upc_barrier;
	//inputTime += gettime();



	/** Graph construction **/
	//constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	upc_barrier;
	//constrTime += gettime();

	/** Graph traversal **/
	//traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	upc_barrier;
	//traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
