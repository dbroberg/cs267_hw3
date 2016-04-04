#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include "richard-kmer.h"
//#include <upc_io.h>

//#include "packingDNAseq.h"
//#include "kmer_hash.h"

int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

  current_smer = 0;

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
    unsigned char first_chars[200];
    if(fread(first_chars,sizeof(unsigned char),200,fin)!=200){
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

    fseek(fin, sizeof(unsigned char), SEEK_END);
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
  smers_size      = line_count/100;
  kmers           = upc_all_alloc(THREADS, bucket_size*sizeof(kmer));
  smers           = upc_all_alloc(THREADS, sizeof(shared kmer*)*line_count/1000);

  //Declare all buckets unused
  if(MYTHREAD==0){
    printf("Initializing hash table\n");
  }
  upc_forall(int i=0;i<hashtable_size;i++;&kmers[i])
    kmers[i].l_ext = 0;

  upc_barrier;

  if(MYTHREAD==0){
    printf("Reading file\n");
    FILE *fin = fopen(argv[1],"rb");
    for(int i=0;i<line_count;i++){
      ksym_t kstr[KMER_LENGTH];
      ksym_t l_ext;
      ksym_t r_ext;
      //fscanf(fin,"%" XSTR(KMER_LENGTH) "c %c%c ",kstr,&l_ext,&r_ext);

      if(fread(kstr,sizeof(ksym_t),KMER_LENGTH,fin)!=KMER_LENGTH){
        printf("Didn't read enough characters!\n");
      }
      fgetc(fin);
      l_ext = fgetc(fin);
      r_ext = fgetc(fin);
      fgetc(fin);

      if(ferror(fin))
        printf("Error reading file: %s\n", strerror(ferror(fin)));
      AddKmer(kstr,l_ext,r_ext);
    }
    fclose(fin);

    //Print kmers for debugging purposes
    // for(int i=0;i<hashtable_size;i++)
    //   if(kmers[i].l_ext!=0)
    //     PrintKmer(kmers[i]);
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

  //Distributed graph traversal!
  if(MYTHREAD==0){
    printf("Generating contigs\n");
  }
  upc_forall(int i=0;i<current_smer;i++;&smers[i]){
    GenContig(smers[i]);
  }

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
