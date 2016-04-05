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

shared int lines_read[THREADS];

int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, traversalTime=0.0;

  current_smer = 0;

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();

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
  smers_size      = line_count/100/THREADS;
  kmers           = upc_all_alloc(THREADS, bucket_size*sizeof(kmer));
  smers           = (kmer_ptr*) calloc(line_count/100/THREADS, sizeof(kmer_ptr));

  //Declare all buckets unused
  if(MYTHREAD==0){
    printf("Initializing hash table\n");
  }
  upc_forall(int i=0;i<hashtable_size;i++;&kmers[i])
    kmers[i].l_ext = 0;

  upc_barrier;

  if(MYTHREAD==0){
    printf("Reading file\n");
  }

  //Open file
  FILE *fin = fopen(argv[1],"rb");
  //Seek location to start reading from
  int start_line = line_count/THREADS*MYTHREAD;
  int stop_line  = line_count/THREADS*(MYTHREAD+1);
  fseek(fin, start_line*line_len, SEEK_SET);
  //Determine when to stop reading
  if(MYTHREAD==THREADS-1)
    stop_line = line_count;

  lines_read[MYTHREAD] = 0;

  for(int i=start_line;i<stop_line;i++){
    lines_read[MYTHREAD]++;

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

	inputTime += gettime();

  double kmlist_offload=0;

  // for(int i=0;i<THREADS;i++){           
  //   upc_barrier;
  //   ListCount();
  // }

  kmlist_offload -= gettime();
  printf("Loading kmlist\n");
  for(int j=0;j<2;j++)                  //Do this twice to account for hash targets moving across boundaries
  for(int i=0;i<THREADS;i++){           //Move linked list into table
    upc_barrier;
    LoadKmList( (MYTHREAD+i)%THREADS );
  }
  upc_barrier;

  kmlist_offload += gettime();

  // for(int i=0;i<THREADS;i++){           
  //   upc_barrier;
  //   ListCount();
  // }

  // for(int i=0;i<THREADS;i++){           
  //   upc_barrier;
  //   ListPrint();
  // }
  

  /** Graph traversal **/
  // if(MYTHREAD==0){
  //   int total_lines_read=0;
  //   for(int i=0;i<THREADS;i++)
  //     total_lines_read+=lines_read[i];
  //   int total_nodes_inspected=0;
  //   for(int i=0;i<THREADS;i++)
  //     total_nodes_inspected+=nodes_inspected[i];
  //   int total_kmers_inserted=0;
  //   for(int i=0;i<THREADS;i++)
  //     total_kmers_inserted+=kmers_inserted[i];
  //   int total_kmers_added=0;
  //   for(int i=0;i<THREADS;i++)
  //     total_kmers_added+=kmers_added[i];
  //   int total_contigs=0;
  //   for(int i=0;i<THREADS;i++)
  //     total_contigs+=contig_count[i];

  //   printf("Lines read:      %d\n", total_lines_read);
  //   printf("Kmers inserted:  %d\n", total_kmers_inserted);
  //   printf("Kmers inspected: %d\n", total_nodes_inspected);
  //   printf("Kmers added:     %d\n", total_kmers_added);
  //   printf("Contigs found:   %d\n", total_contigs);
  //   printf("Current smer:    %d\n",  current_smer);
  //   printf("Generating contigs\n");
  // }

  // for(int i=0;i<THREADS;i++){
  //   upc_barrier;
  //   if(MYTHREAD==i){
  //     printf("THREAD %d\n",MYTHREAD);
  //     printf("Lines read: %d\n",lines_read[MYTHREAD]);
  //     printf("Kmers inserted: %d\n", kmers_inserted[MYTHREAD]);
  //     printf("Kmers inspected: %d\n", nodes_inspected[MYTHREAD]);
  //     printf("Kmers added: %d\n", kmers_added[MYTHREAD]);
  //   }
  // }

  assert(kmlist==NULL);

  traversalTime -= gettime();
  char output_name[100];
  sprintf(output_name,"/z/pgen-%d.out",MYTHREAD);
  
  FILE *fout = fopen(output_name, "wb");
  if(!fout){
    printf("Could not open output file!");
    upc_global_exit(-9);
  }

  if(MYTHREAD==0){
    printf("Generating contigs...\n");
  }
  for(int i=0;i<current_smer;i++)
    GenContig(fout,smers[i]);
	upc_barrier;
	traversalTime += gettime();

  if(MYTHREAD==0){
    int total_contigs_generated=0;
    for(int i=0;i<THREADS;i++)
      total_contigs_generated += contigs_generated[i];
    printf("Contigs generated: %d\n",total_contigs_generated);
  }

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
    printf("Kmlist offload: %f seconds\n", kmlist_offload);
    printf("Total time: %fs\n", (inputTime+traversalTime+kmlist_offload));
	}
	return 0;
}
