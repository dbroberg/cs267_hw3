//gcc test.c -DKMER_LENGTH=19 -DTESTCOMPILE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include "richard-kmer.h"

unsigned char NumToChar(int i){
  switch(i%4){
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
  }
  assert(0);
  return 0;
}

int CheckRawKstr(ksym_t *k_raw_str){
  for(int i=0;i<KMER_LENGTH;i++)
    if(SymbolToBibit(k_raw_str[i])==255)
      return 0;
  return 1;
}

int main(int argc, char *argv[]){
  current_smer = 0;

  printf("Kmer length:        %d\n",KMER_LENGTH);
  printf("Kmer packed length: %d\n",KMER_PACKED_LENGTH);
  printf("Kmer last offset:   %d\n",KMER_LAST_OFFSET);

  printf("Testing packing routine:\n");
  ksym_t packed[KMER_PACKED_LENGTH];
  {
    ksym_t unpacked[KMER_LENGTH];
    for(int i=0;i<KMER_LENGTH;i++)
      unpacked[i] = NumToChar(i);
    PackSequence(unpacked,packed);

    PrintPackedKstr(packed);
  }

  {
    printf("Testing left shift:\n");
    ksym_t testseq[KMER_PACKED_LENGTH];
    memcpy(testseq,packed,KMER_PACKED_LENGTH);
    for(int i=0;i<10;i++){
      PrintPackedKstr(testseq);
      ShiftKmerLeft(testseq);
    }
  }

  {
    printf("Testing right shift:\n");
    ksym_t testseq[KMER_PACKED_LENGTH];
    memcpy(testseq,packed,KMER_PACKED_LENGTH);
    for(int i=0;i<10;i++){
      PrintPackedKstr(testseq);
      ShiftKmerRight(testseq);
    }
  }

  {
    printf("Testing set left:\n");
    ksym_t testseq[KMER_PACKED_LENGTH];
    for(int i=0;i<10;i++){
      memset(testseq,0,KMER_PACKED_LENGTH);
      SetKstrLeft(testseq,NumToChar(i));
      PrintPackedKstr(testseq);
    }
  } 

  {
    printf("Testing set right:\n");
    ksym_t testseq[KMER_PACKED_LENGTH];
    for(int i=0;i<10;i++){
      memset(testseq,0,KMER_PACKED_LENGTH);
      SetKstrRight(testseq,NumToChar(i));
      PrintPackedKstr(testseq);
    }
  } 

  {
    printf("Testing shift left and add right:\n");
    ksym_t testseq[KMER_PACKED_LENGTH];
    memcpy(testseq,packed,KMER_PACKED_LENGTH);    
    for(int i=0;i<10;i++){
      PrintPackedKstr(testseq);
      ShiftAndAdd(testseq,'l','C','C');
    }
  } 

  {
    printf("Testing shift right and add left:\n");
    ksym_t testseq[KMER_PACKED_LENGTH];
    memcpy(testseq,packed,KMER_PACKED_LENGTH);      
    for(int i=0;i<10;i++){
      PrintPackedKstr(testseq);
      ShiftAndAdd(testseq,'r','C','C');
    }
  } 





  if(argc!=2){
    printf("Need to specify program file!\n");
    return -1;
  }

  FILE *fin = fopen(argv[1],"r");
  if(!fin){
    printf("Failed to open kmers file!\n");
    return -2;
  }

  char first_chars[200];
  if(fread(first_chars,sizeof(char),200,fin)!=200){
    printf("Failed to read first part of kmers file!\n");
    return -3;
  }

  for(int i=0;i<200;i++){
    if(first_chars[i]=='\n'){
      line_len = i+1;
      break;
    } else if(first_chars[i]==0x09){
      kmer_len = i;
    } else if(i==200){
      printf("Control characters not found in kmers file!\n");
      return -4;
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
    return -5;
  }
  fclose(fin);

  int bucket_size = 3*line_count/THREADS;
  hashtable_size  = bucket_size*THREADS;
  printf("Setting hash table size: %d\n",hashtable_size);
  kmers           = calloc(hashtable_size, sizeof(kmer));
  smers           = calloc(line_count/1000, sizeof(kmer_ptr));

  //Declare all buckets unused
  printf("Initializing hashtable\n");
  for(int i=0;i<hashtable_size;i++)
    kmers[i].l_ext = 0;

  printf("Reading file\n");
  fin = fopen(argv[1],"r");
  for(int i=0;i<line_count;i++){
    long line_start = ftell(fin);
    ksym_t kstr[KMER_LENGTH];
    ksym_t l_ext;
    ksym_t r_ext;
    fread(kstr,sizeof(ksym_t),KMER_LENGTH,fin);
    fgetc(fin);
    l_ext = fgetc(fin);
    r_ext = fgetc(fin);
    fgetc(fin);
    //fscanf(fin,"%" XSTR(KMER_LENGTH) "c %c%c ",kstr,&l_ext,&r_ext);
    if(!CheckRawKstr(kstr) || (1936870<=i && i<=1936880)){
      printf("Problem on line: %d\n",i);
      printf("At file position: %ld\n",line_start);
      printf("Read %." XSTR(KMER_LENGTH) "s %c %c\n",kstr,l_ext,r_ext);
      for(int i=0;i<KMER_LENGTH;i++)
        printf(" %02x", kstr[i]);
      printf("\n");
    }
    AddKmer(kstr,l_ext,r_ext);
  }
  fclose(fin);

  
  printf("Printing beginning of hashtable");
  for(int i=0;i<hashtable_size;i++)
    if(kmers[i].l_ext!=0)
      PrintKmer(kmers[i]);
    else
      printf("-empty slot-\n");


  //Distributed graph traversal!
  //upc_forall(int i=0;i<current_smer;i++;&smers[i]){
  //  if(MYTHREAD==0){
  //    PrintKmer(*smers[i]);
  //  }
 // }

  /** Print timing and output info **/
  /***** DO NOT CHANGE THIS PART ****/
  // if(MYTHREAD==0){
  //   printf("%s: Input set: %s\n", argv[0], argv[1]);
  //   printf("Number of UPC threads: %d\n", THREADS);
  //   printf("Input reading time: %f seconds\n", inputTime);
  //   printf("Graph construction time: %f seconds\n", constrTime);
  //   printf("Graph traversal time: %f seconds\n", traversalTime);
  // }
  // return 0;

  return 0;
}
