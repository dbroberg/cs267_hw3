#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include "packingDNAseq.h"
#include "kmer_hash.h"

int main(int argc, char **argv) {
   /* define current contig and extension variables, input kmer files */
   char cur_contig[MAXIMUM_CONTIG_SIZE] 
   char unpackedKmer[KMER_LENGTH+1] 
   char left_ext, right_ext, *input_UFX_name;
   int64_t contigID = 0, totBases = 0, ptr = 0;
   int64_t posInContig, nKmers, cur_chars_read, total_chars_to_read;
   kmer_t *cur_kmer_ptr; /* (defined in Algorithm 4) */
   start_kmer_t *startKmersList = NULL, *curStartNode; /* (defined in Algorithm 4) */
   unsigned char *working_buffer;
   /* ================================================ */
   /* ============== GRAPH CONSTRUCTION ============== */
   /* ================================================ */
   /* Read in input file input_UFX_name */
   /* Initialize lookup table that will be used for the DNA packing routines */
   init_LookupTable();  /* (defined in Algorithm 2) */
   /* Extract the number of k-mers in the input file */
   nKmers = getNumKmersInUFX(input_UFX_name); /* (defined in Algorithm 4) */
   hash_table_t *hashtable; /* (defined in Algorithm 4) */
   memory_heap_t memory_heap; /* (defined in Algorithm 4) */
   /* Create a hash table */
   hashtable = create_hash_table(nKmers, &memory_heap);  /* (defined in Algorithm 3) */
   /* Read the kmers from the input file and store them in the working_buffer */
   total_chars_to_read = nKmers * LINE_SIZE; 
   working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char)); 
   cur_chars_read = fread(working_buffer, sizeof(unsigned char),
						total_chars_to_read , inputFile); 
   /* Process the working_buffer and store the k-mers in the hash table */
   while (ptr < cur_chars_read) {
      /* Add k-mer to hash table */
      add_kmer(hashtable, &memory_heap, &working_buffer[ptr],
				 left_ext, right_ext);   /* (defined in Algorithm 3) */
      /* Create a list with the "guard" kmers as left (backward) extension */
      if (left_ext == 'F')
         addKmerToStartList(&memory_heap, &startKmersList);/*(defined in Algorithm 3)*/
      ptr += LINE_SIZE; /*Iterate to next start k-mer in working buffer */
   }

   /* ============================================= */
   /* ============== GRAPH TRAVERSAL ============== */
   /* ============================================= */
   curStartNode = startKmersList; 
   while (curStartNode != NULL ) {
      /* unpack the current seed kmer*/
      cur_kmer_ptr = curStartNode->kmerPtr;   
      unpackSequence((unsigned char*) cur_kmer_ptr->kmer,  
				(unsigned char*) unpackedKmer, KMER_LENGTH);
      /* Initialize current contig with the seed content */
      memcpy(cur_contig ,unpackedKmer, KMER_LENGTH * sizeof(char));
      posInContig = KMER_LENGTH;
      right_ext = cur_kmer_ptr->r_ext;
      /* Keep adding bases while not finding a terminal node */
      while (right_ext != 'F') {
         cur_contig[posInContig] = right_ext;
         posInContig++;
         cur_kmer_ptr = lookup_kmer(hashtable, (const unsigned char *) 
         &cur_contig[posInContig-KMER_LENGTH]); /* (defined in Algorithm 3) */
         right_ext = cur_kmer_ptr->r_ext;
      }
      /* Print the contig since we have found the corresponding terminal node */
      cur_contig[posInContig] = '\0';
      fprintf(serialOutputFile,"%s\n", cur_contig);
      contigID++;
      totBases += strlen(cur_contig);
      /* Move to the next start node in the list */
      curStartNode = curStartNode->next;
   }
   return 0;
}
