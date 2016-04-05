#ifndef _richard_kmer_h_
#define _richard_kmer_h_

#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <sys/time.h>

#ifndef TESTCOMPILE
  #include <upc.h>
#endif

#define XSTR(s) STR(s)
#define STR(s) #s

#define KMER_PACKED_LENGTH (KMER_LENGTH/4+1)
#define KMER_LAST_OFFSET   (2*(4-KMER_LENGTH%4))
#define KMER_LAST_MASK     (0xFF << KMER_LAST_OFFSET) //Puts zeros where there should be no data
#define CONTIG_SEQ_MAX     1000000

typedef unsigned char ksym_t;

#ifndef TESTCOMPILE
  #define EXIT(X) upc_global_exit(X)
#else 
  #include <stdlib.h>
  #define THREADS 1
  #define EXIT(X) exit(X)
  #define upc_memget memcpy
#endif

#ifndef TESTCOMPILE
  #pragma message "Compiling with UPC"
  shared int kmer_len;
  shared int line_len;
  shared int line_count;
#else
  #pragma message "Compiling without UPC"
  int kmer_len;
  int line_len;
  int line_count;
#endif

int hashtable_size;
int smers_size;
int locks_size;

typedef struct {
  ksym_t kmer[KMER_PACKED_LENGTH];
  ksym_t l_ext;
  ksym_t r_ext;
} kmer;

struct kmlist_t {
  int64_t hash;
  kmer km;
  struct kmlist_t *next;
};

typedef struct kmlist_t KmListNode;

KmListNode *kmlist = NULL;


#ifndef TESTCOMPILE
  //These are private pointers to shared space
  typedef shared kmer* kmer_ptr;
  typedef upc_lock_t*  lockptr;
  shared kmer_ptr*     smers;        //Start kmers!
  shared lockptr*      locks;
  shared int           current_smer; //Next empty start kmer location
  shared lockptr       current_smer_lock;
#else
  typedef kmer*        kmer_ptr;
  kmer_ptr*            smers;
  int current_smer;
#endif

kmer_ptr kmers;




static double gettime(void) {
    struct timeval tv;
    if (gettimeofday(&tv, NULL)) {
  perror("gettimeofday");
  abort();
    }
   return ((double)tv.tv_sec) + tv.tv_usec/1000000.0;
}



ksym_t BibitToSymbol(ksym_t symbol){
  switch(symbol){
    case 0: return 'A'; 
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
  }
  printf("Unknown sequence: %d",(int)symbol);
  return 255;
}

ksym_t SymbolToBibit(ksym_t symbol){
  switch(symbol){
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
  }
  printf("Unknown character: %c\n",symbol);
  return 255;
}

void PrintCharBits(ksym_t a){
  for(int i=7;i>=0;i--)
    printf("%d",(int)((a&(1<<i))!=0));
  printf(" ");
}

void PrintPackedKstr(const ksym_t *kmer){
  for(int i=0;i<KMER_PACKED_LENGTH;i++)
    PrintCharBits(kmer[i]);
  printf("\n");
}

void PrintKmer(const kmer km){ //TODO
  PrintPackedKstr(km.kmer);
  printf("\t%c %c\n",km.l_ext,km.r_ext);
}

void PackSequence(const ksym_t *unpacked, ksym_t *packed){
  char pack_offset = 6;

  for(int i=0;i<KMER_LENGTH;i++,unpacked++){
    if(pack_offset==6) //Make sure memory is initialized to 0
      *packed = 0;
    ksym_t bibit = SymbolToBibit(*unpacked);
    *packed     |= (bibit<<pack_offset);
    pack_offset -= 2;
    if(pack_offset==-2){
      pack_offset = 6;
      packed++;
    }
  }
}

void ConvertPackedToString(const ksym_t *packed, ksym_t *unpacked){
  for(int i=0;i<KMER_PACKED_LENGTH;i++){
    for(int n=6;n>=0;n-=2){
      if(i==KMER_PACKED_LENGTH-1 && n<KMER_LAST_OFFSET)
        return;
      char temp = packed[i];
      temp    >>= n;
      temp     &= 0x3;
      *unpacked = BibitToSymbol(temp);
      unpacked++;
    }
  }
}

void PrintPackedAsString(const ksym_t *packed){
  ksym_t unpacked[KMER_LENGTH];
  ConvertPackedToString(packed,unpacked);
  printf("%." XSTR(KMER_LENGTH) "s", unpacked);
}

int CompareKmer(const ksym_t *seq1, const ksym_t *seq2){
  return memcmp(seq1, seq2, KMER_PACKED_LENGTH)==0; //Are sequences equal?
}

int64_t HashKmer(int64_t hashtable_size, const ksym_t *kpacked){
  int64_t hashval;
  hashval = 5381;
  for(int i = 0; i < (KMER_PACKED_LENGTH); i++)
    hashval = kpacked[i] + (hashval << 5) + hashval;

  return hashval % hashtable_size;
}

int hash_collisions = 0;

void InsertKmer(int64_t hash, const kmer temp){
  kmers[hash] = temp;
  upc_lock(current_smer_lock);
  if(temp.l_ext=='F'){ // || r_ext=='F'){ //Only start from left-terminating kmers
    assert(current_smer!=smers_size);
    smers[current_smer] = &kmers[hash];
    current_smer++;
  }
  upc_unlock(current_smer_lock);
}

int64_t GetOpenBin(int64_t hash){
  //printf("%ld\n", hash);
  for(;kmers[hash].l_ext!=0;hash++,hash_collisions++){}
  return hash;
}

void AddKmer(const ksym_t *raw_kmer, ksym_t l_ext, ksym_t r_ext){
  kmer kmtemp;
  //printf("Packing: %.19s\n",raw_kmer);  
  PackSequence(raw_kmer,kmtemp.kmer);
  //PrintPackedKmer(kmer_packed);
  int64_t hash = HashKmer(hashtable_size,kmtemp.kmer);
  kmtemp.l_ext = l_ext;
  kmtemp.r_ext = r_ext;

  //Linear probing
  hash = GetOpenBin(hash);

  if(upc_threadof(&kmers[hash])==MYTHREAD){
    InsertKmer(hash,kmtemp);
  } else {
    KmListNode *tempnode = malloc(sizeof(KmListNode));
    tempnode->km         = kmtemp;
    tempnode->hash       = hash;
    tempnode->next       = kmlist;
    kmlist               = tempnode;
  }
}

void ShiftKmerLeft(ksym_t *kpacked){
  ksym_t high_bibit_right = 0;
  ksym_t high_bibit_this  = 0;
  for(int i=KMER_PACKED_LENGTH-1;i >= 0;i--){
    high_bibit_this   = kpacked[i] & 0xC0; //11 00 00 00
    high_bibit_this >>= 6;                 //Shift high bits right
    kpacked[i]      <<= 2;                 //Shift left two
    kpacked[i]       |= high_bibit_right;  //Take high bits from byte to the right
    high_bibit_right  = high_bibit_this;   //Save the high bits
  }
}

void ShiftKmerRight(ksym_t *kpacked){
  ksym_t low_bibit_left = 0;
  ksym_t low_bibit_this = 0;
  for(int i=0;i<KMER_PACKED_LENGTH;i++){
    low_bibit_this   = kpacked[i] & 0x3;   //00 00 00 11
    low_bibit_this <<=6;                   //Shift low bits left
    kpacked[i]     >>= 2;                  //Shift right two
    kpacked[i]      |= low_bibit_left;     //Take low bits from byte to the left
    low_bibit_left   = low_bibit_this;     //Save the low bits
  }
  kpacked[KMER_PACKED_LENGTH-1] &= KMER_LAST_MASK;
}

void SetKstrLeft(ksym_t *kpacked, ksym_t l_ext){
  kpacked[0] |= SymbolToBibit(l_ext)<<6;
}

void SetKstrRight(ksym_t *kpacked, ksym_t r_ext){
  kpacked[KMER_PACKED_LENGTH-1] |= SymbolToBibit(r_ext)<<KMER_LAST_OFFSET;
}

void ShiftAndAdd(ksym_t *kpacked, char direction, ksym_t l_ext, ksym_t r_ext){
  if(direction=='r'){
    ShiftKmerLeft(kpacked);
    SetKstrRight(kpacked,r_ext);
  } else if(direction=='l'){
    ShiftKmerRight(kpacked);
    SetKstrLeft(kpacked,l_ext);
  }
}

kmer_ptr FindKmer(const ksym_t *kstr){
  int64_t offset   = HashKmer(hashtable_size, kstr);
  kmer_ptr hashloc = &kmers[offset];
  for(;hashloc->l_ext!=0;hashloc++){
    ksym_t remote_kstr[KMER_PACKED_LENGTH];
    upc_memget(remote_kstr, hashloc->kmer, KMER_PACKED_LENGTH);
    if(CompareKmer(kstr,remote_kstr))
      return hashloc;
  }
  return NULL;
}

kmer_ptr NextKmer(const kmer_ptr km, char direction){
  assert(direction=='l' || direction=='r');

  ksym_t knextstr[KMER_PACKED_LENGTH];
  upc_memget(knextstr,km->kmer,KMER_PACKED_LENGTH);

  ShiftAndAdd(knextstr, direction, km->l_ext, km->r_ext);

  //printf("\tLooking for: ");PrintPackedAsString(knextstr);printf("\n");

  return FindKmer(knextstr);
}

void LoadKmList(int i){
  int low  = hashtable_size/THREADS*i;
  int high = hashtable_size/THREADS*(i+1);

  KmListNode *lptr = kmlist;
  while(lptr!=NULL){
    if(low<= lptr->hash && lptr->hash<high){ //Do we have permission to write here?
      lptr->hash = GetOpenBin(lptr->hash);
      if( !(low<= lptr->hash && lptr->hash<high) ) //Are we still in the locked range?
        lptr = lptr->next;
      InsertKmer(lptr->hash, lptr->km);
      KmListNode *next = lptr->next;
      free(lptr);
      lptr = next;
    } else {
      lptr = lptr->next;
    }
  }
  kmlist = lptr;
}

void GenContig(FILE *fout, kmer_ptr km){
  char direction;

  ksym_t contigseq[CONTIG_SEQ_MAX];
  int contiseq_len = 0;

  if(km->l_ext=='F'){
    direction = 'r';
  } else if(km->r_ext=='F'){
    direction = 'l';
  } else {
    printf("Cannot generate a contig from a kmer that doesn't start a sequence.\n");
    EXIT(-6);
  }

  kmer lkm = *km; //A local copy

  ConvertPackedToString(lkm.kmer,contigseq);
  contiseq_len = KMER_LENGTH;
  while( km!=NULL && ((direction=='l' && km->l_ext!='F') || (direction=='r' && km->r_ext!='F')) ){
    contigseq[contiseq_len++] = km->r_ext;
    if(contiseq_len==CONTIG_SEQ_MAX){
      printf("Reached maximum contig length!\n");
      EXIT(-7);
    }
    //PrintPackedAsString(km->kmer);
    //printf(" - %c %c\n",km->l_ext,km->r_ext);
    km = NextKmer(km, direction);
  }

  contigseq[contiseq_len] = '\0';
  fprintf(fout, "%s\n", contigseq);

  //printf("%s\n",contigseq);

  //if(km!=NULL){
  //  PrintPackedAsString(km->kmer);
  //  printf(" - %c %c\n",km->l_ext,km->r_ext);
  //}
}

#endif