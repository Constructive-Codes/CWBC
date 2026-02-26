#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define MAXW 16

int n,w,d,s;
int olimit;
int threshold;
int64_t numcandidates;
int64_t *candidates;
int64_t weightvec[MAXW+1]; // for possible overlaps 0...w (w<=16).
int64_t *code;
int64_t *bestcand;

void init_candidates() {
    int64_t i = 0;
    int64_t limit = 1ULL << n;
    int64_t cur = (1ULL << w) - 1ULL;
    while (cur<limit) {
      candidates[i++] = cur;
      candidates[i++] = 0;  // initial dot product
      // Gosper's hack f
      int64_t u = cur & (int64_t)(-(int64_t)cur);
      int64_t v = cur + u;
      cur = v+(((v^cur)/u)>> 2);
    }
}

void randomize_weightvec() {
  int64_t maxweight = 0;
  for(int i=0;i<=MAXW;i++) {
    weightvec[i] = rand();
    if(rand()&1) weightvec[i] = -weightvec[i];
    if(weightvec[i]>maxweight) maxweight = weightvec[i];
  }
  while(weightvec[olimit]<maxweight) weightvec[olimit] = rand();
}


void generate_code() {
  randomize_weightvec();
  code[0] = candidates[0];  // wlog
  int num = 1;
  int branches[s];
  while((numcandidates>0)&&(num<s)) {
    int64_t lastadd = code[num-1];
    int numbest = 0;
    int64_t best = 0;
    for(int64_t i=0;i<(numcandidates<<1);i+=2) {
      int overlap = __builtin_popcountll(lastadd&(candidates[i]));
      if(overlap<=olimit) { // allow
	candidates[i+1] += weightvec[overlap];
	if((candidates[i+1]>best)||(numbest==0)) {
	  best = candidates[i+1];
	  bestcand[0] = candidates[i];	  
	  numbest = 1;
	} else if(candidates[i+1]==best) {
	  bestcand[numbest] = candidates[i];
	  numbest++;
	}
      } else { // remove
	if(i!=(numcandidates<<1)-2) {
	  candidates[i] = candidates[(numcandidates<<1)-2];
	  candidates[i+1] = candidates[(numcandidates<<1)-1];
	}
	numcandidates--;
	i-=2;
      }
    }
    if(numbest>0) {
      code[num] = bestcand[rand()%numbest];
      branches[num] = numbest;
      num++;
    }
  }
  printf("got %d with ",num);
  for(int i=0;i<=olimit;i++) printf("%d(%d):%lld ",i,(w-i)*2,(long long) weightvec[i]); 
  printf("\n");
  if(num>threshold) {
    printf("branching: ");
    for(int i=1;i<num;i++) printf("%d:%d ",i,branches[i]);
    printf("\n");
    printf("code:\n");
    for(int i=0;i<num;i++) {
      for(int j=0;j<n;j++) {
	printf("%d ",(int) ((code[i]>>j)&1LL));
      }
      printf("\n");
    }
  }
  fflush(stdout);
}


int main(int argc, char **argv) {
    if(argc<7) {
      fprintf(stderr,"Usage: %s n w d s seed threshold\n",argv[0]);
      return 1;
    }
    n = atoi(argv[1]);
    w = atoi(argv[2]);
    d = atoi(argv[3]);
    s = atoi(argv[4]);
    int seed = atoi(argv[5]);
    threshold = atoi(argv[6]);
    srand(seed);

    olimit = w-(d/2); // max overlapping 1 bits between two compatible codewords
    
    int64_t n_choose_w = 1;
    for(int64_t i=n;i>(n-w);i--) {
      n_choose_w *= i;
    }
    for(int64_t i=2;i<=w;i++) {
      n_choose_w /= i;
    }
    if(n_choose_w>=(1LL<<31)) {
      printf("unsupported (n choose w) %lld\n",(long long int) n_choose_w);
      fflush(stdout);
      exit(1);
    }
    candidates = malloc(n_choose_w*sizeof(int64_t)*2); // store each codeword together with its current dot product
    bestcand = malloc(n_choose_w*sizeof(int64_t)); // for use in inner loop
    code = malloc(s*sizeof(int64_t));
    while(1) {
      init_candidates();
      numcandidates = n_choose_w;
      generate_code();
    }
}

    
