#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define MAXW 16

int n,w,d,s;
int olimit;
int threshold;
int num;
int64_t numcandidates;
int64_t *candidates;
int64_t numsavecandidates;
int64_t *savecandidates;
int16_t *inithisto;
int64_t weightvec[MAXW+1]; // for possible overlaps 0...w (w<=16).
int64_t saveweights[MAXW+1]; 
int64_t *code;
int64_t *bestcand;

void init_candidates(int slicebits,int slice) {
  numsavecandidates = 0;
  int n2 = n-slicebits;
  int w2 = w-(__builtin_popcountll((int64_t) slice));
  int64_t limit = 1ULL << n2;
  int64_t cur = (1ULL << w2) - 1ULL;
  while (cur<limit) {
    for(int j=0;j<=MAXW;j++) {
      inithisto[numsavecandidates*(MAXW+1)+j] = 0;
    }
    int64_t codeword = (cur|(slice<<n2));
    int acceptflag = 1;
    for(int check=0;check<num;check++) {
      int overlap = __builtin_popcountll(codeword&(code[check]));
      if(overlap>olimit) {
	acceptflag = 0;
	break;
      }
      inithisto[numsavecandidates*(MAXW+1)+overlap]++;
    }
    if(acceptflag) {
      savecandidates[(numsavecandidates<<1)] = (cur|(slice<<n2));
      savecandidates[(numsavecandidates<<1)+1] = 0;  // initial dot product
      numsavecandidates++;
    }
    // Gosper's hack f
    int64_t u = cur & (int64_t)(-(int64_t)cur);
    int64_t v = cur + u;
    cur = v+(((v^cur)/u)>> 2);
  }
}

void reinit_candidates() {
  for(int64_t i=0;i<numsavecandidates;i++) {
    candidates[i<<1] = savecandidates[i<<1];
    candidates[(i<<1)+1] = 0;
    for(int j=0;j<=MAXW;j++) {
      candidates[(i<<1)+1] += inithisto[i*(MAXW+1)+j]*weightvec[j];
    }
  }
  numcandidates = numsavecandidates;
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
  int initnum = num;
  if(numcandidates==0) return;
  while((numcandidates>0)&&(num<s)) {
    int64_t lastadd = code[num-1];
    int numbest = 0;
    int64_t best = 0;
    for(int64_t i=0;i<(numcandidates<<1);i+=2) {
      int overlap = __builtin_popcountll(lastadd&(candidates[i]));
      if((overlap<=olimit)||(num==initnum)) { // allow
	if(num>initnum) { // otherwise already taken care of
	  candidates[i+1] += weightvec[overlap];
	}
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
      num++;
    }
  }
}


int main(int argc, char **argv) {
    if(argc<9) {
      fprintf(stderr,"Usage: %s n w d s seed threshold slicebits slicereps\n",argv[0]);
      return 1;
    }
    n = atoi(argv[1]);
    w = atoi(argv[2]);
    d = atoi(argv[3]);
    s = atoi(argv[4]);
    int seed = atoi(argv[5]);
    threshold = atoi(argv[6]);
    int slicebits = atoi(argv[7]);
    int slicereps = atoi(argv[8]);
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
    savecandidates = malloc(n_choose_w*sizeof(int64_t)*2); 
    bestcand = malloc(n_choose_w*sizeof(int64_t)); // for use in inner loop
    code = malloc(s*sizeof(int64_t));
    int64_t *slicebestcode = malloc(s*sizeof(int64_t));
    inithisto = malloc(n_choose_w*sizeof(int16_t)*(MAXW+1));
    while(1) {
      int slicemax = 1<<slicebits;
      int sliceshuf[slicemax];
      for(int i=0;i<slicemax;i++) {
	sliceshuf[i] = i;
      }
      for(int i=0;i<slicemax-1;i++) {
	int swap = (rand()%(slicemax-i))+i;
	int tmp = sliceshuf[i];
	sliceshuf[i] = sliceshuf[swap];
	sliceshuf[swap] = tmp;
      }
      num = 0;
      for(int i=0;i<slicemax;i++) {
	init_candidates(slicebits,sliceshuf[i]);
	int slicestart = num;
	int slicebest = num;
	for(int j=0;j<slicereps;j++) {
	  randomize_weightvec();
	  for(int i=0;i<=olimit;i++) {
	    saveweights[i] = weightvec[i];
	  }
	  reinit_candidates();
	  num = slicestart;
	  generate_code();  // re-randomizes weightvec
	  if(num>slicebest) {
	    printf("got %d with ",num);
	    for(int i=0;i<=olimit;i++) printf("%d(%d):%lld ",i,(w-i)*2,(long long) weightvec[i]);
	    printf(" | ");
	    for(int i=0;i<=olimit;i++) printf("%d(%d):%lld ",i,(w-i)*2,(long long) saveweights[i]);
	    printf("\n");
	    if(num>threshold) {
	      printf("code:\n");
	      for(int i=0;i<num;i++) {
		for(int j=0;j<n;j++) {
		  printf("%d ",(int) ((code[i]>>j)&1LL));
		}
		printf("\n");
	      }
	    }
	    fflush(stdout);
	    slicebest = num;
	    for(int k=slicestart;k<num;k++) {
	      slicebestcode[k] = code[k];
	    }
	  }
	  if(i==0) break;
	}
	for(int k=slicestart;k<slicebest;k++) {
	  code[k] = slicebestcode[k];
	}
	num = slicebest;
	printf("after slice %d:%d - %ld candidates - %d\n",i,sliceshuf[i],numsavecandidates,num);
      }
      printf("\n");
    }
}

    
