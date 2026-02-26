// cwbc_tabu_pr.c
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <limits.h>

static int n,w,d,s;
static uint64_t maskAll;
static uint64_t *cbCur;
static int *dmat;           // flat s×s distance matrix
static int *conflictlist;
int numconflicts;
static int penaltyCur, bestPenalty;
static int *tabu_remove, *tabu_add;  // size s×n, flat

// Hyper‐parameters (may be overridden on command line)
static int tabu_min = 5;
static int tabu_max = 15;
static int max_no_improve_restart = 10000;

// utility: popcount
static inline int popc(uint64_t x) {
    return __builtin_popcountll(x);
}

// compute all pair‐distances and initial penalty
static void compute_all() {
    penaltyCur = 0;
    numconflicts = 0;
    for(int i=0;i<s;i++){
        for(int j=i+1;j<s;j++){
            int dij = popc(cbCur[i] ^ cbCur[j]);
	    int loc = (i*s+j)<<1;
            dmat[loc] = dij;
            if(dij < d) {
	      penaltyCur += (d - dij);
	      dmat[loc+1] = numconflicts;
	      conflictlist[numconflicts<<1] = i;
	      conflictlist[(numconflicts<<1)+1] = j;
	      numconflicts++;
	    }
        }
    }
}

// initialize a random constant‐weight codebook of size s
static void init_code() {
    memset(tabu_remove, 0, s*n*sizeof(int));
    memset(tabu_add,    0, s*n*sizeof(int));
    // generate s distinct vectors
    int *pos = malloc(n * sizeof(int));
    for(int i=0;i<n;i++) pos[i] = i;
    for(int i=0;i<s;i++){
        uint64_t cw;
        int dup;
        do{
            cw = 0;
            // partial shuffle for w positions
            for(int j=0;j<w;j++){
                int r = j + rand() % (n - j);
                int tmp = pos[j]; pos[j] = pos[r]; pos[r] = tmp;
                cw |= (1ULL << pos[j]);
            }
            dup = 0;
            for(int k=0;k<i;k++){
                if(cbCur[k] == cw){ dup=1; break; }
            }
        } while(dup);
        cbCur[i] = cw;
    }
    free(pos);
    compute_all();
    bestPenalty = penaltyCur;
}


// main tabu‐search loop
static void do_search() {
    long iter=0;
    int noimp_pr=0, noimp_rt=0;
    while(1) {
        if(penaltyCur == 0) break;
        iter++;

        // pick a random violating pair (a,b)
	int pick = rand()%numconflicts;
        int a = conflictlist[pick<<1];
	int b = conflictlist[(pick<<1)+1];
        // search best 2‐bit swap in either codeword a or b
        int best_c=-1, best_p=-1, best_q=-1;
        int best_new_pen = INT_MAX;
        // try each of the two words
        for(int which=0; which<2; which++){
            int c = (which==0? a : b);
            int o = (which==0? b : a);
            uint64_t ck = cbCur[c], ok = cbCur[o];
            // to increase d(c,o), swap within c on positions both=1 -> both=0
            uint64_t Pmask = ck & ok;               // both 1
            uint64_t Qmask = (~ck & ~ok) & maskAll; // both 0
            if(Pmask && Qmask){
                int Ppos[64], Qpos[64], pn=0, qn=0;
                for(int bit=0;bit<n;bit++){
                    if((Pmask>>bit)&1) Ppos[pn++] = bit;
                    if((Qmask>>bit)&1) Qpos[qn++] = bit;
                }
                for(int pi=0;pi<pn;pi++){
                    int p = Ppos[pi];
                    for(int qi=0;qi<qn;qi++){
                        int q = Qpos[qi];
                        // check tabu on c,p→remove and c,q→add
                        int tr = tabu_remove[c*n + p];
                        int ta = tabu_add   [c*n + q];
                        // simulate delta penalty
                        uint64_t ck2 = ck ^ (1ULL<<p) ^ (1ULL<<q);
                        int delta=0;
                        for(int t=0;t<s;t++){
                            if(t==c) continue;
			    int oldd;
			    if(c<t) {
			      oldd = dmat[(c*s + t)<<1];
			    } else {
			      oldd = dmat[(t*s + c)<<1];
			    }
                            int newd = popc(ck2 ^ cbCur[t]);
                            if(oldd < d) delta -= (d-oldd);
                            if(newd < d) delta += (d-newd);
                        }
                        int pen2 = penaltyCur + delta;
                        // if tabu but not aspiration => skip
                        if((tr>iter || ta>iter) && pen2 >= bestPenalty)
                            continue;
                        if(pen2 < best_new_pen){
                            best_new_pen = pen2;
                            best_c = c;
                            best_p = p;
                            best_q = q;
                        }
                    }
                }
            }
        }
        // apply the best move
        if(best_c>=0){
            int c = best_c;
            uint64_t ck = cbCur[c];
            uint64_t ck2 = ck ^ (1ULL<<best_p) ^ (1ULL<<best_q);
            // update distances row c
            for(int t=0;t<s;t++){
                if(t==c) continue;
                int nd = popc(ck2 ^ cbCur[t]);
		int pos;
		if(c<t) {
		  pos = (c*s + t)<<1;
		  if((dmat[pos]>=d)&&(nd<d)) { // add to end of conflictlist
		    conflictlist[numconflicts<<1] = c;
		    conflictlist[(numconflicts<<1)+1] = t;
		    dmat[pos] = nd;
		    dmat[pos+1] = numconflicts;
		    numconflicts++;		    
		  } else if((dmat[pos]<d)&&(nd>=d)) { // remove from conflictlist; move end of conflictlist into its place
		    int oldpos = dmat[pos+1];
		    if(oldpos!=numconflicts-1) {
		      conflictlist[oldpos<<1] = conflictlist[(numconflicts-1)<<1];
		      conflictlist[(oldpos<<1)+1] = conflictlist[((numconflicts-1)<<1)+1];
		      dmat[(((conflictlist[oldpos<<1])*s+(conflictlist[(oldpos<<1)+1]))<<1)+1] = oldpos;
		    }
		    numconflicts--;
		    dmat[pos] = nd;
		  } else {
		    dmat[pos] = nd;
		  }
		} else {
		  pos = (t*s + c)<<1;
		  if((dmat[pos]>=d)&&(nd<d)) { // add to end of conflictlist
		    conflictlist[numconflicts<<1] = t;
		    conflictlist[(numconflicts<<1)+1] = c;
		    dmat[pos] = nd;
		    dmat[pos+1] = numconflicts;
		    numconflicts++;		    
		  } else if((dmat[pos]<d)&&(nd>=d)) {
		    int oldpos = dmat[pos+1];
		    if(oldpos!=numconflicts-1) {
		      conflictlist[oldpos<<1] = conflictlist[(numconflicts-1)<<1];
		      conflictlist[(oldpos<<1)+1] = conflictlist[((numconflicts-1)<<1)+1];
		      dmat[(((conflictlist[oldpos<<1])*s+(conflictlist[(oldpos<<1)+1]))<<1)+1] = oldpos;
		    }		    
		    numconflicts--;
		    dmat[pos] = nd;
		  } else {
		    dmat[pos] = nd;
		  }
		}
            }
            cbCur[c] = ck2;
            // update penaltyCur
            penaltyCur = best_new_pen;
            // set new tabu tenures
            int ten = tabu_min + rand() % (tabu_max - tabu_min + 1);
            tabu_remove[c*n + best_p] = iter + ten;
            tabu_add   [c*n + best_q] = iter + ten;
        }
	//	printf("%ld %d\n",iter,numconflicts);
	//	fflush(stdout);
        // update best/global
        if(penaltyCur < bestPenalty){
            bestPenalty = penaltyCur;
            noimp_pr = noimp_rt = 0;
        } else {
            noimp_pr++;
            noimp_rt++;
        }
        // restart?
        if((best_c<0)||(noimp_rt >= max_no_improve_restart)){
	  int seen = 0;
	  int penaltycheck = 0;
	  for (int i=0;i<s;i++) {
	    for (int j=i+1;j<s;j++) {
	      int dij = popc(cbCur[i] ^ cbCur[j]);
	      if(dij!=dmat[((i*s+j)<<1)]) {
		fprintf(stderr,"Distance failure: %d %d %d %d\n",i,j,dij,dmat[((i*s+j)<<1)]);
	      }
	      if (dmat[((i*s+j)<<1)] < d) {
		seen++;
		penaltycheck += (d - dij);
		int pos = dmat[((i*s+j)<<1)+1];
		if((conflictlist[(pos<<1)]!=i)||(conflictlist[(pos<<1)+1]!=j)) {
		  fprintf(stderr, "Mismatch: %d %d %d %d %d\n",i,j,pos,conflictlist[(pos<<1)],conflictlist[(pos<<1)+1]);
		}
	      }
	    }
	  }
	  if (seen != numconflicts) {
	    fprintf(stderr, "Inconsistency: seen=%d numconflicts=%d at iter=%ld\n",
		    seen, numconflicts, iter);
	  }
	  if(penaltycheck!=penaltyCur) {
	    fprintf(stderr,"Penalty mismatch: %ld: %d %d (%d)\n",iter,penaltycheck,penaltyCur,numconflicts);
	  }
	  init_code();
	  noimp_rt = noimp_pr = 0;
	  iter = 0;
	  //	  printf("restarted at %ld\n",iter);
	  //	  fflush(stdout);
        }
    }
}

int main(int argc, char **argv) {
    if(argc < 6){
        fprintf(stderr,"Usage: %s n w d s seed [tabu_min tabu_max max_no_improve_restart]\n", argv[0]);
        return 1;
    }
    n = atoi(argv[1]);
    w = atoi(argv[2]);
    d = atoi(argv[3]);
    s = atoi(argv[4]);
    unsigned seed = (unsigned)strtoul(argv[5],0,0);
    srand(seed);
    // optional hyperparams
    if(argc >= 9){
        tabu_min = atoi(argv[6]);
        tabu_max = atoi(argv[7]);
        max_no_improve_restart = atoi(argv[8]);
    }
    maskAll = (n==64 ? ~0ULL : ((1ULL<<n)-1ULL));
      
    cbCur = malloc(s * sizeof(uint64_t));
    tabu_remove = malloc(s * n * sizeof(int));
    tabu_add    = malloc(s * n * sizeof(int));
    dmat = malloc(s * s * sizeof(int) * 2);
    conflictlist = malloc(s * s * sizeof(int));

    init_code();
    do_search();
    // output final code (cbCur)
    for(int i=0;i<s;i++){
        for(int b=0;b<n;b++){
            int bit = (cbCur[i] >> b) & 1U;
            printf("%d%c", bit, (b+1<n?' ':'\n'));
        }
    }
    return 0;
}
