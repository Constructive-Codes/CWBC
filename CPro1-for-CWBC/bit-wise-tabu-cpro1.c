// cwbc_tabu_pr.c
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <limits.h>

static int n,w,d,s;
static uint64_t maskAll;
static uint64_t *cbCur, *cbBest;
static int *dmat;           // flat s×s distance matrix
static int penaltyCur, bestPenalty;
static int *tabu_remove, *tabu_add;  // size s×n, flat

// Elite pool for path relinking
#define MAX_ELITE 5
static uint64_t *Ecb[MAX_ELITE];
static int Epen[MAX_ELITE], Esize=0;

// Hyper‐parameters (may be overridden on command line)
static int tabu_min = 5;
static int tabu_max = 15;
static int max_no_improve_pr = 1000;
static int max_no_improve_restart = 10000;

// utility: popcount
static inline int popc(uint64_t x) {
    return __builtin_popcountll(x);
}

// copy current codebook to best
static void copy_cb(uint64_t *dst, uint64_t *src) {
    memcpy(dst, src, s * sizeof(uint64_t));
}

// add current cbCur/penaltyCur into elite pool
static void elite_add_current() {
    if (Esize < MAX_ELITE) {
        copy_cb(Ecb[Esize], cbCur);
        Epen[Esize] = penaltyCur;
        Esize++;
    } else {
        // replace worst if current is better
        int wi=0;
        for(int i=1;i<MAX_ELITE;i++)
            if(Epen[i] > Epen[wi]) wi = i;
        if(penaltyCur < Epen[wi]) {
            copy_cb(Ecb[wi], cbCur);
            Epen[wi] = penaltyCur;
        }
    }
}

// compute all pair‐distances and initial penalty
static void compute_all() {
    penaltyCur = 0;
    for(int i=0;i<s;i++){
        dmat[i*s + i] = 0;
        for(int j=i+1;j<s;j++){
            int dij = popc(cbCur[i] ^ cbCur[j]);
            dmat[i*s + j] = dmat[j*s + i] = dij;
            if(dij < d) penaltyCur += (d - dij);
        }
    }
}

// initialize a random constant‐weight codebook of size s
static void init_code() {
    cbCur = malloc(s * sizeof(uint64_t));
    cbBest = malloc(s * sizeof(uint64_t));
    tabu_remove = malloc(s * n * sizeof(int));
    tabu_add    = malloc(s * n * sizeof(int));
    memset(tabu_remove, 0, s*n*sizeof(int));
    memset(tabu_add,    0, s*n*sizeof(int));
    dmat = malloc(s * s * sizeof(int));
    // allocate elites
    for(int i=0;i<MAX_ELITE;i++){
        Ecb[i] = malloc(s * sizeof(uint64_t));
    }
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
    copy_cb(cbBest, cbCur);
    elite_add_current();
}

// path‐relinking between two elites Ecb[ip], Ecb[iq]
static void path_relink(int ip, int iq) {
    // local copy
    uint64_t *cbl = malloc(s*sizeof(uint64_t));
    int *dml = malloc(s*s*sizeof(int));
    memcpy(cbl, Ecb[ip], s*sizeof(uint64_t));
    // compute local dmat and penalty
    int penl=0;
    for(int i=0;i<s;i++){
        for(int j=i+1;j<s;j++){
            int dij = popc(cbl[i] ^ cbl[j]);
            dml[i*s+j] = dml[j*s+i] = dij;
            if(dij < d) penl += (d - dij);
        }
        dml[i*s+i] = 0;
    }
    // measure = sum HammingDist(cbl[i],Ecb[iq][i])
    int measure = 0;
    for(int i=0;i<s;i++){
        measure += popc(cbl[i] ^ Ecb[iq][i]);
    }
    // track best along path
    int bestLpen = penl;
    uint64_t *bestLcb = malloc(s*sizeof(uint64_t));
    memcpy(bestLcb, cbl, s*sizeof(uint64_t));

    // iteratively reduce measure by intra‐word 2‐bit swaps
    while(measure > 0){
        int best_k=-1, best_p=-1, best_q=-1;
        int best_pen2 = INT_MAX;
        // for each codeword k that still differs from target at some bits
        for(int k=0;k<s;k++){
            uint64_t ck = cbl[k];
            uint64_t tgt = Ecb[iq][k];
            uint64_t Pmask = ck & ~tgt;    // bits=1 here but target=0
            uint64_t Qmask = ~ck & tgt & maskAll; // bits=0 here but target=1
            if(Pmask==0 || Qmask==0) continue;
            // extract positions
            int Ppos[64], Qpos[64], pn=0, qn=0;
            for(int b=0;b<n;b++){
                if((Pmask>>b)&1) Ppos[pn++] = b;
                if((Qmask>>b)&1) Qpos[qn++] = b;
            }
            // scan all p,q
            for(int pi=0;pi<pn;pi++){
                int p = Ppos[pi];
                for(int qi=0;qi<qn;qi++){
                    int q = Qpos[qi];
                    // simulate swap in cbl[k]
                    uint64_t ck2 = ck ^ (1ULL<<p) ^ (1ULL<<q);
                    // compute delta penalty
                    int delta=0;
                    for(int t=0;t<s;t++){
                        if(t==k) continue;
                        int oldd = dml[k*s + t];
                        int newd = popc(ck2 ^ cbl[t]);
                        if(oldd < d) delta -= (d-oldd);
                        if(newd < d) delta += (d-newd);
                    }
                    int pen2 = penl + delta;
                    if(pen2 < best_pen2){
                        best_pen2 = pen2;
                        best_k=k; best_p=p; best_q=q;
                    }
                }
            }
        }
        if(best_k<0) break; // no path‐move found
        // apply the best move
        uint64_t ck = cbl[best_k];
        uint64_t ck2 = ck ^ (1ULL<<best_p) ^ (1ULL<<best_q);
        penl = best_pen2;
        // update dml for that row/col
        for(int t=0;t<s;t++){
            if(t==best_k) continue;
            int nd = popc(ck2 ^ cbl[t]);
            dml[best_k*s + t] = dml[t*s + best_k] = nd;
        }
        cbl[best_k] = ck2;
        measure -= 2;
        if(penl < bestLpen){
            bestLpen = penl;
            memcpy(bestLcb, cbl, s*sizeof(uint64_t));
        }
    }
    // if we improved the global best, accept bestLcb
    if(bestLpen < bestPenalty){
        // copy into global
        memcpy(cbCur, bestLcb, s*sizeof(uint64_t));
        compute_all();
        penaltyCur = penl = bestLpen;
        bestPenalty = penl;
        copy_cb(cbBest, cbCur);
        elite_add_current();
    }
    free(bestLcb);
    free(cbl);
    free(dml);
}

// main tabu‐search loop
static void do_search() {
    long iter=0;
    int noimp_pr=0, noimp_rt=0;
    while(1) {
        if(penaltyCur == 0) break;
        iter++;
        // pick a random violating pair (a,b)
        int viol_count=0;
        for(int i=0;i<s;i++){
            for(int j=i+1;j<s;j++){
                if(dmat[i*s+j] < d) viol_count++;
            }
        }
        if(viol_count == 0) {
            // should be caught by penaltyCur==0
            break;
        }
        int pick = rand() % viol_count;
        int a=-1,b=-1;
        for(int i=0;i<s;i++){
            for(int j=i+1;j<s;j++){
                if(dmat[i*s+j] < d){
                    if(pick-- == 0){
                        a=i; b=j;
                        goto FOUNDPAIR;
                    }
                }
            }
        }
        FOUNDPAIR: ;
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
                            int oldd = dmat[c*s + t];
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
        // if no non‐tabu found, do a raw best ignoring tabu
        if(best_c<0){
            for(int which=0; which<2; which++){
                int c = (which==0? a : b);
                int o = (which==0? b : a);
                uint64_t ck = cbCur[c], ok = cbCur[o];
                uint64_t Pmask = ck & ok;
                uint64_t Qmask = (~ck & ~ok) & maskAll;
                if(!Pmask||!Qmask) continue;
                int Ppos[64], Qpos[64], pn=0, qn=0;
                for(int bit=0;bit<n;bit++){
                    if((Pmask>>bit)&1) Ppos[pn++] = bit;
                    if((Qmask>>bit)&1) Qpos[qn++] = bit;
                }
                for(int pi=0;pi<pn;pi++){
                    int p = Ppos[pi];
                    for(int qi=0;qi<qn;qi++){
                        int q = Qpos[qi];
                        uint64_t ck2 = ck ^ (1ULL<<p) ^ (1ULL<<q);
                        int delta=0;
                        for(int t=0;t<s;t++){
                            if(t==c) continue;
                            int oldd = dmat[c*s + t];
                            int newd = popc(ck2 ^ cbCur[t]);
                            if(oldd < d) delta -= (d-oldd);
                            if(newd < d) delta += (d-newd);
                        }
                        int pen2 = penaltyCur + delta;
                        if(pen2 < best_new_pen){
                            best_new_pen = pen2;
                            best_c=c; best_p=p; best_q=q;
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
                dmat[c*s + t] = dmat[t*s + c] = nd;
            }
            cbCur[c] = ck2;
            // update penaltyCur
            penaltyCur = best_new_pen;
            // set new tabu tenures
            int ten = tabu_min + rand() % (tabu_max - tabu_min + 1);
            tabu_remove[c*n + best_p] = iter + ten;
            tabu_add   [c*n + best_q] = iter + ten;
        }
        // update best/global
        if(penaltyCur < bestPenalty){
            bestPenalty = penaltyCur;
            copy_cb(cbBest, cbCur);
            elite_add_current();
            noimp_pr = noimp_rt = 0;
        } else {
            noimp_pr++;
            noimp_rt++;
        }
        // intensify?
        if(noimp_pr >= max_no_improve_pr && Esize>=2){
            int i1 = rand()%Esize, i2;
            do{ i2 = rand()%Esize; } while(i2==i1);
            path_relink(i1, i2);
            noimp_pr = 0;
        }
        // diversify?
        if(noimp_rt >= max_no_improve_restart){
            // random shake: one random swap per word
            for(int rep=0; rep<s; rep++){
                int c = rand()%s;
                uint64_t ck = cbCur[c];
                // pick any two positions to swap
                int ones[64], zeros[64], on=0, zn=0;
                for(int b=0;b<n;b++){
                    if((ck>>b)&1) ones[on++] = b;
                    else                zeros[zn++] = b;
                }
                if(on>0 && zn>0){
                    int p = ones[rand()%on];
                    int q = zeros[rand()%zn];
                    uint64_t ck2 = ck ^ (1ULL<<p) ^ (1ULL<<q);
                    // update row c
                    for(int t=0;t<s;t++){
                        if(t==c) continue;
                        int nd = popc(ck2 ^ cbCur[t]);
                        dmat[c*s + t] = dmat[t*s + c] = nd;
                    }
                    cbCur[c] = ck2;
                }
            }
            compute_all();
            penaltyCur = 0;
            for(int i=0;i<s;i++){
                for(int j=i+1;j<s;j++){
                    int dij = dmat[i*s+j];
                    if(dij<d) penaltyCur += (d-dij);
                }
            }
            noimp_rt = noimp_pr = 0;
        }
    }
}

int main(int argc, char **argv) {
    if(argc < 6){
        fprintf(stderr,"Usage: %s n w d s seed [tabu_min tabu_max max_no_improve_pr max_no_improve_restart]\n", argv[0]);
        return 1;
    }
    n = atoi(argv[1]);
    w = atoi(argv[2]);
    d = atoi(argv[3]);
    s = atoi(argv[4]);
    unsigned seed = (unsigned)strtoul(argv[5],0,0);
    srand(seed);
    // optional hyperparams
    if(argc >= 10){
        tabu_min = atoi(argv[6]);
        tabu_max = atoi(argv[7]);
        max_no_improve_pr = atoi(argv[8]);
        max_no_improve_restart = atoi(argv[9]);
    }
    maskAll = (n==64 ? ~0ULL : ((1ULL<<n)-1ULL));
    init_code();
    do_search();
    // output final code (cbCur or cbBest)
    for(int i=0;i<s;i++){
        for(int b=0;b<n;b++){
            int bit = (cbCur[i] >> b) & 1U;
            printf("%d%c", bit, (b+1<n?' ':'\n'));
        }
    }
    return 0;
}