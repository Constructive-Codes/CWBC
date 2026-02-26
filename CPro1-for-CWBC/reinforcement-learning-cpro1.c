#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// count 1-bits in a 64-bit word
static inline int popcount64(uint64_t x) {
    return __builtin_popcountll(x);
}

// draw a uniform [0,1) double
static inline double rand_uniform() {
    return rand() / (RAND_MAX + 1.0);
}

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr,
            "Usage: %s n w d s seed [alpha] [beta]\n"
            "  n,w,d,s,seed  integer parameters of CWBC\n"
            "  alpha (opt.)  policy-gradient learning rate [default=0.01]\n"
            "  beta  (opt.)  softmax temperature          [default=1.0]\n",
            argv[0]);
        return 1;
    }

    // parse required parameters
    int n   = atoi(argv[1]);
    int w   = atoi(argv[2]);
    int d   = atoi(argv[3]);
    int s   = atoi(argv[4]);
    int seed= atoi(argv[5]);

    // optional hyperparameters
    double alpha = 0.01;
    double beta  = 1.0;
    if (argc >= 7) alpha = atof(argv[6]);
    if (argc >= 8) beta  = atof(argv[7]);

    if (n <= 0 || n >= 64 || w <= 0 || w > n || d < 0 || d > n || s <= 0) {
        fprintf(stderr, "Invalid parameters.\n");
        return 1;
    }

    srand(seed);

    // 1) Generate all weight-w bitmasks of length n via Gosper's hack
    size_t cap = 1024;
    uint64_t *P0 = malloc(cap * sizeof(uint64_t));
    if (!P0) { fprintf(stderr, "Out of memory\n"); return 1; }
    size_t M = 0;
    uint64_t limit = (n == 64 ? 0ULL : (1ULL << n));
    uint64_t comb = ((w>=64)?0ULL:((1ULL << w) - 1ULL));
    while (comb && (n==64 || comb < limit)) {
        if (M >= cap) {
            cap <<= 1;
            P0 = realloc(P0, cap * sizeof(uint64_t));
            if (!P0) { fprintf(stderr, "Out of memory\n"); return 1; }
        }
        P0[M++] = comb;
        // Gosper's hack for next combination
        uint64_t u = comb & (uint64_t)(-(int64_t)comb);
        uint64_t v = comb + u;
        if (v == 0) break;
        comb = v + (((v ^ comb) / u) >> 2);
    }
    if (M == 0) {
        fprintf(stderr, "No weight-%d vectors in length %d\n", w, n);
        return 1;
    }

    // B = number of distance bins = n+1 (distances 0..n)
    int B = n + 1;

    // allocate core buffers
    uint64_t *code = malloc(s * sizeof(uint64_t));
    int *valid_ids = malloc(M * sizeof(int));
    int *features_mem = malloc(s * B * sizeof(int));
    double *Efeatures_mem = malloc(s * B * sizeof(double));
    int *f_all = malloc(M * B * sizeof(int));
    double *Theta = calloc(B, sizeof(double));
    double *E_f = malloc(B * sizeof(double));
    double *score = malloc(M * sizeof(double));
    double *exp_score = malloc(M * sizeof(double));

    if (!code || !valid_ids || !features_mem || !Efeatures_mem
        || !f_all || !Theta || !E_f || !score || !exp_score) {
        fprintf(stderr, "Out of memory\n");
        return 1;
    }

    // Reinforcement-learning loop: repeat episodes until success
    for (;;) {
        int code_count = 0;
        int step_count = 0;

        // Build one episode
        for (int t = 0; t < s; t++) {
            // 2) Filter valid candidates w.r.t. current code
            int valid_count = 0;
            for (size_t j = 0; j < M; j++) {
                uint64_t x = P0[j];
                int ok = 1;
                for (int u = 0; u < code_count; u++) {
                    int dist = popcount64(x ^ code[u]);
                    if (dist < d) { ok = 0; break; }
                }
                if (ok) {
                    valid_ids[valid_count++] = (int)j;
                }
            }
            if (valid_count == 0) {
                // failure: no valid extension
                break;
            }

            // 3) Compute feature histograms f_all[k] for each valid candidate
            for (int k = 0; k < valid_count; k++) {
                int *f = &f_all[k*B];
                memset(f, 0, B * sizeof(int));
                uint64_t x = P0[(size_t)valid_ids[k]];
                for (int u = 0; u < code_count; u++) {
                    int dist = popcount64(x ^ code[u]);
                    f[dist]++;
                }
            }

            // 4) Compute scores = Theta ⋅ f and softmax weights
            double max_s = -1e300;
            for (int k = 0; k < valid_count; k++) {
                double sc = 0.0;
                int *f = &f_all[k*B];
                for (int i = 0; i < B; i++) {
                    sc += Theta[i] * f[i];
                }
                score[k] = sc;
                if (sc > max_s) max_s = sc;
            }
            double sumexp = 0.0;
            for (int k = 0; k < valid_count; k++) {
                double ex = exp(beta * (score[k] - max_s));
                exp_score[k] = ex;
                sumexp += ex;
            }

            // 5) Compute expected feature vector E_f
            for (int i = 0; i < B; i++) E_f[i] = 0.0;
            for (int k = 0; k < valid_count; k++) {
                double p = exp_score[k] / sumexp;
                int *f = &f_all[k*B];
                for (int i = 0; i < B; i++) {
                    E_f[i] += p * f[i];
                }
            }

            // 6) Sample an action k according to softmax
            double r = rand_uniform();
            double csum = 0.0;
            int chosen_k = valid_count - 1;
            for (int k = 0; k < valid_count; k++) {
                csum += exp_score[k] / sumexp;
                if (r < csum) { chosen_k = k; break; }
            }
            int pick_j = valid_ids[chosen_k];
            code[code_count++] = P0[(size_t)pick_j];

            // 7) Store features and expected features for update
            for (int i = 0; i < B; i++) {
                features_mem[step_count*B + i] = f_all[chosen_k*B + i];
                Efeatures_mem[step_count*B + i] = E_f[i];
            }
            step_count++;
        }

        // 8) Check for success
        if (code_count == s) {
            // Print final code: one line per codeword
            for (int u = 0; u < s; u++) {
                uint64_t x = code[u];
                for (int b = 0; b < n; b++) {
                    int bit = (int)((x >> (n-1-b)) & 1ULL);
                    printf("%d%s", bit, (b+1<n) ? " " : "");
                }
                printf("\n");
            }
            fflush(stdout);
            return 0;
        }

        // 9) Episode finished with failure => policy update
        //    Use REINFORCE: reward = 1 on success, 0 on failure
        //    Here failure => we *subtract* the gradients
        for (int t = 0; t < step_count; t++) {
            int *f_t = &features_mem[t*B];
            double *E_t = &Efeatures_mem[t*B];
            for (int i = 0; i < B; i++) {
                double diff = (double)f_t[i] - E_t[i];
                Theta[i] -= alpha * diff;
            }
        }
        // Loop again for next episode
    }

    // unreachable
    return 0;
}