# CPro1 configuration for constant weight binary codes.

Download the CPro1 code from the CPro1 repository, and follow the instructions there to set it up:
https://github.com/Constructive-Codes/CPro1

It was run using the OpenAI o4-mini model: MODEL="o4-mini-2025-04-16"

These are the two successful CPro1-produced programs, that each succeeded on one Open instance during the CPro1 run:

- `reinforcement-learning-cpro1.c` has the "reinforcement learning" approach from CPro1 that led to RSDH as described in the paper.  After compiling, this can be run as follows:

`./reinforcement-learning-cpro1 22 9 8 281 1000 0.08117297297297298 1.0`

This runs it on instance n=22 w=9 d=8 s=281, with random seed 1000, and hyperparameters 0.08117297297297298 1.0 (which are the hyperparameters as tuned by CPro1).

- `bit-wise-tabu-cpro1.c` is the original version of Bit-Wise Tabu from CPro1.  After compiling, it can be run like this:

`./bit-wise-tabu-cpro1 18 6 6 133 1000`

This runs it on instance n=18 w=6 d=6, with random seed 1000.

