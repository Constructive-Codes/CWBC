# CPro1 configuration for constant weight binary codes.

Download the CPro1 code from the CPro1 repository, and follow the instructions there to set it up:
https://github.com/Constructive-Codes/CPro1

It was run using the OpenAI o4-mini model: MODEL="o4-mini-2025-04-16"

Put the following two files in the directory where you will run, replacing the versions you got from the CPro1 repository:

`conf.py` here has the setup for constant weight binary codes.  It matches defaults, except note that:
- It generates only 200 programs per run, rather than the default of 1000 (STRATEGY_REPS = 10 instead of 50).
- The code optimization step is turned off ("prog_opt_prompt_eval": None), to save on tokens.
- PROMPT_NOVELTY = True.  Half the runs in the paper were run with this True, half with it False.
- PROMPT_NUM_HYPERPARAMETERS = -1.  This turns off the prompt limiting hyperparameters, matching the older behavior of CPro1 that was used for the constant weight binary codes runs.

These are the two successful CPro1-produced programs; each one succeeded on one Open instance during the CPro1 run:

- `reinforcement-learning-cpro1.c` has the "reinforcement learning" approach from CPro1 that led to RSDH as described in the paper.  After compiling, this can be run as follows:

`./reinforcement-learning-cpro1 22 9 8 281 1000 0.08117297297297298 1.0`

This runs it on instance n=22 w=9 d=8 s=281, with random seed 1000, and hyperparameters 0.08117297297297298 1.0 (which are the hyperparameters as tuned by CPro1).

- `bit-wise-tabu-cpro1.c` is the original version of Bit-Wise Tabu from CPro1.  After compiling, it can be run like this:

`./bit-wise-tabu-cpro1 18 6 6 133 1000`

This runs it on instance n=18 w=6 d=6, with random seed 1000.

