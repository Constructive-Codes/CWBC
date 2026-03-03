# Data and code for "Automated Discovery of Improved Constant Weight Binary Codes"

This repository has the data and code for:

C.D. Rosin. (March 2026) “Automated Discovery of Improved Constant Weight Binary Codes.” [https://arxiv.org/pdf/2603.00174](https://arxiv.org/pdf/2603.00174)

The directory **codes/** has the new constant weight binary codes.  Files are named code-n-d-w-s.txt where n is the number of bits, d is the distance, w is the weight, and s is the number of codewords.  The codewords are one per line of the file.

The C programs that generated these codes are below.  Under Linux, they can be compiled with (for example):
`gcc -O3 -march=native -o rsdh rsdh.c -lm`
Each program is single-threaded, but multiple runs can be started on a single machine and proceed simultaneously without interference.

- `bit-swap-tabu.c`: Bit-Swap Tabu Search.  Command-line parameters are `n w d s seed tmin tmax max_no_improve_restart`.  Runs for the paper used tmin=5, tmax=15, max_no_improve_restart=1000000.  Example:  
`./bit-swap-tabu 18 6 6 134 123 5 15 1000000 > log`

- `rsdh.c` - Random-Score Distance Histograms (RSDH).  Command-line parameters are `n w d s seed threshold`.  Use threshold\<s.  When more than threshold codewords are found, the code is output into the log.  The run can then continue, and find codes up with up to s codewords.  Example:  
`./rsdh 22 8 9 300 123 280 > log`

- `s-rsdh.c` - Sliced Random-Score Distance Histograms (S-RSDH).  Command-line parameters are `n w d s seed threshold slicebits trials`.  slicebits is the number of bits that define the slice; runs in the paper use 1-3.  trials is the number of completions after the first slice; runs in the paper use trials=1000.  Example:  
`./s-rsdh 24 11 8 1400 123 1288 1 1000 > log`

- `ms-rsdh.c` - Multi-Sliced Random-Score Distance Histograms (MS-RSDH).  Command-line parameters are the same as for s-rsdh.c.  Use slicebits 2-3 (slicebits 1 is equivalent to s-rsdh with slicebgits 1).  Runs in the paper use trials=1000.  Example:  
`./ms-rsdh 26 8 8 800 123 760 2 1000 > log`

The above programs are based on the automated results from CPro1, but have been further developed by hand, as described in the paper.

The directory **CPro1-for-CWBC/** has the configuration files for CPro1, instructions for running it on this problem, and the original raw code generated automatically by CPro1.

The directory **distance-map-details/** has the code and details for the distance map figure in the paper.
