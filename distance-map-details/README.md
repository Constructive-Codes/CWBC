# Details for the distance map code: n=25 d=8 w=11 code with 1666 codewords

The distance map figure in the paper is based on an n=25 d=8 w=11 codeword code with 1666 codewords, obtained by RSDH.  This code is in `code-25-8-11-1666.txt`.

The distance map image is in distmap.png.

The red/marooon/black upper-left block corresponds to the first 960 codewords.  Each of these has the 3rd column equal to 1.  Taking these first 960 codewords and removing the third column gives `code960-from-1666.txt`.

The file `code960-24-8-10.txt` is the 960-codeword code with n=24 d=8 n=10, derived from the Golay code, linked from [https://aeb.win.tue.nl/codes/Andw.html](https://aeb.win.tue.nl/codes/Andw.html).  The specific code is [https://aeb.win.tue.nl/codes/cwc/d8/a24.8.10.960H](https://aeb.win.tue.nl/codes/cwc/d8/a24.8.10.960H) which here has been converted to binary.

Applying this permutation (numbering from 0):  
`perm = [0, 1, 2, 21, 11, 13, 5, 8, 17, 10, 20, 7, 14, 6, 12, 22, 15, 19, 23, 9, 4, 16, 3, 18]`  
to the columns of `code960-from-1666.txt` generates the same codewords as `code960-24-8-10.txt`, showing they are isomorphic (as noted in the paper).  The short Python script `permute.py` verifies this.
