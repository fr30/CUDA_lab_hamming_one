# CUDA_lab_hamming_one

Solution for Hamming one problem. It works in O(n*l*logn) time complexity.

Few thing worth mentioning:
1. Hashing function doesn't use modulo operation, since I calculate remainder by overflowing long long integer.
2. Program returns results only for different pairs, thus if you need to get the exact answer you need to parse the output yourself.
3. In order to better understand the code It's probably easier to figure out CPU version first, then GPU.
