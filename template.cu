/*
Jakub Frąc 298 795

Hamming One solved in O(n*k*logn)

It could be solved in O(nk) if instead of using binsearch I used 3 small hashing constants for modulo
and then calculate 3 different hashes for single word. Chance of a collision of such a technique
would be p_1 * p_2 * p_3 where p_i is our i-th prime hash constant.

Both binsearch and calculating prefix sums (or prefix hashes) could be implemented more efficiently
with use of GPU-scan algorithm.
*/

#include <math.h>
#include <cuda_runtime.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <chrono>
#include <numeric>

#define N 100*1001
#define M 1050

using namespace std;

// Hashing prime number constant
__constant__ int d_p = 2137;


// Function to calculate array with hash prefixes for given word
__device__ void calculatePref(unsigned long long *pref, unsigned long long *pPow, const char* word, int strlen) {
	pref[0] = word[0];

	for (int i = 1; word[i] != '\0'; i++)
		pref[i] = pref[i - 1] + word[i] * pPow[i];
}

// Performes binsearch on array of hashes created during preprocessing. 
// If hash given as parameter is found in the array - it returns its index.
__device__ int fIndex(unsigned long long hash, unsigned long long *idx, int n) {
	int l = 0, r = n - 1;
	while (l != r) {
		int mid = (l + r) / 2;
		if (idx[mid] < hash)
			l = mid + 1;
		else
			r = mid;
	}
	if (idx[r] != hash)
		return -1;
	return r;
}

// Function to find index in the hashmap.
__device__ int find(unsigned long long hash, int *hMap, unsigned long long *idx, int n) {
	int id = fIndex(hash, idx, n);
	if (id == -1)
		return -1;

	return hMap[id];
}

// Function to add hash to the hashmap.
__device__ void add(unsigned long long hash, int index, int *hMap, unsigned long long *idx, int n) {
	int id = fIndex(hash, idx, n);
	atomicExch(hMap + id, index);
}

// Each thread receives it's word. First it calculates prefix hashes for given word. It then adds the word to the hashmap.
// Then it starts iterating throughout each letter of a word, calculating hash of word with i-th bit flipped (with hamming dist 1
// It then checks if given hash exists in the hashmap - if it does it adds it as a found pair.
__global__ void Kernel(int *n, char* words, unsigned long long *pPow, int *hMap, int *wlens, int *results, unsigned long long *idx) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if (tid >= *n)
		return;

	unsigned long long pref[M];
	unsigned long long curHash = 0;
	int wlen = wlens[tid];

	const char *word = (words + tid * M);

	calculatePref(pref, pPow, word, wlen);

	curHash = 0;

	add(pref[wlen - 1], tid, hMap, idx, *n);

	for (int i = wlen - 1; i >= 0; i--) {
		unsigned long long prefHash = 0;
		if (i > 0)
			prefHash = pref[i - 1];

		unsigned long long hammHash = prefHash + ((word[i] - '0' + 1) % 2 + '0') * pPow[i] + curHash * pPow[i + 1];


		int index = find(hammHash, hMap, idx, *n);

		if (index != -1) {
			results[tid] = index;
		}

		curHash = curHash * d_p + word[i];
	}
}

const int h_p = 2137;

unsigned long long h_pPow[M], h_idx[N], h_wHash[N];
int h_hMap[N], h_results[N], h_wlens[N];
char h_words[N][M], h_temp[N][M];

// Calculating hash of a word using honer's algorithm
unsigned long long calcHash(char *word) {
	unsigned long long hash = 0;

	for (int i = strlen(word) - 1; i >= 0; i--)
		hash = hash * h_p + word[i];

	return hash;
}

// Pretty self explanatory
void readInput(int *h_n) {
	ifstream myfile("test.big");

	if (myfile.is_open())
	{
		string line;
		getline(myfile, line);
		*h_n = stoi(line);
		int n = *h_n;

		for (int i = 0; i < n; i++)
			myfile.getline(h_words[i], M);

		myfile.close();
	}
	else {
		cout << "Unable to open file";
		exit(0);
	}
}
// Calculates hash of every word, then sorts the words by their hash (hash with smaller value first)
// We can then do binsearch on these hashes and grab index the word associated with it.
// We also calculate word lengths and array containing next powers of prime number p (declared above)
void preprocessInput(int h_n) {
	for (int i = 0; i < h_n; i++)
		h_wHash[i] = calcHash(h_words[i]);

	iota(h_idx, h_idx + h_n, 0);
	sort(h_idx, h_idx + h_n, [](unsigned long long l, unsigned long long r) { return h_wHash[l] < h_wHash[r]; });

	for (int i = 0; i < h_n; i++)
		memmove(h_temp[i], h_words[h_idx[i]], M);

	for (int i = 0; i < h_n; i++) {
		memmove(h_words[i], h_temp[i], M);
		h_idx[i] = h_wHash[h_idx[i]];
	}

	for (int i = 0; i < h_n; i++)
		h_wlens[i] = strlen(h_words[i]);

	h_pPow[0] = 1;
	for (int i = 1; i < M; i++)
		h_pPow[i] = h_pPow[i - 1] * h_p;

	for (int i = 0; i < N; i++)
		h_hMap[i] = -1;
}
// Pretty self explanatory
void prepareMemory(int h_n, int *&d_results, int *&d_wlens, int *&d_hMap, int *&d_n, unsigned long long *&d_pPow, unsigned long long *&d_idx, char *&d_words) {
	memset(h_results, -1, h_n * sizeof(int));

	cudaMalloc(&d_results, h_n * sizeof(int));
	cudaMalloc(&d_wlens, h_n * sizeof(int));
	cudaMalloc(&d_words, h_n * M * sizeof(char *));
	cudaMalloc(&d_n, sizeof(int));
	cudaMalloc(&d_pPow, M * sizeof(long long));
	cudaMalloc(&d_hMap, h_n * sizeof(int));
	cudaMalloc(&d_idx, h_n * sizeof(long long));

	cudaMemcpy(d_pPow, h_pPow, h_n * sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(d_hMap, h_hMap, h_n * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_idx, h_idx, h_n * sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(d_n, &h_n, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_results, h_results, h_n * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_words, h_words, h_n * M * sizeof(char), cudaMemcpyHostToDevice);
	cudaMemcpy(d_wlens, h_wlens, h_n * sizeof(int), cudaMemcpyHostToDevice);
}

// Calculates number of threads, starts kernel and copies result to cpu memory
void runAlgorithm(int h_n, int *d_n, int *d_hMap, int *d_wlens, int *d_results, unsigned long long *d_pPow, unsigned long long *d_idx, char *d_words) {
	int threadsPerBlock = 256;
	int blocksPerGrid = (h_n + threadsPerBlock - 1) / threadsPerBlock;

	Kernel << <blocksPerGrid, threadsPerBlock >> > (d_n, d_words, d_pPow, d_hMap, d_wlens, d_results, d_idx);

	int r = cudaDeviceSynchronize();
	if (r != cudaSuccess) {
		cout << "error" << " " << r << endl;
		exit(0);
	}
	cudaMemcpy(h_results, d_results, h_n * sizeof(int), cudaMemcpyDeviceToHost);
}


// Pretty self explanatory
void freeMemory(unsigned long long *d_pPow, unsigned long long *d_idx, int *d_hMap, int *d_wlens, int *d_n, int *d_results, char *d_words) {
	cudaFree(d_pPow);
	cudaFree(d_idx);
	cudaFree(d_hMap);
	cudaFree(d_wlens);
	cudaFree(d_n);
	cudaFree(d_results);
	cudaFree(d_words);
}

int main() {
	unsigned long long *d_pPow, *d_idx;
	int h_n, *d_hMap, *d_n, *d_wlens, *d_results;
	char *d_words;

	readInput(&h_n);

	auto gStart = std::chrono::high_resolution_clock::now();
	auto start = std::chrono::high_resolution_clock::now();
	preprocessInput(h_n);
	auto end = std::chrono::high_resolution_clock::now();

	cout << "Time spent for preprocessing input: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6 << endl;

	start = std::chrono::high_resolution_clock::now();
	prepareMemory(h_n, d_results, d_wlens, d_hMap, d_n, d_pPow, d_idx, d_words);
	end = std::chrono::high_resolution_clock::now();
	cout << "Time spent for preparation of memory: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6 << endl;

	start = std::chrono::high_resolution_clock::now();
	runAlgorithm(h_n, d_n, d_hMap, d_wlens, d_results, d_pPow, d_idx, d_words);
	end = std::chrono::high_resolution_clock::now();
	cout << "Time spent for running the algorithm: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6 << endl;

	int ans = 0;
	for (int i = 0; i < h_n; i++) {
		if (h_results[i] > i) {
			ans++;
		}
	}

	auto gEnd = std::chrono::high_resolution_clock::now();
	cout << "Overall time spent: " << std::chrono::duration_cast<std::chrono::microseconds>(gEnd - gStart).count() / 1e6 << endl;

	cout << "Answer: " << ans << endl;
	freeMemory(d_pPow, d_idx, d_hMap, d_wlens, d_n, d_results, d_words);
}

