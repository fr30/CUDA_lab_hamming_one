//#include <iostream>
//#include <vector>
//#include <algorithm>
//#include <string>
//#include <fstream>
//#include <numeric>
//#include <chrono>
//
//#define maxHash 1000*1001
//#define N 100*1001
//#define M 3*1000
//
//using namespace std;
//
//const unsigned long long p = 2137;
//
//int n;
//int hashmap[maxHash];
//unsigned long long wHash[N];
//unsigned long long idx[N];
//unsigned long long pref[M];
//unsigned long long pPow[M];
//string words[N];
//string temp[N];
//
//void calculatePref(unsigned long long *pref, unsigned long long *pPow, const string& word) {
//	memset(pref, 0, word.size() + 1);
//	pref[0] = word[0];
//
//	for (int i = 1; i < word.size(); i++)
//		pref[i] = (pref[i - 1] + word[i] * pPow[i]);
//}
//
//int fIndex(unsigned long long hash) {
//	int l = 0, r = n - 1;
//	while (l != r) {
//		int mid = (l + r) / 2;
//		if (idx[mid] < hash)
//			l = mid + 1;
//		else
//			r = mid;
//	}
//	if (idx[r] != hash)
//		return -1;
//	return r;
//}
//
//int find(unsigned long long hash) {
//	int id = fIndex(hash);
//	if (id == -1)
//		return -1;
//
//	return hashmap[id];
//}
//
//void add(unsigned long long hash, int index) {
//	int id = fIndex(hash);
//	hashmap[id] = index;
//}
//
//unsigned long long calcHash(string word) {
//	unsigned long long hash = 0;
//
//	for (int i = word.size() - 1; i >= 0; i--)
//		hash = hash * p + word[i];
//
//	return hash;
//}
//
//int main() {
//	vector <pair <string, string>> result;
//
//	ifstream myfile("test.big");
//
//	if (myfile.is_open())
//	{
//		string line;
//		getline(myfile, line);
//		n = stoi(line);
//		
//		for (int i = 0; i < n; i++) {
//			getline(myfile, words[i]);
//			wHash[i] = calcHash(words[i]);
//		}
//	}
//	auto start = std::chrono::high_resolution_clock::now();
//
//	iota(idx, idx + n, 0);
//	sort(idx, idx + n, [](unsigned long long l, unsigned long long r) { return wHash[l] < wHash[r]; });
//
//	for (int i = 0; i < n; i++)
//		temp[i] = move(words[idx[i]]);
//	for (int i = 0; i < n; i++) {
//		words[i] = move(temp[i]);
//		idx[i] = wHash[idx[i]];
//	}
//
//	pPow[0] = 1;
//	for (int i = 1; i < M; i++)
//		pPow[i] = (pPow[i - 1] * p);
//
//	for (int i = 0; i <= maxHash; i++)
//		hashmap[i] = -1;
//
//	unsigned long long curHash = 0;
//
//	for (int w = 0; w < n; w++) {
//		const string& word = words[w];
//		calculatePref(pref, pPow, word);
//		curHash = 0;
//
//		for (int i = word.size() - 1; i >= 0; i--) {
//			unsigned long long prefHash = 0;
//
//			if (i > 0)
//				prefHash = pref[i - 1];
//
//			unsigned long long hammHash = (prefHash + ((word[i] - '0' + 1) % 2 + '0') * pPow[i] + curHash * pPow[i + 1]);
//
//
//			int index = find(hammHash);
//
//			if (index != -1) {
//				result.push_back(make_pair(words[index], word));
//			}
//
//			curHash = (curHash * p + word[i]);
//		}
//		add(pref[word.size() - 1], w);
//	}
//
//	sort(result.begin(), result.end());
//	auto last = unique(result.begin(), result.end());
//	result.erase(last, result.end());
//
//	cout << result.size() << endl;
//	//sort(result.begin(), result.end());
//	//for (auto r : result)
//	//	cout << "(" << r.first << ", " << r.second << ")" << endl;
//	auto end = std::chrono::high_resolution_clock::now();
//	cout << "Overall time spent: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1e6 << endl;
//}