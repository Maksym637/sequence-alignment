/*
This is an implementation of the global alignment algorithm 
of the DNA molecule.
DNA consists with 4 nucleotides :
- adenine (A)
- cytosine (C)
- guanine (G)
- thymine (T)
*/
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// S(xi, yj) = 1 if match or 0 if mismatch
const int MATCH = 1;
const int MISMATCH = 0;

// d = 0 is our penalty function
const int GAP = 0;

int maximum(int a, int b, int c) {
    return max(a, max(b, c));
}

int globalAlignmentScore(const string& s1, const string& s2, vector<vector<int>>& dp) {
    int n = s1.size(), m = s2.size();

    dp[0][0] = 0;

    for (int i = 1; i <= n; ++i) {
        dp[i][0] = i * GAP;
    }

    for (int j = 1; j <= m; ++j) {
        dp[0][j] = j * GAP;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            /*
            F(i, j) = max {
                F(i - 1, j - 1) + S(xi, yj),
                F(i - 1, j) + d,
                F(i, j - 1) + d
            }
            */
            int match = dp[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? MATCH : MISMATCH);
            int deleteChar = dp[i - 1][j] + GAP;
            int insertChar = dp[i][j - 1] + GAP;
            dp[i][j] = maximum(match, deleteChar, insertChar);
        }
    }

    return dp[n][m];
}

void traceback(const string& s1, const string& s2, const vector<vector<int>>& dp, string& align1, string& align2) {
    /*
    Direct sequence alignment
    */
    int n = s1.size(), m = s2.size();
    int i = n, j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? MATCH : MISMATCH)) {
            align1 = s1[i - 1] + align1;
            align2 = s2[j - 1] + align2;
            --i;
            --j;
        } else if (i > 0 && dp[i][j] == dp[i - 1][j] + GAP) {
            align1 = s1[i - 1] + align1;
            align2 = "-" + align2;
            --i;
        } else {
            align1 = "-" + align1;
            align2 = s2[j - 1] + align2;
            --j;
        }
    }
}

int main() {
    string sequence1 = "ATGGCCTC", sequence2 = "ACGGCTC";
    string alignment1, alignment2;

    vector<vector<int>> dp(sequence1.size() + 1, vector<int>(sequence2.size() + 1, 0));

    int score = globalAlignmentScore(sequence1, sequence2, dp);
    traceback(sequence1, sequence2, dp, alignment1, alignment2);

    cout << "------------------------------------" << endl;
    cout << "Total number of matches (score) : " << score << endl;
    cout << "------------------------------------" << endl;
    cout << alignment1 << endl << alignment2 << endl;
    cout << "------------------------------------" << endl;

    return 0;
}
