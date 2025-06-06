#include "needleman_wunsch.hpp"
#include <vector>
#include <string>
#include <algorithm>

int needleman_wunsch(const std::string& seq1, const std::string& seq2, std::vector<std::pair<std::string, std::string>>& alignment) {
    size_t m = seq1.size(), n = seq2.size();
    std::vector<std::vector<int>> score(m + 1, std::vector<int>(n + 1, 0));

    for (size_t i = 0; i <= m; ++i) score[i][0] = i * GAP;
    for (size_t j = 0; j <= n; ++j) score[0][j] = j * GAP;

    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            int match = score[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH : MISMATCH);
            int del = score[i - 1][j] + GAP;
            int ins = score[i][j - 1] + GAP;
            score[i][j] = std::max({match, del, ins});
        }
    }

    std::string aln1, aln2;
    size_t i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && score[i][j] == score[i - 1][j] + GAP) {
            aln1 = seq1[i - 1] + aln1;
            aln2 = '-' + aln2;
            --i;
        } else if (j > 0 && score[i][j] == score[i][j - 1] + GAP) {
            aln1 = '-' + aln1;
            aln2 = seq2[j - 1] + aln2;
            --j;
        } else {
            aln1 = seq1[i - 1] + aln1;
            aln2 = seq2[j - 1] + aln2;
            --i; --j;
        }
    }

    alignment.push_back({aln1, aln2});
    return score[m][n];
}
