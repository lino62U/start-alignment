#include "align_utils.hpp"
#include "needleman_wunsch.hpp"
#include <fstream>
#include <algorithm>

std::vector<std::pair<std::string, std::string>>
align_with_center(const std::string& center, const std::vector<std::string>& sequences, int center_index) {
    std::vector<std::pair<std::string, std::string>> alignments;
    for (size_t i = 0; i < sequences.size(); ++i) {
        if ((int)i == center_index) {
            // Devuelve la secuencia central alineada consigo misma sin gaps
            alignments.push_back({center, center});
        } else {
            std::vector<std::pair<std::string, std::string>> result;
            needleman_wunsch(center, sequences[i], result);
            alignments.push_back(result[0]);
        }
    }
    return alignments;
}

std::vector<std::string> merge_alignments(const std::vector<std::pair<std::string, std::string>>& alignments) {
    std::string base = alignments[0].first;
    std::vector<std::string> aligned;

    for (size_t i = 0; i < alignments.size(); ++i)
        aligned.push_back(alignments[i].second);

    aligned.insert(aligned.begin(), base);
    size_t max_len = 0;
    for (const auto& s : aligned)
        max_len = std::max(max_len, s.size());

    for (auto& s : aligned)
        s += std::string(max_len - s.size(), '-');

    return aligned;
}

void write_to_file(const std::string& filename,
                   const std::vector<std::vector<int>>& score_matrix,
                   const std::vector<int>& sum_row,
                   const std::vector<std::pair<std::string, std::string>>& alignments,
                   const std::vector<std::string>& multiple_alignment) {
    std::ofstream f(filename);
    int n = score_matrix.size();

    f << "Matriz de Scores:\n";
    for (int i = 0; i < n; ++i) {
        for (int val : score_matrix[i])
            f << val << "  ";
        f << "| suma: " << sum_row[i] << "\n";
    }

    f << "\nAlineamientos pareados:\n";
    for (size_t i = 0; i < alignments.size(); ++i)
        f << "\nAlineamiento " << i + 1 << ":\n" << alignments[i].first << "\n" << alignments[i].second << "\n";

    f << "\nAlineamiento Múltiple Final:\n";
    for (const auto& s : multiple_alignment)
        f << s << "\n";
}
