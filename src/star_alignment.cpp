#include "star_alignment.hpp"
#include "score_matrix.hpp"
#include "align_utils.hpp"
#include <iostream>

std::vector<std::string> star_alignment(const std::vector<std::string>& sequences, const std::string& output_file) {
    
    auto [score_matrix, sum_row] = build_score_matrix(sequences);
    int center_index = print_score_table(sequences, score_matrix, sum_row);
    auto center_seq = sequences[center_index];

    auto pairwise_aligns = align_with_center(center_seq, sequences, center_index);

    std::vector<std::string> multiple;
    multiple.push_back(pairwise_aligns[center_index].first); // Secuencia central con gaps

    for (int i = 0; i < (int)sequences.size(); ++i) {
        if (i == center_index) continue;
        multiple.push_back(pairwise_aligns[i].second); // Secuencias alineadas con gaps
    }

    // Igualar longitud de todas las secuencias agregando gaps al final
    size_t max_len = 0;
    for (const auto& s : multiple) {
        if (s.size() > max_len) max_len = s.size();
    }

    for (auto& s : multiple) {
        if (s.size() < max_len) {
            s += std::string(max_len - s.size(), '-');
        }
    }

    std::cout << "Alineamiento múltiple final:\n";
    for (const auto& s : multiple) std::cout << s << "\n";

    write_to_file(output_file, score_matrix, sum_row, pairwise_aligns, multiple);

    return multiple;
}