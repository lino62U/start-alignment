#pragma once
#include <string>
#include <vector>
#include <utility>

std::vector<std::pair<std::string, std::string>>
align_with_center(const std::string& center, const std::vector<std::string>& sequences, int center_index) ;

std::vector<std::string>
merge_alignments(const std::vector<std::pair<std::string, std::string>>& alignments);

void write_to_file(
    const std::string& filename,
    const std::vector<std::vector<int>>& score_matrix,
    const std::vector<int>& sum_row,
    const std::vector<std::pair<std::string, std::string>>& alignments,
    const std::vector<std::string>& multiple_alignment
);
