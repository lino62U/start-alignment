#pragma once
#include <vector>
#include <string>

std::tuple<std::vector<std::vector<int>>, std::vector<int>>
build_score_matrix(const std::vector<std::string>& sequences);

int print_score_table(
    const std::vector<std::string>& sequences,
    const std::vector<std::vector<int>>& score_matrix,
    const std::vector<int>& sum_row
);
