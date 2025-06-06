#pragma once
#include <string>
#include <vector>
#include <tuple>

constexpr int MATCH = 1;
constexpr int MISMATCH = -1;
constexpr int GAP = -2;

int needleman_wunsch(
    const std::string& seq1,
    const std::string& seq2,
    std::vector<std::pair<std::string, std::string>>& alignment
);
