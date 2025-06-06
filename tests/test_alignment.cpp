#include <gtest/gtest.h>
#include "star_alignment.hpp"  // Tu archivo con las funciones declaradas
#include "needleman_wunsch.hpp"

#include <fstream>
#include <string>
#include <vector>

TEST(StarAlignmentTest, FinalMultipleAlignmentIsCorrect) {
    std::vector<std::string> input = {
        "ATTGCCATT",
        "ATGGCCATT",
        "ATCCAATTTT",
        "ATCTTCTT",
        "ACTGACC"
    };

    std::string dummy_output = "unused_output.txt"; // No se usará en este test
    std::vector<std::string> result = star_alignment(input, dummy_output);

    std::vector<std::string> expected = {
        "ATTGCCATT--",
        "ATGGCCATT--",
        "ATC-CAATTTT",
        "ATCTTC-TT--",
        "ACTGACC----"
    };

    ASSERT_EQ(result.size(), expected.size()) << "Cantidad incorrecta de secuencias alineadas.";

    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(result[i], expected[i]) << "Diferencia en la línea " << i;
    }
}

TEST(StarAlignmentTest, IdenticalSequencesNoGaps) {
    std::vector<std::string> input = {
        "AGCT",
        "AGCT",
        "AGCT"
    };

    std::string dummy_output = "unused.txt";
    auto result = star_alignment(input, dummy_output);

    ASSERT_EQ(result.size(), input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        EXPECT_EQ(result[i], input[i]);
    }
}

TEST(StarAlignmentTest, ContainsEmptySequence) {
    std::vector<std::string> input = {
        "AGCT",
        "",
        "AG-T"
    };

    std::string dummy_output = "unused.txt";
    auto result = star_alignment(input, dummy_output);

    ASSERT_EQ(result.size(), input.size());
    for (const auto& seq : result) {
        EXPECT_EQ(seq.size(), result[0].size());  // Todas igual longitud
    }
}


TEST(NeedlemanWunschTest, PerfectMatch) {
    std::string seq1 = "GATTACA";
    std::string seq2 = "GATTACA";
    std::vector<std::pair<std::string, std::string>> alignment;

    int score = needleman_wunsch(seq1, seq2, alignment);

    EXPECT_EQ(score, seq1.length() * MATCH);
    ASSERT_EQ(alignment.size(), 1);
    EXPECT_EQ(alignment[0].first, seq1);
    EXPECT_EQ(alignment[0].second, seq2);
}

TEST(NeedlemanWunschTest, GapsAndMismatches) {
    std::string seq1 = "GATTACA";
    std::string seq2 = "GCATGCU";
    std::vector<std::pair<std::string, std::string>> alignment;

    int score = needleman_wunsch(seq1, seq2, alignment);

    ASSERT_EQ(alignment.size(), 1);
    std::string aligned1 = alignment[0].first;
    std::string aligned2 = alignment[0].second;

    // Mismo tamaño con gaps insertados
    ASSERT_EQ(aligned1.size(), aligned2.size());

    // Al menos alguna diferencia
    EXPECT_NE(aligned1, aligned2);
}

TEST(NeedlemanWunschTest, CompleteMismatch) {
    std::string seq1 = "AAAA";
    std::string seq2 = "TTTT";
    std::vector<std::pair<std::string, std::string>> alignment;

    int score = needleman_wunsch(seq1, seq2, alignment);

    ASSERT_EQ(alignment.size(), 1);
    EXPECT_EQ(alignment[0].first.size(), 4);
    EXPECT_EQ(alignment[0].second.size(), 4);
    EXPECT_EQ(score, 4 * MISMATCH);
}
