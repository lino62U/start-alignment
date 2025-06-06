#include "score_matrix.hpp"
#include "needleman_wunsch.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>

std::tuple<std::vector<std::vector<int>>, std::vector<int>>
build_score_matrix(const std::vector<std::string>& sequences) {
    int n = sequences.size();
    std::vector<std::vector<int>> score_matrix(n, std::vector<int>(n, 0));
    std::vector<int> sum_row(n, 0);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j) {
                std::vector<std::pair<std::string, std::string>> align;
                score_matrix[i][j] = needleman_wunsch(sequences[i], sequences[j], align);
                sum_row[i] += score_matrix[i][j];
            }

    return {score_matrix, sum_row};
}

int print_score_table(const std::vector<std::string>& sequences,
                      const std::vector<std::vector<int>>& score_matrix,
                      const std::vector<int>& sum_row) {


    int n = sequences.size();
    int width = 6;

    std::cout << "Matriz de Alineamientos Globales (Needleman-Wunsch):\n";
    
    // Imprimir encabezados
    std::cout << "      ";
    for (int i = 0; i < n; ++i) {
        std::cout << std::setw(width) << ("S" + std::to_string(i)) << "  ";
    }
    std::cout << std::setw(width) << "Suma" << "\n\n";

    // Imprimir filas con scores y suma por fila
    for (int i = 0; i < n; ++i) {
        std::cout << "S" << i << " | ";
        for (int j = 0; j < n; ++j) {
            std::cout << std::setw(width) << score_matrix[i][j] << "  ";
        }
        std::cout << "| " << std::setw(width) << sum_row[i] << "\n";
    }

    // Calcular suma por columnas
    std::vector<int> col_sum(n, 0);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            col_sum[j] += score_matrix[i][j];
        }
    }

    // Imprimir suma por columnas
    std::cout << "Suma| ";
    for (int j = 0; j < n; ++j) {
        std::cout << std::setw(width) << col_sum[j] << "  ";
    }

    int total = 0;
    for (int val : col_sum) total += val;

    std::cout << "| " << std::setw(width) << total << "\n\n";

    // Secuencia central
    int center = std::distance(sum_row.begin(), std::max_element(sum_row.begin(), sum_row.end()));
    std::cout << "🟊 Secuencia central: S" << center << " -> " << sequences[center] << "\n\n";

    return center;
}
