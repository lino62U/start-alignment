#include "fasta_reader.hpp"
#include <fstream>

std::vector<std::string> leer_fasta(const std::string& filepath) {
    std::ifstream file(filepath);
    std::vector<std::string> secuencias;
    std::string line, current;

    while (getline(file, line)) {
        // Eliminar caracteres de retorno de carro (Windows)
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current.empty()) {
                secuencias.push_back(current);
                current.clear();
            }
        } else {
            // Concatenar sin espacios ni saltos extraños
            for (char c : line) {
                if (!isspace(c)) current += c;
            }
        }
    }

    if (!current.empty()) {
        secuencias.push_back(current);
    }

    return secuencias;
}
