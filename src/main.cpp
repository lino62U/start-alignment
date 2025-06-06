#include <iostream>
#include <vector>
#include <string>
#include "star_alignment.hpp"
#include "fasta_reader.hpp"


int main() {
    std::vector<std::string> secuencias = {
        "ATTGCCATT",
        "ATGGCCATT",
        "ATCCAATTTT",
        "ATCTTCTT",
        "ACTGACC"
    };

    star_alignment(secuencias);

    return 0;
}




/*
int main(int argc, char* argv[]) {
   

    std::string input_fasta = "data/frataxin_mamiferos.fasta";

    std::string output_file = (argc >= 3) ? argv[2] : "alineamiento_estrella.txt";

    std::cout << "Leyendo secuencias desde: " << input_fasta << "\n";
    auto sequences = leer_fasta(input_fasta);

    if (sequences.empty()) {
        std::cerr << "No se encontraron secuencias en el archivo.\n";
        return 1;
    }

     std::cout << "Secuencias leídas:\n";
    for (size_t i = 0; i < sequences.size(); ++i) {
        std::cout << "Secuencia " << i + 1 << ": " << sequences[i] << "\n";
    }


    std::cout << "Secuencias leídas: " << sequences.size() << "\n";
    star_alignment(sequences, output_file);

    std::cout << "Resultado guardado en: " << output_file << "\n";

    return 0;
}
*/