#include <iostream>
#include <filesystem>
#include <string>

using namespace std;

class SCCGD {
  public:
    SCCGD(const char* reference_genome_file, const char* input_file,
          const char* output_directory)
        : reference_genome_file_(reference_genome_file),
          input_file_(input_file),
          output_directory_(output_directory) {}
  
    void run();
  
  private:
    const string reference_genome_file_;
    const string input_file_;
    const string output_directory_;
};

int main(int argc, char** argv) {
  // check number of arguments
  if (argc < 4) {
    std::cout << "Usage: " << argv[0]
              << " <reference genome file> <input file> <output_directory>"
              << std::endl;
    return 1;
  }

  // check reference genome file exists
  if (!std::filesystem::exists(argv[1])) {
    std::cout << "Error: Reference genome file does not exist: " << argv[1]
              << std::endl;
    return 1;
  }

  // check input file exists
  if (!filesystem::exists(argv[2])) {
    std::cout << "Error: Input file does not exist: " << argv[2] << std::endl;
    return 1;
  }

  // check output directory exists
  if (!filesystem::exists(argv[3])) {
    std::cout << "Error: Output directory does not exist: " << argv[3]
              << std::endl;
    return 1;
  }

  SCCGD sccgd(argv[1], argv[2], argv[3]);

  sccgd.run();
  return 0;
}

void SCCGD::run() {
  std::cout << "Running SCCGD" << std::endl;


}