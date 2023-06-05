#include <iostream>
#include <filesystem>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

class SCCGD {
  public:
    SCCGD(std::string referenceGenomePath, std::string inputFilePath,
        std::string outputDirPath)
      : referenceGenomePath(referenceGenomePath),
        inputFilePath(inputFilePath),
        outputDirPath(outputDirPath){};
  
    void run();
  
  private:
    const string referenceGenomePath;
    const string inputFilePath;
    const string outputDirPath;
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

  std::string referenceSeq = readReferenceGenome(referenceGenomePath);

  std::ofstream interimFile(outputDirPath + "/interim.txt");

  std::ofstream outputFile(outputDirPath + "/output.txt");

  // read input file
  std::ifstream inputFile(inputFilePath);
  // check file opened successfully
  if (!inputFile.is_open()) {
    std::cout << "Error: Failed to open input file" << std::endl;
    std::exit(1);
  }
  std::string line;

  // read header
  std::string targetHeader;
  getline(inputFile, targetHeader);
  outputFile << targetHeader << std::endl;

  // read lowercase positions
  std::string lowercasePositions;
  getline(inputFile, lowercasePositions);
  std::istringstream iss(lowercasePositions);
  int start, length;
  std::vector<std::pair<int, int>> lpos;
  while (iss >> start >> length) {
    lpos.push_back(std::make_pair(start, length));
  }

  // read target sequence
  std::string targetSeq;
  getline(inputFile, targetSeq);


}

std::string readReferenceGenome(std::string referenceGenomePath) {
  ifstream referenceGenomeFile(referenceGenomePath);
  // check file opened successfully
  if (!referenceGenomeFile.is_open()) {
    std::cout << "Error: Failed to open referencreade genome file" << std::endl;
    std::exit(1);
  }

  string referenceGenome = "";
  string line;
  // skip first line
  getline(referenceGenomeFile, line);

  while (getline(referenceGenomeFile, line)) {
    referenceGenome += line;
  }

  return referenceGenome;
}