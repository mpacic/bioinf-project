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
    string targetHeader;
    int lineLength;

    std::string readReferenceGenome(std::string referenceGenomePath);
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

  if (!outputFile.is_open()) {
    std::cout << "Error: Failed to open output file" << std::endl;
    std::exit(1);
  }

  // read header
  getline(inputFile, targetHeader);
  string line;
  getline(inputFile, line);
  lineLength = stoi(line);
  outputFile << targetHeader << std::endl;

  // read lowercase positions
  cout << "Reading lowercase positions..." << endl;
  std::string lowercasePositions;
  getline(inputFile, lowercasePositions);
  std::istringstream iss(lowercasePositions);
  int start, length;
  std::vector<std::pair<int, int>> lpos;
  while (iss >> start >> length) {
    lpos.push_back(std::make_pair(start, length));
  }

  // read N positions
  cout << "Reading N positions..." << endl;
  std::string NPositions;
  getline(inputFile, NPositions);

  // read target sequence
  cout << "Reading target sequence..." << endl;
  std::string targetUncompressed = "";
  std::string targetSeq;
  int prevEnd = 0;
  while (getline(inputFile, targetSeq)) {
    if (targetSeq.find(',') != std::string::npos) {
      int start = stoi(targetSeq.substr(0, targetSeq.find(',')));
      int subseq_len = stoi(targetSeq.substr(targetSeq.find(',') + 1));
      targetUncompressed += referenceSeq.substr(prevEnd + start, subseq_len);
      prevEnd += start + subseq_len;
    } else {
      targetUncompressed += targetSeq;
    }
  }

  // to lowercase
  cout << "Updating lowercase positions..." << endl;
  int offset = 0;
  for (auto pos : lpos) {
    for (int i = pos.first; i < pos.first + pos.second; i++) {
      targetUncompressed[offset + i] = tolower(targetUncompressed[offset + i]);
    }
    offset += pos.first + pos.second;
  }

  for (int i = 0; i < targetUncompressed.length(); i += lineLength) {
    outputFile << targetUncompressed.substr(i, lineLength) << std::endl;
  }
}

std::string SCCGD::readReferenceGenome(std::string referenceGenomePath) {
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