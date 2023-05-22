#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>
#include <unordered_map>
using namespace std;

class SCCGC {

 public:
  SCCGC(
    std::string referenceGenomePath,
    std::string inputFilePath,
    std::string outputDirPath
  ) : referenceGenomePath(referenceGenomePath), inputFilePath(inputFilePath), outputDirPath(outputDirPath) {};
  ~SCCGC(){};
  void run();
 private:
  std::string referenceGenomePath;
  std::string inputFilePath;
  std::string outputDirPath;
  std::string referenceSeq;
  int kmer_size;
};

int main(int argc, char **argv) {

  // check number of arguments
  if(argc < 4) {
    std::cout << "Usage: " << argv[0] << " <reference genome file> <input file> <output_directory>" << std::endl;
    return 1;
  }

  // check reference genome file exists
  if(!filesystem::exists(argv[1])) {
    std::cout << "Error: Reference genome file does not exist: " << argv[1] << std::endl;
    return 1;
  }

  // check input file exists
  if(!filesystem::exists(argv[2])) {
    std::cout << "Error: Input file does not exist: " << argv[2] << std::endl;
    return 1;
  }

  // check output directory exists
  if(!filesystem::exists(argv[3])) {
    std::cout << "Error: Output directory does not exist: " << argv[3] << std::endl;
    return 1;
  }

  SCCGC sccgc(argv[1], argv[2], argv[3]);

  sccgc.run();
}

std::string parseReferenceGenome(std::string referenceGenomePath) {
  ifstream referenceGenomeFile(referenceGenomePath);
  // check file opened successfully
  if(!referenceGenomeFile.is_open()) {
    std::cout << "Error: Failed to open reference genome file" << std::endl;
  }

  string referenceGenome = "";
  string line;
  // skip first line
  getline(referenceGenomeFile, line);

  while(getline(referenceGenomeFile, line)) {
    for(int i = 0; i < line.length(); i++) {
      char c = toupper(line[i]);
      if (c != 'N') {
        referenceGenome += toupper(line[i]);
      }
    }
  }
  return referenceGenome;
}

std::string parseTargetGenome(std::string targetGenomePath) {
  ifstream targetGenomeFile(targetGenomePath);
  // check file opened successfully
  if(!targetGenomeFile.is_open()) {
    std::cout << "Error: Failed to open target genome file" << std::endl;
  }

  string targetGenome = "";
  string line;
  // skip first line
  getline(targetGenomeFile, line);

  while(getline(targetGenomeFile, line)) {
    for(int i = 0; i < line.length(); i++) {
      char c = toupper(line[i]);
      if (c != 'N') {
        targetGenome += toupper(line[i]);
      }
    }
  }
  return targetGenome;
}

void SCCGC::run() {
  cout << "Running SCCGC" << endl;

  kmer_size = 21;

  // open files
  ifstream inputFile(this->inputFilePath);
  
  if(!inputFile.is_open()) {
    std::cout << "Error: Failed to open input file" << std::endl;
  }

  // tempfile
  string tempfilePath = outputDirPath + "/intermediate.txt";
  ofstream tempfile(tempfilePath);
  if(!tempfile.is_open()) {
    std::cout << "Error: Failed to open intermediate file:" << tempfilePath << std::endl;
  }

  // parse reference genome file
  referenceSeq = parseReferenceGenome(referenceGenomePath);
  cout << referenceSeq << endl;
}