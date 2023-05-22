#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>
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
  void buildGlobalHashTable(string read, int kmer_length);

 private:
  std::string referenceGenomePath;
  std::string inputFilePath;
  std::string outputDirPath;
  std::string referenceSeq;
  int kmer_size;
  static const int maxchar = 67108864;
  static const int maxseq = 268435456; 
  std::vector<int> kmer_location;
  std::vector<int> next_kmer;
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
  return 0;
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

  buildGlobalHashTable(referenceSeq, 21);
}

void SCCGC::buildGlobalHashTable(string reference, int kmer_length) {
  
  int iters = reference.length() - kmer_length + 1;

  // allocate hash table
  kmer_location = std::vector<int>(maxseq);
  next_kmer = std::vector<int>(iters);
  
  for (int i = 0; i < maxseq; i++) {
    kmer_location[i] = -1;
  }

  for (int i = 0; i < iters; i++) {
    cout << i << endl;
    string kmer = reference.substr(i, kmer_length);
    hash<string> hasher;
    long key = labs(hasher(kmer));

    if (key == -2147483648) {
      key = 0;
    }

    while (key > maxseq - 1) {
      key = key / 2;
    }

    cout << key << endl;

    next_kmer[i] = kmer_location[static_cast<int>(key)];
    kmer_location[static_cast<int>(key)] = i;
  }
}
