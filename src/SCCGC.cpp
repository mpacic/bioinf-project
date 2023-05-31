#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

class SCCGC {
 public:
  SCCGC(std::string referenceGenomePath, std::string inputFilePath,
        std::string outputDirPath)
      : referenceGenomePath(referenceGenomePath),
        inputFilePath(inputFilePath),
        outputDirPath(outputDirPath){};
  ~SCCGC(){};
  void run();
  void buildGlobalHashTable(const string reference, int kmer_length);
  void buildLocalHashTable(const string reference, int kmer_length);
  std::vector<std::pair<int, int>> getLowercasePositions(const string input);

 private:
  std::string referenceGenomePath;
  std::string inputFilePath;
  std::string outputDirPath;
  std::string referenceSeq;
  int kmer_size;
  static const int maxchar = 67108864;
  static const int ght_maxlen = 268435456;
  std::vector<int> kmer_location;
  std::vector<int> next_kmer;
};

int main(int argc, char **argv) {
  // check number of arguments
  if (argc < 4) {
    std::cout << "Usage: " << argv[0]
              << " <reference genome file> <input file> <output_directory>"
              << std::endl;
    return 1;
  }

  // check reference genome file exists
  if (!filesystem::exists(argv[1])) {
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

  SCCGC sccgc(argv[1], argv[2], argv[3]);

  sccgc.run();
  return 0;
}

std::string parseReferenceGenome(std::string referenceGenomePath) {
  ifstream referenceGenomeFile(referenceGenomePath);
  // check file opened successfully
  if (!referenceGenomeFile.is_open()) {
    std::cout << "Error: Failed to open reference genome file" << std::endl;
  }

  string referenceGenome = "";
  string line;
  // skip first line
  getline(referenceGenomeFile, line);

  while (getline(referenceGenomeFile, line)) {
    for (int i = 0; i < line.length(); i++) {
      char c = toupper(line[i]);
      if (c != 'N') {
        referenceGenome += c;
      }
    }
  }
  return referenceGenome;
}

std::string readTargetGenome(std::string targetGenomePath) {
  ifstream targetGenomeFile(targetGenomePath);
  // check file opened successfully
  if (!targetGenomeFile.is_open()) {
    std::cout << "Error: Failed to open target genome file" << std::endl;
  }

  string targetGenome = "";
  string line;
  // skip first line
  getline(targetGenomeFile, line);

  while (getline(targetGenomeFile, line)) {
    for (int i = 0; i < line.length(); i++) {
      targetGenome += line[i];
    }
  }
  return targetGenome;
}

void SCCGC::run() {
  cout << "Running SCCGC" << endl;

  kmer_size = 21;

  // open files
  ifstream inputFile(this->inputFilePath);

  if (!inputFile.is_open()) {
    std::cout << "Error: Failed to open input file" << std::endl;
  }

  // tempfile
  string tempfilePath = outputDirPath + "/intermediate.txt";
  ofstream tempfile(tempfilePath);
  if (!tempfile.is_open()) {
    std::cout << "Error: Failed to open intermediate file:" << tempfilePath
              << std::endl;
  }

  // parse reference genome file
  cout << "Parsing reference sequence... " << std::endl;
  referenceSeq = parseReferenceGenome(referenceGenomePath);

  // read target genome file
  cout << "Reading target sequence... " << std::endl;
  string targetSeq = readTargetGenome(inputFilePath);

  cout << "Building global hash table... " << std::endl;
  buildGlobalHashTable(referenceSeq, 21);
}

void SCCGC::buildGlobalHashTable(const string reference, int kmer_length) {
  // there are L - k + 1 kmers in the sequence
  int iters = reference.length() - kmer_length + 1;

  // allocate hash table
  kmer_location = std::vector<int>(ght_maxlen);
  next_kmer = std::vector<int>(iters);

  // set all hash table entries to default value
  std::fill(kmer_location.begin(), kmer_location.end(), -1);

  // calculate hashcode for every kmer
  for (int i = 0; i < iters; i++) {
    string kmer = reference.substr(i, kmer_length);
    hash<string> hasher;
    long key = labs(hasher(kmer)) % ght_maxlen;

    next_kmer[i] = kmer_location[static_cast<int>(key)];
    kmer_location[static_cast<int>(key)] = i;
  }
}

void SCCGC::buildLocalHashTable(const string reference, int kmer_length) {
  // build local hash table
}

std::vector<std::pair<int, int>> SCCGC::getLowercasePositions(
    const string input) {
  std::vector<std::pair<int, int>> positions;
  bool multiple = false;
  int start = 0;
  for (int i = 0; i < input.length(); i++) {
    if (std::islower(input[i])) {
      if (multiple) {
        continue;
      } else {
        start = i;
        multiple = true;
      }
    } else {
      if (multiple) {
        positions.push_back(std::make_pair(start, i));
        multiple = false;
      }
      start = 0;
    }
  }

  // check if the last subsequence is all lowercase
  if (multiple) {
    positions.push_back(std::make_pair(start, input.length() - 1));
  }
  return positions;
}