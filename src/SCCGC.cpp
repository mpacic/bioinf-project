#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using HashTable = std::unordered_map<std::string, std::vector<int>>;

class SCCGC {
 public:
  SCCGC(std::string referenceGenomePath, std::string inputFilePath,
        std::string outputDirPath)
      : referenceGenomePath(referenceGenomePath),
        inputFilePath(inputFilePath),
        outputDirPath(outputDirPath){};
  ~SCCGC(){};
  void run();

 private:
  std::string referenceGenomePath;
  std::string inputFilePath;
  std::string outputDirPath;
  std::string referenceSeq;
  int kmer_size;
  static const int segment_length = 30000;
  static const int maxchar = 67108864;
  static const int ght_maxlen = 268435456;  // max size of the global hash table
  std::vector<int> kmer_location;           // global hash table
  std::vector<int> next_kmer;  // linked list of kmers with the same hashcode

  void buildGlobalHashTable(const string reference, int kmer_length);
  HashTable makeLocalHashTable(const string reference, int kmer_length);
  std::vector<std::pair<int, int>> getLowercasePositions(const string input);
  std::vector<std::pair<int, int>> getNPositions(const string input);

  void matchLocal(const string target, const string reference, int kmer_length);
  void matchGlobal(const string target, const string reference,
                   int kmer_length);
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
    targetGenome += line;
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

  // local matching phase
  cout << "Local matching phase... " << std::endl;
  matchLocal(targetSeq, referenceSeq, kmer_size);

  // Global matching phase
  cout << "Global matching phase... " << std::endl;
  matchGlobal(targetSeq, referenceSeq, kmer_size);

  // cout << "Building global hash table... " << std::endl;
  // buildGlobalHashTable(referenceSeq, 21);
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

// expects preprocessed reference segment
HashTable SCCGC::makeLocalHashTable(const string reference, int kmer_length) {
  HashTable kmer_location_map;
  int length = reference.length();

  for (int i = 0; i < length - kmer_length + 1; i++) {
    std::string current = reference.substr(i, kmer_length);
    if (kmer_location_map.count(current) > 0) {
      kmer_location_map[current].push_back(i);
    } else {
      kmer_location_map[current] = std::vector<int>(1, i);
    }
  }

  return kmer_location_map;
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

std::vector<std::pair<int, int>> SCCGC::getNPositions(const string input) {
  std::vector<std::pair<int, int>> positions;
  bool multiple = false;
  int start = 0;
  for (int i = 0; i < input.length(); i++) {
    if (input[i] == 'N') {
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

  // check if the last subsequence is all N characters
  if (multiple) {
    positions.push_back(std::make_pair(start, input.length() - 1));
  }
  return positions;
}

void SCCGC::matchLocal(const string reference, const string target, int kmer_size) {
  
}

 //global matching
void SCCGC::matchGlobal(const string reference, const string target, int kmer_size) {
  std::vector<std::pair<int, int>> nPositions = getNPositions(target);

  std::ofstream interimFile("tempFile.txt");
  for (const auto& pos : nPositions) {
    interimFile << pos.first << " " << pos.second << std::endl;
  }
  interimFile.close();

  std::string processedRef = reference;
  std::string processedTarg = target;

  for (const auto& pos : nPositions) {
    processedRef.erase(pos.first, pos.second - pos.first + 1);
    processedTarg.erase(pos.first, pos.second - pos.first + 1);
  }

  buildGlobalHashTable(processedRef, kmer_size);
}
