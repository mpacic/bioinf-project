#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <unistd.h>

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
  std::string interimFilePath;
  std::string outputDirPath;
  std::string referenceSeq;
  int kmer_size;
  static const int segment_length = 30000;
  static const int maxchar = 67108864;
  static const int ght_maxlen = 268435456;  // max size of the global hash table
  std::vector<int> kmer_location;           // global hash table
  std::vector<int> next_kmer;  // linked list of kmers with the same hashcode
  float T1 = 0.5; // threshold for local matching
  int T2 = 4; // similarity threshold
  bool global = false;
  std::string targetHeader;
  int lineLength;

  void buildGlobalHashTable(const string reference, int kmer_length);
  HashTable makeLocalHashTable(const string reference, int kmer_length);
  std::vector<std::pair<int, int>> getLowercasePositions(const string input);
  std::vector<std::pair<int, int>> getNPositions(const string input);

  void matchLocal(const string target, const string reference, int kmer_length);
  void matchGlobal(const string target, const string reference,
                   int kmer_length);

  void postprocess();
  void run7zip(const string filename);

  bool allN(const string input);
};

unsigned long long getMemoryUsageInKB() {
    std::ifstream statm("/proc/self/statm");
    unsigned long long size, resident, share, text, lib, data, dt;
    statm >> size >> resident >> share >> text >> lib >> data >> dt;
    return resident * (unsigned long long)getpagesize() / 1024;
}

void printMemoryUsage() {
  long long memusage = getMemoryUsageInKB();
  cout << "Memory usage: " << memusage << " KB" << endl;
}

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
    std::exit(1);
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

std::string readReferenceGenome(std::string referenceGenomePath) {
  ifstream referenceGenomeFile(referenceGenomePath);
  // check file opened successfully
  if (!referenceGenomeFile.is_open()) {
    std::cout << "Error: Failed to open reference genome file" << std::endl;
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

std::string readTargetGenome(std::string targetGenomePath) {
  ifstream targetGenomeFile(targetGenomePath);
  // check file opened successfully
  if (!targetGenomeFile.is_open()) {
    std::cout << "Error: Failed to open target genome file" << std::endl;
    std::exit(1);
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
  interimFilePath = outputDirPath + "/interim.txt";
  ofstream outputStream(outputDirPath + "/output.sccg");

  if (!inputFile.is_open()) {
    std::cout << "Error: Failed to open input file" << std::endl;
    std::exit(1);
  }

  if (!outputStream.is_open()) {
    std::cout << "Error: Failed to open output file" << std::endl;
    std::exit(1);
  }

  std::getline(inputFile, targetHeader);

  std::string line;
  std::getline(inputFile, line);
  int lineLength = line.length();
  inputFile.close();

  // parse reference genome file
  cout << "Parsing reference sequence... " << std::endl;
  referenceSeq = readReferenceGenome(referenceGenomePath);
  std::transform(referenceSeq.begin(), referenceSeq.end(), referenceSeq.begin(), ::toupper);

  // read target genome file
  cout << "Reading target sequence... " << std::endl;
  string targetSeq = readTargetGenome(inputFilePath);

  //write target header and line length to output file
  outputStream << targetHeader << std::endl << lineLength << std::endl;

  std::vector<std::pair<int, int>> lowercasePositions = getLowercasePositions(targetSeq);
  // convert target sequence to uppercase
  std::transform(targetSeq.begin(), targetSeq.end(), targetSeq.begin(), ::toupper);

  // write lowercase positions to file
  int temp_end = 0;
  for (const auto& pos : lowercasePositions) {
    outputStream << pos.first - temp_end << " " << pos.second - pos.first << " ";
    temp_end = pos.second;
  }
  outputStream << std::endl << std::endl;

  // local matching phase
  cout << "Local matching phase... " << std::endl;
  matchLocal(targetSeq, referenceSeq, kmer_size);

  if (!global) {
  // write local matching result to output file
    cout << "Postprocessing... " << std::endl;
    postprocess();

    // delete interim file
    // std::remove(interimFilePath.c_str());

    // close output file
    outputStream.close();

    // run PPMd using 7zip
    cout << "Running PPMd... " << std::endl;
    std::string outputFilePath = outputDirPath + "/output.sccg";
    run7zip(outputFilePath);

    // delete uncompressed output file
    // std::remove(outputFilePath.c_str());

    printMemoryUsage();

    std::exit(0);
  }

  // Global matching phase

  cout << "Global matching phase... " << std::endl;

  // write N positions to file
  std::vector<std::pair<int, int>> nPositions = getNPositions(targetSeq);
  for (const auto& pos : nPositions) {
    outputStream << pos.first << " " << pos.second << " ";
  }
  outputStream << std::endl;
  
  matchGlobal(targetSeq, referenceSeq, kmer_size);

  cout << "Postprocessing... " << std::endl;
  postprocess();

  // delete interim file
  // std::remove(interimFilePath.c_str());

  // close output file
  outputStream.close();

  // run PPMd using 7zip
  cout << "Running PPMd... " << std::endl;
  std::string outputFilePath = outputDirPath + "/output.sccg";
  run7zip(outputFilePath);

  printMemoryUsage();

  // delete uncompressed output file
  // std::remove(outputFilePath.c_str());
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
HashTable SCCGC::makeLocalHashTable(
    const string reference, int kmer_length) {
  HashTable kmer_location_map(reference.length() - kmer_length + 1, std::hash<std::string>{});
  int length = reference.length();
  kmer_location_map.reserve(length - kmer_length + 1);

  for (int i = 0; i < length - kmer_length + 1; i++) {
    std::string current = reference.substr(i, kmer_length);
    kmer_location_map[current].push_back(i);
  }

  return kmer_location_map;
}

// returns a vector of pairs of start and end positions of lowercase subsequences
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

// local matching
void SCCGC::matchLocal(const string target, const string reference,
                      int kmer_length) {
  std::ofstream interimStream(interimFilePath);
  long total_length = std::min(target.length(), reference.length());
  long num_segments = total_length / segment_length;
  int unmatched_segments = 0;

  if (num_segments < 5) {
    global = true;
    return;
  }

  // debug info
  cout << "total_length:" << total_length << endl;
  cout << "target.length():" << target.length() << endl;
  cout << "reference.length():" << reference.length() << endl;
  cout << "num_segments:" << num_segments << endl;

  for (int i = 0; i < num_segments; i++) {
    string t_seg = target.substr(i * segment_length, segment_length);
    string r_seg = reference.substr(i * segment_length, segment_length);
    HashTable hashtable = makeLocalHashTable(r_seg, kmer_length);
    std::ostringstream ss;
    for (int j = 0; j < segment_length; j++) {
      if (segment_length - j < kmer_length) {
        ss << t_seg[j];
        continue;
      }
      string kmer = t_seg.substr(j, kmer_length);
      if (hashtable.count(kmer) > 0) { // check if key exists
        vector<int> kmer_positions = hashtable[kmer];
        int longest_len = -1;
        int longest_pos = -1;
        for (int pos : kmer_positions) {
          if (t_seg.substr(j, kmer_length) == r_seg.substr(pos, kmer_length)) {
            int ext = 0; // match extension length
            // find longest match between target and reference 
            while (j + kmer_length + ext < t_seg.length()
                && pos + kmer_length + ext < r_seg.length()
                && r_seg[pos + kmer_length + ext] == t_seg[j + kmer_length + ext]) {
              ext++;
            }
            // if current match is longer than previous longest match, update
            if (kmer_length + ext - 1 > longest_len) {
              longest_len = kmer_length + ext - 1;
              longest_pos = pos;
            }
          }
        }

        if (longest_len == -1) {
          ss << t_seg[j];
          continue;
        }
        
        // write to file
        if (j != 0) {
          ss << endl;
        }
        int start = i*segment_length + longest_pos;
        int end = start + longest_len;
        ss << start << "," << end << endl;
        // write the unmatched character at the end
        if (longest_pos+longest_len + 1 < segment_length) {
          ss << t_seg[longest_pos+longest_len+1];
        }

        // update index to skip over longest match and unmatched character
        j += longest_len + 1;
      } else {
        // write unmatched character to file
        ss << t_seg[j];
      }
    }

    // check ratio of directly stored characters
    if (ss.str().length() > segment_length * T1 && !allN(t_seg)) {
      unmatched_segments++;
    }

    // check number of unmatched segments
    if (unmatched_segments > T2) {
      global = true;
      interimStream.close();
      std::remove(interimFilePath.c_str());
      cout << "exited from local matching on segment " << i << endl;
      return;
    }

    // write newline to file after segment
    // interimStream.write(ss.str().c_str(), ss.str().length());
    interimStream << ss.str();
    interimStream.flush();
  }

  // last segment
  string t_seg = target.substr(num_segments * segment_length, target.length() - num_segments * segment_length);
  string r_seg = reference.substr(num_segments * segment_length, reference.length() - num_segments * segment_length);
  HashTable hashtable = makeLocalHashTable(r_seg, kmer_length);
  for (int j = 0; j < t_seg.length() - kmer_length + 1; j++) {
    string kmer = t_seg.substr(j, kmer_length);
    if (hashtable.count(kmer) > 0) { // check if key exists
      vector<int> kmer_positions = hashtable[kmer];
      int longest_len = 0;
      int longest_pos = 0;
      for (int pos : kmer_positions) {
        int ext = 0; // match extension length
        // find longest match between target and reference 
        while (j + kmer_length + ext < t_seg.length()
            && pos + kmer_length + ext < r_seg.length()
            && r_seg[pos + kmer_length + ext] ==  t_seg[j + kmer_length + ext]) {
          ext++;
        }
        // if current match is longer than previous longest match, update
        if (kmer_length + ext - 1 > longest_len) {
          longest_len = kmer_length + ext - 1;
          longest_pos = pos;
        }
      }

      // write to file
      if (j != 0) {
        interimStream << endl;
      }
      int start = num_segments*segment_length + longest_pos;
      int end = start + longest_len;
      interimStream <<  start << "," << end << endl;

      // update index to skip over longest match
      j += longest_len - 1;
    } else {
      interimStream << t_seg[j];
    }
  }
  interimStream << endl;
  interimStream.close();
}

void SCCGC::run7zip(const string filename) {
  string command = "../7za a -m0=PPMd " + filename + ".7z " + filename + " > /dev/null";
  system(command.c_str());
}

// merging of continuous matches and delta encoding
void SCCGC::postprocess() {
  
  std::ifstream interimStream(interimFilePath);
  std::string line;
  std::vector<int> positions;
  std::stringstream ss;

  // merge continuous matches
  while (std::getline(interimStream, line)) {
    if (line.find(',') != std::string::npos) {
      int start = stoi(line.substr(0, line.find(',')));
      int end = stoi(line.substr(line.find(',') + 1));
      if (positions.size() > 0) {
        if (start != positions.back() + 1) {
          ss << start << "," << end << std::endl;
          positions.clear();
        }
      }
      positions.push_back(start);
      positions.push_back(end);
    } else if (line.length() > 0) {
      if (positions.size() > 0) {
        ss << positions.front() << "," << positions.back() << std::endl;
        positions.clear();
      }
      ss << line << std::endl;
    }
  }
  if (positions.size() > 0) {
    ss << positions.front() << "," << positions.back() << std::endl;
  }
  interimStream.close();

  // delta encoding
  std::stringstream oss;
  bool successive = false;
  int prev = 0;
  while (std::getline(ss, line)) {
    if (line.find(',') != std::string::npos) {
      int start = stoi(line.substr(0, line.find(',')));
      int end = stoi(line.substr(line.find(',') + 1));
      oss << start - prev << "," << end - start << std::endl;
      prev = end;
    } else if (line.length() > 0) {
      oss << line << std::endl;
    }
  }

  // write to file
  std::ofstream outputStream(outputDirPath + "/output.sccg", std::ios::app);
  outputStream << oss.str() << endl;
  outputStream.close();
}

bool SCCGC::allN(const string input) {
  for (int i = 0; i < input.length(); i++) {
    if (input[i] != 'N') {
      return false;
    }
  }
  return true;
}
