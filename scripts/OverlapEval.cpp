#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <assert.h>

#include <stdlib.h>

// typedefs
typedef unsigned int uint;

// namespaces
using namespace std;
namespace CONST {
  char colon= ':';
  const uint MAXD = 50220423;
  const string groundTruth("/home/qw376/complex_ranker_scores/groundTruth");
}

// methods
void ReportQualityMetrics(const string groundTruth,
                          const string candidates,
                          const vector<size_t> &cutoffLevels);
vector<uint> ParseMetadataUint(string str); 
uint StrToUint(const string number_in_string);
void ComputeOverlap(const map<uint, vector<uint> > &groundTruthMap,
                    const map<uint, vector<uint> > &candidatesMap,
                    const vector<size_t> &cutoffLevels,
                    vector<vector<uint> > &overlaps);
void ReportOverlap(const vector<size_t> &cutoffLevels,
                   const vector<vector<uint> > &overlaps);
uint GetOverlapAtK(const vector<uint> &truth,
                   const vector<uint> &cand,
                   const size_t cutoff);

// main driver.
int main() {
  // cutoff levels.
  vector<size_t> cutoffLevels;
  //cutoffLevels.push_back(1);
  //cutoffLevels.push_back(2);
  //cutoffLevels.push_back(5);
  cutoffLevels.push_back(10);
  cutoffLevels.push_back(50);
  cutoffLevels.push_back(100);
  cutoffLevels.push_back(200);
  cutoffLevels.push_back(300);
  cutoffLevels.push_back(500);
  cutoffLevels.push_back(600);
  cutoffLevels.push_back(800);
  cutoffLevels.push_back(1000); 
  cutoffLevels.push_back(1200); 
  cutoffLevels.push_back(1500); 
  cutoffLevels.push_back(1800); 
  cutoffLevels.push_back(2000);

  // string candidatesFile = "/home/qw376/Candidates_Trec06/Candidates_T6";
  // string candidatesFile = "/home/qw376/Candidates_M9/NoSpaceConstraint_9000";
  // string candidatesFile = "/home/qw376/Candidates_M9/Lv16_Class5kall_1000_realcase_rewrite";
  // string candidatesFile = "/home/qw376/Candidates_M9/Lv16_Class5kall_1000_realcase_rewrite_rightrange"; //Lv16_Class5kall_1000_realcase_rewrite_3m_2012
  // string candidatesFile = "/home/qw376/Candidates_M9/M9_Size6_2k_0.5k";
  // string candidatesFile = "/home/qw376/WSDM16/Candidates_M9/Exp5_2k_500_S6";
  string candidatesFile = "/home/qw376/Data_for_SIGIR2016/Greedy_Candidates/Exp6_5k_3k";
  // string candidatesFile = "/home/qw376/Data_for_SIGIR2016/Prio_Candidates/Prio2K";
  ReportQualityMetrics(CONST::groundTruth,
                       candidatesFile,
                       cutoffLevels);

  return 0;
}

// Obtain the overlap @ given k and return the overlap @ cutoff.
// Assuming both input vector of same size.
uint GetOverlapAtK(const vector<uint> &truth,
                   const vector<uint> &cand,
                   const size_t cutoff) {
  assert(truth.size() >= cutoff);
  assert(cand.size() >= cutoff);
  uint num_overlaps = 0;
  for (size_t i = 0; i < cutoff; ++i) {  // foreach cand elem, check if it appears in truth.
    //for (size_t j = 0; j < cutoff; ++j) {
    for (size_t j = 0; j < 10; ++j) {
      if (cand[i] == truth[j]) {
//	cout<<cand[i]<<endl;
        ++num_overlaps;
        break;
      }
    }
  }
  assert(num_overlaps <= cutoff);
  return num_overlaps;
}

// Given the vector with cutoff levels and the overlaps at each level, report the overlap
// at various level across queries.
void ReportOverlap(const vector<size_t> &cutoffLevels,
                   const vector<vector<uint> > &overlaps) {
  // For each cutoff Level, compute the average overlap @ cutoff level and report it.
  for (size_t i = 0; i < cutoffLevels.size(); ++i) {
    uint num_overlaps = 0;
    for (size_t j = 0; j < overlaps.size(); ++j) { // foreach query
      num_overlaps += overlaps[j][i];
    }
    double avgOverlap = 0;
    if (overlaps.size() > 0) {
      avgOverlap = (double) num_overlaps / (double) overlaps.size();  // total overlaps / # queries
    }
    // Reporting.
    cout << "Average overlaps @ " << cutoffLevels[i] << " is " << setprecision(19) << avgOverlap << endl;
  }
}

void ReportQualityMetrics(const string groundTruth,
                          const string candidates,
                          const vector<size_t> &cutoffLevels) {
  // Map <qid -> vecOfTopResults> forming the ground truth (complex ranker top-k results).
  map<uint, vector<uint> > groundTruthMap;
  map<uint, vector<uint> >::iterator grIt;
 
  // Map <qid -> vecOfTopResults> forming the proposed method's results.
  map<uint, vector<uint> > candidatesMap;
  map<uint, vector<uint> >::iterator cIt;
 
  // Load the complex rankers top results.
  ifstream fH(groundTruth.c_str());
  if (fH.good()) {
    while (fH.good()) {
      string strQID, topResultsStr;
      getline(fH, strQID, CONST::colon);
      fH >> ws;
      getline(fH, topResultsStr);
      fH >> ws;
 
      // Tokenize string and convert it to uints.
      vector<uint> tmp = ParseMetadataUint(topResultsStr);
      uint qId = StrToUint(strQID);
 
      // Corner case.
      if (tmp.size() == 0) {
        cout << "Results from complex ranker for qid: " << qId << " = 0!! Exiting." << endl;
        exit(0);
      }
 
      // Check for duplicates
      grIt = groundTruthMap.find(qId);
      if (grIt != groundTruthMap.end()) { 
        cout << "QueryId: " << qId << " already exists in groundTruth map! Exit!" << endl;
        exit(0);
      }
 
      // Add to Map.
      groundTruthMap[ qId ] = tmp;
    }
  }
  fH.close();

  // Reporting.
  cout << "The top results under CR for "<< groundTruthMap.size()
       << " queries were loaded from file " << groundTruth << endl;

  // Load the proposed candidates.
  ifstream c_fH(candidates.c_str());
  if (c_fH.good()) {
    while (c_fH.good()) {
      string strQID, topResultsStr;
      getline(c_fH, strQID, CONST::colon);
      c_fH >> ws;
      getline(c_fH, topResultsStr);
      c_fH >> ws;
 
      // Tokenize string and convert it to uint.
      vector<uint> tmp = ParseMetadataUint(topResultsStr);
      uint qId = StrToUint(strQID);
 
      // Corner case.
      if (tmp.size() == 0) {
        cout << "Results from candidate method for qid: " << qId << " = 0!! Exiting." << endl;
        exit(0);
      }
 
      // Check for duplicates
      cIt = candidatesMap.find(qId);
      if (cIt != candidatesMap.end()) {
        cout << "QueryId: " << qId << " already exists in candidates map! Exit!" << endl;
        exit(0);
      }
 
      // Add to candidates map.
      // Note that there is no resizing in the candidate pool now, since this parameter is
      // provided when we conduct the experiments.
      candidatesMap[ qId ] = tmp;
    }
  }
  c_fH.close();
 
  // Reporting.
  cout << "The top results under our candidate approach for " << candidatesMap.size()
       << " queries were loaded from file " << candidates << endl;

  // Compute overlaps at various cutoffs and update the vector of vectors which contains
  // the overlap per query for various cutoffs.
  vector<vector<uint> > overlaps;
  ComputeOverlap(groundTruthMap, candidatesMap, cutoffLevels, overlaps);
  ReportOverlap(cutoffLevels, overlaps);
}

// Given a string, tokenize it and return the values as a vector of uints.
vector<uint> ParseMetadataUint(string str) {
  istringstream iss(str);
  vector<uint> metadata;
  do {
    string token;
    iss >> token;
    if (token=="") continue;
      metadata.push_back( StrToUint(token) );
    } while (iss);
  return metadata;
}

uint StrToUint(const string number_in_string) {
  uint number = 0;
  istringstream(number_in_string) >> number;
  return number;
}

// Compute overlap @ various cutoffs, given maps maintaining the top results for each method
// ground truth and the proposed method. The first vector is the ground truth and the second
// is the proposed one. The result is pushed in the given vector of vector.
void ComputeOverlap(const map<uint, vector<uint> > &groundTruthMap,
                    const map<uint, vector<uint> > &candidatesMap,
                    const vector<size_t> &cutoffLevels,
                    vector<vector<uint> > &overlaps) {
  // Largest cutoff level.
  const size_t maxCutoffLevel = cutoffLevels[ cutoffLevels.size() - 1 ];

  // For each query in candidateMap.
  map<uint, vector<uint> >::const_iterator it;
  map<uint, vector<uint> >::const_iterator grIt;
  for (it = candidatesMap.begin(); it != candidatesMap.end(); ++it) {
    // Find if the corresponding map exists.
    grIt = groundTruthMap.find(it->first);
    if (grIt == groundTruthMap.end()) {
      cout << "Queryid: " << it->first << " does not exist in the groundTruthMap. Exit!" << endl;
      exit(0);
    }
 
    // Copy to new vectors so we can resize.
    vector<uint> truth = grIt->second;;
    //assert(truth.size() >= maxCutoffLevel);
    truth.resize(maxCutoffLevel);
    vector<uint> cand = it->second;
 
    // Corner case where we return less than the ones in ground truth, resize to same with MAXD.
    if (cand.size() < maxCutoffLevel) {
      cand.resize(maxCutoffLevel, CONST::MAXD + 1);  // default value is MAXD + 1
    } else {
      cand.resize(maxCutoffLevel);
    }
    assert(cand.size() >= maxCutoffLevel);
    assert(truth.size() == cand.size());
 
    // Get Overlap @ various cutoffs.
    vector<uint> queryOverlaps;
    for (size_t i = 0; i < cutoffLevels.size(); ++i) {
    //  cout << "Computing overlap @ " << cutoffLevels[i] << " truthSz: "
    //       << truth.size() << " candSz:" << cand.size() << endl;
       queryOverlaps.push_back( GetOverlapAtK(truth, cand, cutoffLevels[i]) );
    }

    // Sanity checks.
    assert(queryOverlaps.size() == cutoffLevels.size());
    assert(queryOverlaps.size() > 1);
    for (size_t i = 1; i < queryOverlaps.size() - 1; ++i) {
      assert(queryOverlaps[i] <= queryOverlaps[i+1]);
    }
    overlaps.push_back(queryOverlaps);  // finally push back.
    
    cout<<"qid: "<<it->first<<endl; 
    // Report per query.
    for (size_t i = 0; i < cutoffLevels.size(); ++i) {
      cout << "Per query overlap @ " << cutoffLevels[i] << " " << queryOverlaps[i] << endl;
    }
  }
}
