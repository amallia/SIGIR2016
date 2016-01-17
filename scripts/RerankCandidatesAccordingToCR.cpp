#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

namespace CONST{
	char Colon = ':';
	// const string CandidatesDir("/home/qw376/Candidates_SpecialQuery/Lookup0.5k_2012_0.5index_realcase_1000");
	const string CandidatesDir("/home/qw376/WSDM16/Candidates_M9/Exp4_2k_0.5k");
	const string CRResultsDir("/home/qw376/HumanJudgedQueries/ComplexRankerResult_10_11_12");
	const string GAResultsDir("/home/qw376/WSDM16/Candidates_M9/Exp4_2k_0.5k_Rerank");
}

typedef unsigned int uint;

vector<uint> ParseMetadataUint(string str);
uint StrToUint (const string number_in_string);
void LoadCandidates(map <uint, vector<uint> > &CandidatesMap);
void LoadCRMap(map <uint, vector<uint> > &CRMap);
void OutputGAMap(map <uint, vector<uint> > &CandidatesMap, map <uint, vector<uint> > &CRMap, map <uint, vector<uint> > &GAMap);


int main(){

	map<uint, vector<uint> > CandidatesMap;
	map<uint, vector<uint> > CRMap;
	map<uint, vector<uint> > GAMap;

	LoadCandidates(CandidatesMap);
	LoadCRMap(CRMap);
	OutputGAMap(CandidatesMap, CRMap, GAMap);

	return 0;
}


void LoadCandidates(map <uint, vector<uint> > &CandidatesMap){
	ifstream inputf(CONST::CandidatesDir.c_str());
	if(inputf.good()){
		while(inputf.good()){
			string StrQid, StrCandidates;
			getline(inputf, StrQid, CONST::Colon);
			inputf >> ws;
			getline(inputf, StrCandidates);
			inputf >> ws;
			uint qid = StrToUint(StrQid);
			vector<uint> tmp = ParseMetadataUint(StrCandidates);
			CandidatesMap[qid] = tmp;
		}
	}
	inputf.close();
}

void LoadCRMap(map <uint, vector<uint> > &CRMap){
	ifstream inputf(CONST::CRResultsDir.c_str());
	if(inputf.good()){
		while(inputf.good()){
			string StrQid, StrCR;
			getline(inputf, StrQid, CONST::Colon);
			inputf >> ws;
			getline(inputf, StrCR);
			inputf >> ws;
			uint qid = StrToUint(StrQid);
			vector<uint> tmp = ParseMetadataUint(StrCR);
			CRMap[qid] = tmp;
		}
	}
	inputf.close();
}

void OutputGAMap(map <uint, vector<uint> > &CandidatesMap, map <uint, vector<uint> > &CRMap, map <uint, vector<uint> > &GAMap){

	map <uint, vector<uint> >::iterator it;
	for(it = CRMap.begin(); it!=CRMap.end(); ++it){
		vector <uint> tmp;
		for(int i = 0; i < it->second.size(); ++i){
			for(int j = 0; j < CandidatesMap[it->first].size(); ++j){
				if( it->second[i] == CandidatesMap[it->first][j] ){
					tmp.push_back(it->second[i]);
					break;
				}
			}
		}
		GAMap[it->first] = tmp;
		tmp.clear();
	}

	// for(it = GAMap.begin(); it!=GAMap.end(); ++it){	
	// 	cout << it->first << ":" << it->second.size() << endl;
	// 	cout << it->first << ":";
	// 	for(int i = 0; i < it->second.size(); ++i){
	// 		cout<<it->second[i]<<" ";
	// 	}
	// 	cout<<endl;
	// }

	ofstream outputf(CONST::GAResultsDir.c_str());
	if(outputf.good()) {
		for(it = GAMap.begin(); it!=GAMap.end(); ++it){
			outputf << it->first << ":";
			for(int i = 0; i < it->second.size(); ++i){
				outputf<<it->second[i]<<" ";
			}
			outputf<<endl;
		}
	}
	outputf.close();
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