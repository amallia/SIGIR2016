#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

using namespace std;

namespace CONST{
	char Colon = ':';
	char Space = ' ';
	string Start = "51";
	// const string CandidatesDir("/home/qw376/Candidates_SpecialQuery/Lookup3k_5543_0.5index_realcase_1000_CRRerank"); //ComplexRankerResult_10_11_12
	const string CandidatesDir("/home/qw376/WSDM16/Candidates_M9/Exp4_2k_0.5k_Rerank"); //ComplexRankerResult_10_11_12
	// const string CandidatesDir("/home/qw376/Candidates_SpecialQuery/ComplexRankerResult_10_11_12");
	const string CwidToTrecidMappingFile("/home/qw376/Info_Clueweb/basicDocFeatures");
	const string QidToTrecQidMappingFile("/home/qw376/HumanJudgedQueries/QidtoTrecQid_10_11_12_Takingout2");
	const string RelevantInfoFile("/home/qw376/HumanJudgedQueries/10_11_12.qrels");
	int k = 1; //1, 5, 10, 20, 100, 200
}

struct ReleInfo{
	uint did;
	int judge;
};

bool mysortfuc(ReleInfo a, ReleInfo b){ return (a.judge > b.judge); }

typedef unsigned int uint;

vector<uint> ParseMetadataUint(string str);
uint StrToUint(const string number_in_string);
int StrToint(const string number_in_string);
void LoadCandidates(map<uint, vector<uint> > &CandidatesMap);
void LoadRelevantInfo(map<uint, uint> &QidToTrecQidMap, map<string, uint> &CwidToTrecidMap, map<uint, map<uint, int> > &RelevantInfoMap);
void LoadCwidToTrecidMapping(map<string, uint> &CwidToTrecidMap);
void LoadQidToTrecQidMapping(map<uint, uint> &QidToTrecQidMap);
void ComputeDCG(map<uint, vector<uint> > &CandidatesMap, map<uint, map<uint, int> > &RelevantInfoMap, map<uint, float> &DCGMap);
void ComputeIDCG(map<uint, map<uint, int> > &RelevantInfoMap, map<uint, float> &IDCGMap);
void ComputeNDCG(map<uint, float> &DCGMap, map<uint, float> &IDCGMap);

int main(){

	map<uint, vector<uint> > CandidatesMap;
	map<string, uint> CwidToTrecidMap;
	map<uint, map<uint, int> > RelevantInfoMap;
	map<uint, float> DCGMap;
	map<uint, float> IDCGMap;
	map<uint, uint> QidToTrecQidMap;

	LoadQidToTrecQidMapping(QidToTrecQidMap);
	LoadCwidToTrecidMapping(CwidToTrecidMap);
	LoadCandidates(CandidatesMap);
	LoadRelevantInfo(QidToTrecQidMap, CwidToTrecidMap, RelevantInfoMap);
	ComputeDCG(CandidatesMap, RelevantInfoMap, DCGMap);
	ComputeIDCG(RelevantInfoMap, IDCGMap);

	map<uint, float>::iterator it;
	for(it = DCGMap.begin(); it!=DCGMap.end(); ++it){
		cout<< it->first << " " <<it->second << endl;
	}

	for(it = IDCGMap.begin(); it!=IDCGMap.end(); ++it){
		cout<< it->first << " " <<it->second << endl;
	}

	ComputeNDCG(DCGMap, IDCGMap);

	return 0;
}

void LoadQidToTrecQidMapping(map<uint, uint> &QidToTrecQidMap){
	ifstream inputf(CONST::QidToTrecQidMappingFile.c_str());
	if(inputf.good()){
		while(inputf.good()){
			string Qids, TrecQids;
			getline(inputf, Qids, CONST::Space);
			inputf >> ws;
			getline(inputf, TrecQids);
			inputf >> ws;
			uint qid = StrToUint(Qids);
			uint TrecQid = StrToUint(TrecQids);
			QidToTrecQidMap[qid] = TrecQid;
		}
	}

}

void LoadCwidToTrecidMapping(map<string, uint> &CwidToTrecidMap){
	ifstream inputf(CONST::CwidToTrecidMappingFile.c_str());
	if(inputf.good()){
		while(inputf.good()){
			string cwids, dids, redund;
			getline(inputf, dids, CONST::Space);
			inputf >> ws;
			getline(inputf, cwids, CONST::Space);
			inputf >> ws;
			getline(inputf, redund);
			inputf >> ws;
			uint did = StrToUint(dids);
			CwidToTrecidMap[cwids] = did;
		}
	}
	inputf.close();

}

void LoadRelevantInfo(map<uint, uint> &QidToTrecQidMap, map<string, uint> &CwidToTrecidMap, map<uint, map<uint, int> > &RelevantInfoMap){
	string qidrs = CONST::Start;
	map<uint, int> CwidToJugMap;
	map<string, uint>::iterator it;
	map<uint, uint>::iterator it1;

	ifstream inputf(CONST::RelevantInfoFile.c_str());
	if(inputf.good()){
		while(inputf.good()){
			string qids, cwid, judge, redund;
			uint qid, qidr, trecqid;
			int judgeint;
			getline(inputf, qids, CONST::Space);
			if(qids.empty()){
				qid = StrToUint(qidrs);
				it1  = QidToTrecQidMap.find(qid);
				if(it1 != QidToTrecQidMap.end()){
					trecqid = it1->second; 
				}
				RelevantInfoMap[trecqid] = CwidToJugMap;
				CwidToJugMap.clear();
				break;
			}

			qid = StrToUint(qids);
			qidr = StrToUint(qidrs);

			if(qid > qidr){
				qid = StrToUint(qidrs);
				it1  = QidToTrecQidMap.find(qid);
				if(it1 != QidToTrecQidMap.end()){
					trecqid = it1->second; 
				}
				RelevantInfoMap[trecqid] = CwidToJugMap; 
				CwidToJugMap.clear();
			}
			qidrs = qids;
			getline(inputf, redund, CONST::Space);
			getline(inputf, cwid, CONST::Space);
			getline(inputf, judge);

			judgeint = StrToint(judge);
			// cout<<cwid<<" "<<judgeint<<endl;
			it = CwidToTrecidMap.find(cwid);
			if(it!=CwidToTrecidMap.end()){
				CwidToJugMap[it->second] = judgeint; 
			}

		}
	}
	inputf.close();
}

void ComputeDCG(map<uint, vector<uint> > &CandidatesMap, map<uint, map<uint, int> > &RelevantInfoMap, map<uint, float> &DCGMap){
	
	// map<uint, map<uint, int> > ::iterator it;
	// map<uint, int>:: iterator it1;
	// for(it = RelevantInfoMap.begin(); it!=RelevantInfoMap.end(); it++){
	// 	cout<<it->first<<": ";
	// 	for(it1 = it->second.begin(); it1!=it->second.end(); it1++){
	// 		cout<<it1->first<<" <-> "<<it1->second<<" ";
	// 	}
	// 	cout<<endl;
	// }

 // 	map<uint, vector<uint> >::iterator it;
	// for(it = CandidatesMap.begin(); it!=CandidatesMap.end(); ++it){
	// 	// cout<<it->first<<": "<<RelevantInfoMap[it->first].size()<<endl;
	// 	cout<<it->first<<": ";
	// 	for(int i = 0; i < 10; ++i){
	// 		cout<< it->second[i] << "<->" << RelevantInfoMap[it->first][it->second[i]] << " ";
	// 	}	
	// 	cout<<endl;
	// }


	map<uint, vector<uint> >::iterator it;
	for(it = CandidatesMap.begin(); it!=CandidatesMap.end(); ++it){
		float DCG = 0;
		cout<<it->first<<": ";
		for(int i = 0; i < CONST::k; ++i){
			float JudInfo = RelevantInfoMap[it->first][it->second[i]];
			if(JudInfo < 0) JudInfo = 0;
			// if(JudInfo = 4) JudInfo = 3;
			if(i == 0){
				DCG = DCG + JudInfo;
				cout<< it->second[i] << "," << JudInfo << "," << DCG << " ";
			}else{
				DCG = DCG + JudInfo/(log2(i+1));
				cout<< it->second[i] << "," << JudInfo/(log2(i+1)) << "," << DCG << " ";
			}
		}	
		cout<<endl;
		DCGMap[it->first] = DCG;
	}

}


void ComputeIDCG(map<uint, map<uint, int> > &RelevantInfoMap, map<uint, float> &IDCGMap){
	map<uint, vector<ReleInfo> > ReleInfoMap; 
	map<uint, map<uint, int> >::iterator it;
	map<uint, int>::iterator it1;
	for(it = RelevantInfoMap.begin(); it!=RelevantInfoMap.end(); ++it){
		vector<ReleInfo> TmpVector;
		for(it1 = it->second.begin(); it1!=it->second.end(); ++it1){
			ReleInfo Tmp;
			Tmp.did = it1->first;
			Tmp.judge = it1->second;
			TmpVector.push_back(Tmp);
		}

		sort(TmpVector.begin(), TmpVector.end(), mysortfuc);
		ReleInfoMap[it->first] = TmpVector;
	}

	// map<uint, vector<ReleInfo> >::iterator it2;
	// for(it2 = ReleInfoMap.begin(); it2 != ReleInfoMap.end(); ++it2){
	// 	cout<<it2->first<<" ";
	// 	for(int i = 0; i<it2->second.size(); ++i){
	// 		cout<<it2->second[i].did<<","<<it2->second[i].judge<<" ";
	// 	}
	// 	cout<<endl;
	// }
	map<uint, vector<ReleInfo> >::iterator it2;
	for(it2 = ReleInfoMap.begin(); it2 != ReleInfoMap.end(); ++it2){
		float IDCG = 0;
		cout<<it2->first<<": ";
		for(int i = 0; i< CONST::k; ++i){
			float JudInfo = it2->second[i].judge;
			if(JudInfo < 0) JudInfo = 0;
			// if(JudInfo = 4) JudInfo = 3;
			if(i == 0 ){
				IDCG = IDCG + JudInfo; 
				cout << it2->second[i].did << "," << JudInfo << "," << IDCG << " ";
			}else{
				IDCG = IDCG + JudInfo/(log2(i+1));
				cout << it2->second[i].did << "," << JudInfo/(log2(i+1)) << "," << IDCG << " ";
			}
		}
		cout<<endl;
		IDCGMap[it2->first] = IDCG;
	}

}

void ComputeNDCG(map<uint, float> &DCGMap, map<uint, float> &IDCGMap){
	map <uint, float>::iterator it;
	float sum = 0;
	for(it = IDCGMap.begin(); it!=IDCGMap.end(); ++it){
		if(it->second != 0)
			sum = sum + DCGMap[it->first]/it->second;
		else
			sum = sum + 1;
	}
	float NDCG = sum / DCGMap.size();
	cout<< NDCG << " " << sum << " " << DCGMap.size() << endl; 
}


void LoadCandidates(map<uint, vector<uint> > &CandidatesMap){
	ifstream inputf1(CONST::CandidatesDir.c_str());
	if(inputf1.good()){
		while(inputf1.good()){
			string strQID, topResultsStr;
			getline(inputf1, strQID, CONST::Colon);
			inputf1 >> ws;
			getline(inputf1, topResultsStr);
			inputf1 >> ws;
			vector<uint> tmp = ParseMetadataUint(topResultsStr);
			uint Qid = StrToUint(strQID);
			CandidatesMap[Qid] = tmp; 
		}
	}	
	inputf1.close();

	// map<uint, vector<uint> >::iterator it;
	// for(it = CandidatesMap.begin(); it!=CandidatesMap.end(); it++){
	// 	cout<<it->first<<": ";
	// 	for(int i = 0; i<it->second.size(); i++){
	// 		cout<<it->second[i]<<" ";
	// 	}
	// 	cout<<endl;
	// }
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

int StrToint(const string number_in_string) {
  int number = 0;
  istringstream(number_in_string) >> number;
  return number;
}