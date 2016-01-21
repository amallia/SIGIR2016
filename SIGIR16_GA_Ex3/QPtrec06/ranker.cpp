/*
 * ranker.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: Qi
 */

#include "ranker.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include <cstring>
#include <string>
#include <fstream>
#include <algorithm>
#include <set>
using namespace std;

void ranker::operator() (int qn, int qid, toplayers& tls, toplayers& otls, pairlists& pls, pairlists& opls) {

	// const int p_l = opls.cutoffs.size();
	// const int t_l = otls.cutoffs.size();
	// const int a_l = opls.cutoffs.size() + otls.cutoffs.size();

	const int p_l = pls.cutoffs.size();
	const int t_l = tls.cutoffs.size();
	const int a_l = pls.cutoffs.size() + tls.cutoffs.size();


	int p_e[10][p_l]; 
	int t_e[10][t_l]; 

	const int buffer_size = 1024*100;
	string queryline;
	ulong results[10];
	int known_l[10];
	for(int i = 0; i < 10; i++){
		known_l[i] = 0;
	}

    char * docid_s;
	// int docid;
	int count = 0;
	int test;

	// Get the top-10 results from Complicated ranking func

	ifstream input_stream;
	// input_stream.open("../../../Dropbox/WSDM_Index_Script/result_log");  //for or
	// input_stream.open("/data/qw376/WSDM_Index_Script/500_cr");  //for cr
	input_stream.open(CONSTS::cr_results.c_str());  //for cr


	for(int i=0; i < qn-1; ++i){
        input_stream.ignore(numeric_limits<streamsize>::max(),'\n');
    }

	getline(input_stream, queryline);

	char * ql = new char [queryline.length()+1];
    std::strcpy (ql, queryline.c_str());

	// cout<<queryline<<endl;


	docid_s = strtok (ql," ");

	while (docid_s != NULL)
  	{
    	// printf ("%s\n",docid_s);
 		results[count] = atol(docid_s);   	
 		count++;
    	docid_s = strtok (NULL, " ,.-");
    	if(count == 10)
    		break;
  	}

  	input_stream.close();
    // Done, get the top-10 results from Complicated ranking func

  	// for(int b=0; b<count; b++){
  	// 	cout<<results[b]<<endl;
  	// }
  	/*-----------------------*/

 //  	ofstream did_index;
	// string did_dir = "/home/qw376/depth_index/dids";
 //  	did_index.open(did_dir.c_str(), ofstream::app);
 //  	did_index << qid << ":";
 //  	for(int b=0; b<count; b++){
 //  		did_index << results[b] << " ";
 //  	}
 //  	did_index << endl;
 //  	did_index.close();

    //Output the number of avaliable structures
  	ofstream did_index;
	string did_dir = "/home/qw376/NumOfValidStructures_M9/NumOfValidStructures_M9";
  	did_index.open(did_dir.c_str(), ofstream::app);

  	int NumOfValidStructures = 0;

  	if(pls.cutoffs.size()>0){
  		for(int k = 0; k<pls.cutoffs.size(); k++){
  			if(pls.cutoffs[k]>0)
  				NumOfValidStructures++;
  		}
  	}

  	if(tls.cutoffs.size()>0){
  		for(int k = 0; k<tls.cutoffs.size(); k++){
  			if(tls.cutoffs[k]>0)
  				NumOfValidStructures++;
  		}
  	}

  	cout << NumOfValidStructures <<endl;
  	did_index << NumOfValidStructures <<endl;
  	did_index.close();
  	//Done, Output the number of avaliable structures

  	cout<<qid<<": "<<endl;

  	/*---load the pair lists into map and outputs the depth----*/
  	// this is for the pairs in ideal case
  	// if(opls.cutoffs.size()>0){

	  // for(int k = 0; k<opls.cutoffs.size(); k++){

			// string did_dir = CONSTS::intersection_index_did + opls.pairnames[k];

		 //    const int length = opls.cutoffs[k];

		 //    // cout<<opls.pairnames[k]<<": "<<length<<endl;

		 //    FILE * loaddids = fopen(did_dir.c_str(), "r");
		 //    if(loaddids==NULL){
		 //            cout<<"can't find file for did"<<endl;
		 //    }
		 //    vector<uint> t_a(length, 0);
		 //    fread(&t_a[0], sizeof(uint), length, loaddids);
		 //    fclose(loaddids);


		 //    for(int i = 0; i < 10; i++){
		 //                    int j;
		 //                    for(j = 0; j<length; j++){
		 //                            if (results[i] == t_a[j])
		 //                                    break;
		 //                    }

		 //            if(j!=length){
		 //                    // cout<<results[i]<<" is at depth: "<<j+1<<" of "<<opls.pairnames[k]<<endl;
		 //                    p_e[i][k] = j+1;

		 //            }else{  
		 //                    // cout<<results[i]<<" is at depth: -"<<" of "<<opls.pairnames[k]<<endl;
		 //                    p_e[i][k] = 0;
		 //            }
		 //    }
   //      }
   //  }

  	// this is for the pairs in real case
  	if(pls.cutoffs.size()>0){

	  for(int k = 0; k<pls.cutoffs.size(); k++){

  		string did_dir = CONSTS::intersection_index_did + pls.pairnames[k];

                const int length = pls.cutoffs[k];

                // cout<<lps[k] -> term<<": "<<size<<endl;

                FILE * loaddids = fopen(did_dir.c_str(), "r");
                if(loaddids==NULL){
                        cout<<"can't find file for did"<<endl;
                }
                vector<uint> t_a(length, 0);
                fread(&t_a[0], sizeof(uint), length, loaddids);
                fclose(loaddids);


                for(int i = 0; i < 10; i++){
                                int j;
                                for(j = 0; j<length; j++){
                                        if (results[i] == t_a[j])
                                                break;
                                }

                        if(j!=length){
                                // cout<<results[i]<<" is at depth: "<<j+1<<" of "<<lps[k]->term<<endl;
                                p_e[i][k] = j+1;

                        }else{  
                                // cout<<results[i]<<" is at depth: -"<<" of "<<lps[k]->term<<endl;
                                p_e[i][k] = 0;
                        }
                }
        }
    }

//    cout<<"pair done"<<endl;

	/*---load the single term lists into map and outputs the depth (From Binary lists)----*/
    // this is for the single in ideal case
	// if(otls.cutoffs.size()>0){

	// 	for(int k = 0; k<otls.cutoffs.size(); k++){

	// 		string did_dir = CONSTS::top_layer_index_did.c_str() + otls.terms[k];

	// 		const int length = otls.cutoffs[k];
	// 		// cout<<otls.terms[k]<<": "<<length<<endl;

	// 		FILE * loaddids = fopen(did_dir.c_str(), "r");
	// 		if(loaddids==NULL){
	// 			cout<<"can't find file for did"<<endl;
	// 		}
	// 		vector<uint> t_a(length, 0);
	// 		fread(&t_a[0], sizeof(uint), length, loaddids);
	// 		fclose(loaddids);

	// 		for(int i = 0; i < 10; i++){
	// 				int j;
	// 				for(j = 0; j<length; j++){
	// 					if (results[i] == t_a[j])
	// 						break;
	// 				}

	// 			if(j!=length){
	// 				// cout<<results[i]<<" is at depth: "<<j+1<<" of "<<otls.terms[k]<<endl;
	// 				t_e[i][k] = j+1;
				
	// 			}else{
	// 				// cout<<results[i]<<" is at depth: -"<<" of "<<otls.terms[k]<<endl;
	// 				t_e[i][k] = 0;
	// 			}

	// 		}

	//   	}
	// }

	// this is for the single in real case
	if(tls.cutoffs.size()>0){

		for(int k = 0; k<tls.cutoffs.size(); k++){

			string did_dir = CONSTS::top_layer_index_did.c_str() + tls.terms[k];

			const int length = tls.cutoffs[k];
//			cout<<otls.terms[k]<<": "<<length<<endl;

			FILE * loaddids = fopen(did_dir.c_str(), "r");
			if(loaddids==NULL){
				cout<<"can't find file for did"<<endl;
			}
			vector<uint> t_a(length, 0);
			fread(&t_a[0], sizeof(uint), length, loaddids);
			fclose(loaddids);

			for(int i = 0; i < 10; i++){
					int j;
					for(j = 0; j<length; j++){
						if (results[i] == t_a[j])
							break;
					}

				if(j!=length){
//					cout<<results[i]<<" is at depth: "<<j+1<<" of "<<otls.terms[k]<<endl;
					t_e[i][k] = j+1;
				
				}else{
//					cout<<results[i]<<" is at depth: -"<<" of "<<otls.terms[k]<<endl;
					t_e[i][k] = 0;
				}

			}

	  	}
	}
	/*---------*/

//	cout<<"single done"<<endl;


	// ofstream index_depth;
	// // const string dir = "/data/qw376/experiments/depth_index";
	// // index_depth.open(CONSTS::index_depth.c_str(), ofstream::app);
	// const string dir = "/home/qw376/depth_index/m9_real_all_mix";
	// // index_depth.open(dir.c_str(), ofstream::app);
	// // cout<<qid<<": "<<endl;
	// // index_depth<<qid<<": "<<endl;
	// for(int i = 0; i < 10; i++){
	// 	for(int n = 0; n < t_l; n++){
	// 		// cout<<otls.terms[n]<<" "<<otls.cutoffs[n]<<" "<<t_e[i][n]<<" ";
	// 		// cout<<tls.terms[n]<<" "<<tls.cutoffs[n]<<" "<<t_e[i][n]<<" ";
	// 		// index_depth<<otls.terms[n]<<" "<<otls.cutoffs[n]<<" "<<t_e[i][n]<<" ";
	// 		// index_depth<<tls.terms[n]<<" "<<tls.cutoffs[n]<<" "<<t_e[i][n]<<" ";
	// 		// index_depth<<t_e[i][n]<<" ";
	// 	}
	// 	for(int n = 0; n < p_l; n++){
	// 		// cout<<opls.pairnames[n]<<" "<<opls.cutoffs[n]<<" "<<p_e[i][n]<<" ";
	// 		// cout<<pls.pairnames[n]<<" "<<pls.cutoffs[n]<<" "<<p_e[i][n]<<" ";
	// 		// index_depth<<opls.pairnames[n]<<" "<<opls.cutoffs[n]<<" "<<p_e[i][n]<<" ";
	// 		// index_depth<<pls.pairnames[n]<<" "<<pls.cutoffs[n]<<" "<<p_e[i][n]<<" ";
	// 		// index_depth<<p_e[i][n]<<" ";
	// 	}
	// 	// cout<<endl;
	// 	// index_depth<<endl;
	// }
	// // index_depth.close();

}
