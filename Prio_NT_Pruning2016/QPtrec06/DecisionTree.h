#include <string>
#include <vector>
#include <iostream>
#include <cstring>
#include "ListIterator.h"

using namespace std;

struct MinHeapItem
{
	unsigned docID;
	float score;
};

class MinHeap
{
public:
	unsigned size;
	float smallest;
	unsigned n;

	MinHeap(unsigned s);
	unsigned push(unsigned docID,float score);
	bool remove(unsigned docID);
	bool increase(unsigned docID,float score);
	bool increaseAdd(unsigned docID,float score);
	unsigned pop(float &score);
	MinHeap* copy();
	~MinHeap();
private:
	struct MinHeapItem *items;
	
	void inline heapify(unsigned t);
};

class Tier
{
public:
	vector<pair<unsigned,unsigned*> > record;

	Tier(unsigned n) {num=n;}
	void addRecord(unsigned docID,unsigned tf[]);
	inline unsigned size();
	~Tier();
private:
	unsigned num;//number of terms in this node
};

struct termReter
{
	unsigned curDocID;
	lptr *PR;
	unsigned termID;
	// float QF;
	float grade;
	float threshold;
    unsigned blockNum;
    // unsigned curBlock;
    // unsigned* blockDocID;
    // float* blockMaxScore;
    // BlockInfo *BI;
};

struct DecisionTreeNode
{
	DecisionTreeNode *parent;
	DecisionTreeNode *left;
	DecisionTreeNode *right;
	Tier *content;
	vector<float> termScores;
	int blockNum;
	float score;
};

struct BlockNode
{
	Tier *content;
	vector<float> termScores;
	//int blockNum;
	unsigned code;
	float score;
};

inline bool comReter(termReter t1,termReter t2) {return t1.grade>t2.grade;}
inline bool comDecisionTreeNode(DecisionTreeNode *n1,DecisionTreeNode *n2) {return n1->score > n2->score;}

class DecisionTree
{
public:
	DecisionTree(termReter *r,unsigned n,unsigned rn);
	inline void putDoc(unsigned curDoc,termReter *r);
	void getGrade(vector<pair<unsigned,float> > &v);
	inline DecisionTreeNode** getBlockList();
	inline int getBlockNum();
	~DecisionTree();
private:
	unsigned num,retNum;
	int blockNum;
	DecisionTreeNode **blockList;
	unsigned *blockSize;
	DecisionTreeNode *head;//only a pointer to root is enough to reverse the tree

	void deleteNode(DecisionTreeNode *n);
	void deleteLeaf(DecisionTreeNode *n);
};

class BlockManager
{
public:
	BlockManager(termReter *r, unsigned n,unsigned rn);
	inline void putDoc(unsigned curDoc,termReter *r);
	void getGrade(vector<pair<unsigned,float> > &v);
	inline BlockNode** getBlockList() {return blockList;}
	inline int getBlockNum() {return blockNum;}
	//inline bool decideCut(unsigned n) {return (cutNum[n]<=0);}
	inline float getTh() {return threshold;}
	~BlockManager();
private:
	unsigned num,retNum;
	int blockNum;
	BlockNode **blockList;
	unsigned *blockSize;
	float threshold;
	//int *cutNum;

	int findScore(float score);
};

class RetManager
{
public:
	RetManager(vector<string>& word_l, termsMap& lex, termsCache& Cache, termsMap& termsInCache);
	uint retrieval(unsigned retNum, unsigned* pages, profilerC& p);
	inline int getEvalCounter() {return evalCounter;}
	~RetManager();
	unsigned* retDocID;
private:
	// IndexReader *theIndex;
	unsigned num;
	string topicNum;

	termReter *r;
	
	float* retDocScore;
	int retN;
	int curDoc;
	MinHeap* topDocs;
	int evalCounter;
	BlockManager *BM;

	inline unsigned findNextDoc();
	inline float grade(unsigned* pages);
};