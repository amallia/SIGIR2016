#include "DecisionTree.h"

using namespace std;

RetManager::RetManager(vector<string>& word_l, termsMap& lex, termsCache& Cache, termsMap& termsInCache)
{
	int i;
	// theIndex = IR;
	num = word_l.size();// num = q->term.size();
	// topicNum = q->topicNum;
	// Initial arrays
	r = new termReter[num];

	std::vector<std::string> peers( word_l.size());
	for(size_t i=0; i<word_l.size(); i++) {
		if(termsInCache.find(word_l[i]) == termsInCache.end()){ //no term! try on-demand processing from rawindex
			// onDemandCall(word_l[i]);
			cout<<"The baby index doesn't contain this term"<<endl;
		}
		peers.push_back(Cache[ termsInCache[word_l[i]]].cpool.getName() ); //This is actually data[term].cpool.getName(), this is just adding the termpath before term. Note: termsInCache is termsMap(std::unordered_map<std::string, int>)
	}//This part of code is really confusing, but the cache part is smart

	// set initial values
	for(i=0;i<num;i++)
	{
		// r[i].QF = q->tf[i];
		termsMap::iterator it;
		it = lex.find(word_l[i]); //r[i].termID = theIndex->termLookup(theIndex->stem(q->term[i]));
		if(it != lex.end()){
			r[i].termID = it->second; 
			lptr* lpsi = &(Cache[termsInCache[word_l[i]]]);
			lpsi->open(peers,lex.find((word_l[i]))->second);
			r[i].PR = lpsi; //r[i].PR = theIndex->getPosting(r[i].termID);
			//r[i].maxScore = theIndex->getMaxScore(r[i].termID);
			float docN = CONSTS::MAXD;//float docN = theIndex->documentN;
			float DF = r[i].PR->unpadded_list_length;//float DF = r[i].PR->docCount;
			r[i].curDocID = r[i].PR->did;
			r[i].grade = log((docN-DF+0.5)/(DF+0.5)); //this is IDF in BM25
			cout << word_l[i] << " listlength: "<< r[i].PR->unpadded_list_length << endl;
		}
		else
		{	
			cout << word_l[i] << ": term not in the lex" <<endl;
			r[i].PR = NULL;
			r[i].curDocID = CONSTS::MAXD;
			r[i].grade = 0;
		}
	}
	sort(r,r+num,comReter); //sort the terms in termReter according to IDF
	// cout <<"test1"<<endl;
	// if(num>0) r[num-1].threshold = r[num-1].grade;
	// cout <<"test2"<<endl;
	// cout << num << endl;
	// for(i=num-2;i>=0;i--) {
	// 	cout << i << endl;
	// 	cout << r[i].grade << endl;
	// 	cout << r[i+1].threshold << endl;
	// 	r[i].threshold = r[i+1].threshold+r[i].grade;
	// }
	// cout <<"test3"<<endl;

	// set other values
	retDocID = NULL;
	retDocScore = NULL;
	retN = 0;
	curDoc = -1;
	topDocs = NULL;
	evalCounter=0;
}

RetManager::~RetManager()
{
	unsigned i;
	for(i=0;i<num;i++) if(r[i].PR!=NULL) delete r[i].PR;
	delete[] r;
	if(retDocID != NULL) delete[] retDocID;
	if(retDocScore!=NULL) delete[] retDocScore;
}

uint RetManager::retrieval(unsigned retNum, unsigned* pages, profilerC& p)
{
	// p.start(CONSTS::ALLQS);
	cout << retNum << endl;
	retDocID = new unsigned[retNum];
	retDocScore = new float[retNum];
	// p.start(CONSTS::STEP1);
	DecisionTree *DT = new DecisionTree(r,num,retNum);
	theHead = DT->getHead();
	// p.end(CONSTS::STEP1);
	// cout << "qp processing" << endl;
	int i,l,n;
	// p.start(CONSTS::STEP2);
	curDoc = findNextDoc();
	uint Totaldocs_decompressed = 1;
	while(curDoc < CONSTS::MAXD)
	{	
		for(i=0;i<num;i++) if(r[i].curDocID<curDoc) r[i].curDocID = r[i].PR->nextGEQ(curDoc);
		DT->putDoc(curDoc,r);;//DT->putDoc(curDoc,r);
		curDoc = findNextDoc();
		Totaldocs_decompressed ++;
	}
	// cout <<"Totaldocs_decompressed: "<< Totaldocs_decompressed <<endl;
	// p.end(CONSTS::STEP2);
	// cout << "qp processing done" << endl;

	// cout << "look up" << endl;
	const float okapiK1=1.2;
	const float okapiB=0.2;
	DecisionTreeNode **blockList = DT->getBlockList();
	vector<pair<unsigned,unsigned*> >::iterator it;
	// p.start(CONSTS::STEP3);
	cout << "Number of Blocks: " <<DT->getBlockNum()<<endl;
	for(i=0;i<DT->getBlockNum();i++)
	{	
		cout << "i: "<< i << endl;
		n=blockList[i]->termScores.size();
		cout << "n: "<< n << endl;
		int theSize = blockList[i]->content->record.size();
		cout<< "theSize of current lock: "<<theSize << endl;
		if(theSize+retN>retNum) theSize = retNum-retN;
		cout<< "TopK: "<< retNum <<" retN: "<<retN<<" theSize needed to reach topK: " <<theSize<<endl;
		if(theSize==0) continue;
		topDocs = new MinHeap(theSize);
		cout<< "start doing look up in block: "<<i<< endl;
		for(it=blockList[i]->content->record.begin();it!=blockList[i]->content->record.end();it++)
		{														
			unsigned* theTF=it->second;
			float score=0;
			float docLength = pages[it->first];//float docLength = theIndex -> getDocLength(it->first);
			for(l=0;l<n;l++)
			{	
				float tf = theTF[l];
				float weight = ((okapiK1+1.0)*tf) / (okapiK1 * (1.0-okapiB+okapiB * docLength /CONSTS::AVGD)+tf);//float weight = ((okapiK1+1.0)*tf) / (okapiK1*(1.0-okapiB+okapiB*docLength/theIndex->docLengthAvg)+tf);
				score+=weight*blockList[i]->termScores[l];
			}
			if(score > topDocs->smallest) topDocs->push(it->first,score);
			evalCounter++;
		}	
		cout<< "Doc evaluated by far: "<<evalCounter<<endl;
		cout<< "doing look up in block: "<<i<<" done"<< endl;
		// p.end(CONSTS::STEP3);
		for(l=retN+theSize-1;l>=retN;l--)
		{
			// cout << l << endl;
			retDocID[l] = topDocs->pop(retDocScore[l]);
		}
		retN+=theSize;
		cout << "retN: "<<retN<<endl;
		delete(topDocs);
	}
	delete(DT);
	cout << "look up done" << endl;
	// p.end(CONSTS::ALLQS); 
	cout << "docs decompressed: "<< Totaldocs_decompressed << endl;
	return Totaldocs_decompressed;
}

unsigned inline RetManager::findNextDoc()
{
	unsigned minDoc = CONSTS::MAXD;
	unsigned i,l;
	float s=0.0f;

	DecisionTreeNode* theNode = theHead;

	for(i=0;i<num;i++)
	{
		if(theNode == NULL) break;
		// cout<<i<<" "<<r[i].curDocID<<" "<<curDoc<<endl;
		if(r[i].curDocID<=curDoc) 
		{
			// r[i].curDocID = r[i].PR->nextGEQ(curDoc+1);
			r[i].curDocID = r[i].PR->nextGEQ(r[i].curDocID+1);
			// if(r[i].PR->nextRecord()) r[i].curDocID = r[i].PR->getDocID();
			// else r[i].curDocID = CONSTS::MAXD;
		}
		// cout<<i<<" "<<r[i].curDocID<<" "<<curDoc<<endl;
		if(minDoc>r[i].curDocID) minDoc = r[i].curDocID;
		theNode = theNode->right;
	}
	// exit(0);
	return minDoc;
}

float inline RetManager::grade(unsigned* pages)
{
	const float okapiB = 0.2;
	const float okapiK1 = 1.2;
	const float  okapiK3 = 7;
	float score = 0;
	float docN = CONSTS::MAXD;//float docN = theIndex -> documentN;
	float docLength = pages[curDoc];//float docLength = theIndex -> getDocLength(curDoc);
	float docLengthAvg = CONSTS::AVGD;//float docLengthAvg = theIndex -> docLengthAvg;
	unsigned i;
	for(i=0;i<num;i++)
	{
		if(curDoc == r[i].curDocID)
		{
			float DF = r[i].PR->unpadded_list_length;//float DF = r[i].PR->docCount;
			float tf = r[i].PR->getFreq();//float tf = r[i].PR->getTF();
			float idf = log((docN-DF+0.5)/(DF+0.5));
			//float idf = log((1.0+docN)/DF);
			//float weight = (1.0+log(1+log(tf)))/(1.0-okapiB+okapiB*docLength/docLengthAvg);
			float weight = ((okapiK1+1.0)*tf) / (okapiK1*(1.0-okapiB+okapiB*docLength/docLengthAvg)+tf);
			//float tWeight = ((okapiK3+1)*r[i].QF)/(okapiK3+r[i].QF);
			//score+=idf*weight*tWeight;
			score+=idf*weight;
		}
	}
	/*unsigned i;
	float score=0;
	const float mu = 2500;
	float docLength = theIndex -> getDocLength(curDoc);
	for(i=0;i<num;i++)
	{
		if(curDoc == curDocID[i])
		{
			float termPro = PR[i]->totalTF/theIndex->termNTotal;
			float tf = PR[i]->getTF();
			score+=log((tf+mu*termPro)/(docLength+mu));
		}
		else if(PR[i] != NULL)
		{
			float termPro = PR[i]->totalTF/theIndex->termNTotal;
			score+=log((mu*termPro)/(docLength+mu));
		}
	}*/
	//float score = 1;
	evalCounter++;
	return score;
}

DecisionTree::DecisionTree(termReter *r,unsigned n,unsigned rn)
{
	num=n;//number of terms
	retNum=rn;//topK
	vector<DecisionTreeNode*> parents,parentsNew;
	int i;
	vector<DecisionTreeNode*>::iterator it;

	//Generate the tree structure	
	if(num<1) {head=NULL;return;}
	head = new DecisionTreeNode;
	head->parent = NULL;
	head->score = 0;
	parents.push_back(head); // head is an empty node
	for(i=0;i<num;i++)
	{
		for(it=parents.begin();it!=parents.end();it++)
		{
			DecisionTreeNode *node = new DecisionTreeNode;
			(*it)->left = node;
			node->parent = *it;
			node->score = (*it)->score + r->grade;// here is to add the IDF to the left nodes in this level
			node->termScores = (*it)->termScores;
			node->termScores.push_back(r->grade);
			parentsNew.push_back(node);
			node = new DecisionTreeNode;
			(*it)->right = node;
			node->parent = *it;
			node->score = (*it)->score;
			node->termScores = (*it)->termScores;
			parentsNew.push_back(node);
			(*it)->content=NULL;
			(*it)->blockNum = -1;
		}
		parents=parentsNew; //renew parents with new leaf nodes, in the end, parenets are all the leaf nodes
		parentsNew.clear();
		r++;// it stops when r increases up to num(num of terms in a query)
	}
	//Generate the blockList
	blockList = new DecisionTreeNode*[parents.size()];
	blockSize = new unsigned[parents.size()];
	blockNum=0;
	for(it=parents.begin();it!=parents.end();it++)
	{
		(*it)->left=NULL;
		(*it)->right=NULL;
		(*it)->content = new Tier((*it)->termScores.size()); //Tier is a vector of <did, tf>
		blockList[blockNum++]=*it;//blocklist store the pointer of all leaf nodes 
	}
	sort(blockList,blockList+blockNum,comDecisionTreeNode); //sort the nodes according accumulated IDF
	for(i=0;i<blockNum;i++) {blockList[i]->blockNum=i;blockSize[i]=0;}//blocksize means number of dids in the block
}
	
void DecisionTree::deleteNode(DecisionTreeNode *n)
{
	if(n->left!=NULL) deleteNode(n->left);
	if(n->right!=NULL) deleteNode(n->right);
	if(n->content!=NULL) delete n->content;
	delete n;
}

void DecisionTree::deleteLeaf(DecisionTreeNode *n)
{
	if(n->content!=NULL) delete n->content;
	if(n->parent!=NULL)
	{
		DecisionTreeNode *parent = n->parent;
		if(parent->left == n) parent->left = NULL;
		if(parent->right == n) parent ->right = NULL;
		if(parent->left == NULL && parent->right == NULL) deleteLeaf(parent);
	}
	delete n;
}

inline DecisionTreeNode** DecisionTree::getBlockList() {return blockList;}
inline int DecisionTree::getBlockNum() {return blockNum;}
inline DecisionTreeNode* DecisionTree::getHead() {return head;}

inline void DecisionTree::putDoc(unsigned curDoc,termReter *r)
{
	DecisionTreeNode* curNode = head;
	unsigned tf[num];
	int n=0;
	while(curNode!=NULL && curNode->content==NULL)
	{	
		if(curDoc==r->curDocID)
		{
			curNode = curNode->left;
			tf[n++] = r->PR->getFreq(); //tf[n++] = r->PR->getTF();	
			// cout <<"TermID: "<<r->termID<<" did: "<<curDoc<< " freq:" <<r->PR->getFreq()<<endl;	
		}
		else
		{
			curNode = curNode->right;
		}
		r++;
	}//it will go down to one leaf node
	if(curNode!=NULL) 
	{
		curNode->content->addRecord(curDoc,tf);// whichever term contains the did, record the freq
		n=curNode->blockNum;
		int i;
		for(i=n;i<blockNum;i++) blockSize[i]++;//this blockSize is accumulated blocksize, any node after this node, size increase by 1
		while(n<=blockNum-2 && blockSize[blockNum-2]>=retNum)// this is pruning here, blockNum-2 means start from the third last one. You can't start from the second last one, since it may contain useful did. The last one has no dids.
		{
			blockNum--;
			deleteLeaf(blockList[blockNum]);
		}
	}	
}		

DecisionTree::~DecisionTree()
{
	//cout<<"blockSize = "<<blockSize[0]<<endl;
	deleteNode(head);
	delete[] blockSize;
	delete[] blockList;
}

void Tier::addRecord(unsigned docID,unsigned tf[])
{
	unsigned *theTF = new unsigned[num];
	memcpy(theTF,tf,sizeof(unsigned)*num);
	record.push_back(pair<unsigned,unsigned*>(docID,theTF));
}

inline unsigned Tier::size() {return record.size();}

Tier::~Tier()
{
	vector<pair<unsigned,unsigned*> >::iterator it;
	for(it=record.begin();it!=record.end();it++) delete[] (it->second);
}

MinHeap::MinHeap(unsigned s)
{
	size=s;
	items=new MinHeapItem[size+1];
	n=0;
	smallest=-2000000000;
}

void inline MinHeap::heapify(unsigned t)
{
	unsigned left= t*2;
	unsigned right= t*2+1;
	unsigned minT=t;
	if(left<=n && items[left].score < items[minT].score) minT=left;
	if(right<=n && items[right].score < items[minT].score) minT=right;
	if(minT!=t)
	{
		unsigned tempD=items[t].docID;
		float tempS=items[t].score;
		items[t].docID=items[minT].docID;
		items[t].score=items[minT].score;
		items[minT].docID=tempD;
		items[minT].score=tempS;
		heapify(minT);
	}
}

unsigned MinHeap::push(unsigned docID,float score)
{
	unsigned minK=0xFFFFFFFF;
	if(n>=size) //heap full
	{
		minK = items[1].docID;
		items[1].docID=docID;
		items[1].score=score;
		heapify(1);
	}
	else
	{	
		n++;
		int p=n;
		while(p > 1 && items[p/2].score > score)
		{
			items[p].docID=items[p/2].docID;
			items[p].score=items[p/2].score;
			p/=2;
		}
		items[p].docID = docID;
		items[p].score=score;
	}
	if(n>=size) {smallest=items[1].score;}
	return minK;
}

bool MinHeap::increase(unsigned docID,float score)
{
	unsigned i;
	for(i=1;i<=n;i++)
	{
		if(items[i].docID==docID)
		{
			items[i].score=score;
			heapify(i);
			if(n>=size) smallest=items[1].score;
			return true;
		}
	}
	return false;
}

bool MinHeap::increaseAdd(unsigned docID,float score)
{
	unsigned i;
	for(i=1;i<=n;i++)
	{
                if(items[i].docID==docID)
                {
                        items[i].score += score;
                        heapify(i);
                        if(n>=size) smallest=items[1].score;
                        return true;
                }
        }
        return false;
}

MinHeap* MinHeap::copy()
{
	MinHeap* MH=new MinHeap(size);
	MH->n=n;
	MH->smallest=smallest;
	int i;
	for(i=1;i<=n;i++)
	{
		MH->items[i].docID=items[i].docID;
		MH->items[i].score=items[i].score;
	}
	return MH;
}

bool MinHeap::remove(unsigned docID)
{
	unsigned i;
	for(i=1;i<=n;i++)
	{
		if(items[i].docID==docID)
		{
			items[i].docID=items[n].docID;
			items[i].score=items[n].score;
			n--;
			heapify(i);
			if(i==1) smallest=items[1].score;
			return true;
		}
	}
	return false;
}

unsigned MinHeap::pop(float &score)
{
	if(n<1) return -1;
	unsigned docID=items[1].docID;
	score=items[1].score;
	if(n>1)
	{
		items[1].docID=items[n].docID;
		items[1].score=items[n].score;
		n--;
		heapify(1);
	}
	else n--;
	return docID;
}

MinHeap::~MinHeap()
{
	delete[] items;
}