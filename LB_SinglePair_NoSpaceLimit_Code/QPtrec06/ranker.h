/*
 * ranker.h
 *
 *  Created on: Feb 13, 2015
 *      Author: Qi
 */

#ifndef RANKER_H_
#define RANKER_H_

#include "PostingOriented_BMW.h"

class ranker{
private:
	unsigned int* pages;

public:
	ranker(unsigned int* pgs) : pages(pgs) {}
	void operator()(int qn, int qid, toplayers& tls,toplayers& otls, pairlists& pls, pairlists& opls);
};

#endif /* RANKER_H_ */