/*
 * exhaustiveOR.h
 *
 *  Created on: May 17, 2012
 *      Author: constantinos
 */

#ifndef AND_H_
#define AND_H_

#include "PostingOriented_BMW.h"

class ExhaustiveAnd{
private:
	unsigned int* pages;

public:
	ExhaustiveAnd(unsigned int* pgs) : pages(pgs) {}
	void operator()(lptrArray& lps, const int topk, QpResult* r);
};

#endif /* AND_H_ */
