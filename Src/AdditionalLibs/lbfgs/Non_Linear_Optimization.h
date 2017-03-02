#pragma once
#include "lbfgs.h"
#include"../prereq.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include"lbfgspre.h"
class ADDITIONALLIBS_CLASS Non_Linear_Optimization
{
protected:
    std::vector<lbfgsfloatval_t> m_variables;
	lbfgs_parameter_t m_parameters;
	bool m_printprogress;
public:
	Non_Linear_Optimization();
	bool& GetPrintProgress();
	//return the number of iterations
	lbfgs_parameter_t& GetParameters();
	
	int  run();
	
protected:
	int progress(
		const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g,
		const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step,
		int n,
		int k,
		int ls
	);
	virtual lbfgsfloatval_t evaluate(
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step
		) = 0;

	static lbfgsfloatval_t _evaluate(
		void *instance,
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step
	);

	static int _progress(
		void *instance,
		const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g,
		const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step,
		int n,
		int k,
		int ls
	);
};