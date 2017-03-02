#include"Non_Linear_Optimization.h"
#include "lbfgs.h"
#include"Non_Linear_Optimization.h"

Non_Linear_Optimization::Non_Linear_Optimization() : m_parameters(_defparam)
{
	m_printprogress = true;
}
bool& Non_Linear_Optimization::GetPrintProgress()
{
	return m_printprogress;
}
//return the number of iterations
lbfgs_parameter_t& Non_Linear_Optimization::GetParameters()
{
	return m_parameters;
}

int  Non_Linear_Optimization::run()
{
	/*
	Start the L-BFGS optimization; this will invoke the callback functions
	evaluate() and progress() when necessary.
	*/
	lbfgsfloatval_t fx = 0;
	int numOfIterations = 0;
	int ret = lbfgs(m_variables.size(), &m_variables[0], &fx, _evaluate, _progress, this, &m_parameters, &numOfIterations);

	/* Report the result. */
	if (m_printprogress)
		std::cerr << "L-BFGS optimization terminated with status code = " << ret << std::endl;

	return numOfIterations;
}

int Non_Linear_Optimization::progress(
	const lbfgsfloatval_t *x,
	const lbfgsfloatval_t *g,
	const lbfgsfloatval_t fx,
	const lbfgsfloatval_t xnorm,
	const lbfgsfloatval_t gnorm,
	const lbfgsfloatval_t step,
	int n,
	int k,
	int ls
)
{
	if (!m_printprogress)
		return 0;
	std::cerr << " Iteration #" << k << ": ";
	std::cerr << " Energy = " << fx << ", ";
	std::cerr << " step = " << step << ", ";
	std::cerr << " Error = " << gnorm / std::max(1.0, xnorm) << std::endl;
	for (int i = 0; i < n; ++i)
	{
		std::cerr << "x[" << i << "] = " << x[i] << ", ";
	}
	std::cerr << std::endl;
	std::cerr << "-----------------------------------------\n";

	return 0;
}

int Non_Linear_Optimization::_progress(
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
)
{
	return reinterpret_cast<Non_Linear_Optimization*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}

 lbfgsfloatval_t Non_Linear_Optimization::_evaluate(
	void *instance,
	const lbfgsfloatval_t *x,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
)
{
	if (reinterpret_cast<Non_Linear_Optimization*>(instance)->m_printprogress)
		std::cerr << "\tCall your evaluation function once..........\n";
	return reinterpret_cast<Non_Linear_Optimization*>(instance)->evaluate(x, g, n, step);
}