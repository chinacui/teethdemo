#include"linear_algebra_alg.h"



Eigen::VectorXd ClinearAlgebraAlg::SolveLdlt(const Eigen::SparseMatrix<double>& CoeffMat, const Eigen::VectorXd& right, double pinvtoler)
{
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(CoeffMat);
	Eigen::VectorXd X_0 = ldlt.permutationP() * right;
	Eigen::VectorXd X_1 = ldlt.matrixL().solve(X_0);
	Eigen::VectorXd X_2(ldlt.vectorD().size());
	X_2.setZero();
	for (int i = 0; i < ldlt.vectorD().size(); ++i)
		if (abs(ldlt.vectorD()(i)) > pinvtoler)
			X_2[i] = X_1[i] / ldlt.vectorD()(i);
	Eigen::VectorXd X_3 = ldlt.matrixU().solve(X_2);
	return ldlt.permutationPinv() * X_3;
}