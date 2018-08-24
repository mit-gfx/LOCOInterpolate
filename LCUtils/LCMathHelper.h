#ifndef OC_MATH_HELPER
#define OC_MATH_HELPER

#include <vector>
#include <Eigen/Dense>

class LCError;

class LCMathHelper
{
public:
	static LCError computeCombinations(int nDims, int nValues,
		std::vector<Eigen::VectorXi> *result);
	static LCError computeCornerCombinations(int nDims, int nValues,
		std::vector<Eigen::VectorXi> *result);
	static void eigen2StdVector(const Eigen::VectorXd &eigenVec, std::vector<double> *stdVec);
	static void std2EigenVector(const std::vector<double> &stdVec, Eigen::VectorXd *eigenVec);
	static Eigen::VectorXd removeDimension(const Eigen::VectorXd &orig, int i);

	static const double EPSILON;
};

#endif

