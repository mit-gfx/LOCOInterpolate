#include "LCMathHelper.h"
#include "LCError.h"

LCError LCMathHelper::computeCombinations(int nDims, int nValues,
	std::vector<Eigen::VectorXi> *result)
{
	
	if ((nDims < 1) || (nValues < 1))
	{
		return LCError("cannot compute combinations for zero dims or values");
	}
	int nCombinations = pow(nValues, nDims);
	result->reserve(nCombinations);
	for (int i = 0; i < nCombinations; i++)
	{
		int num = i;
		Eigen::VectorXi combination(nDims);
		for (int d = 0; d < nDims; d++)
		{
			int current = num % nValues;
			combination(d) = current;
			num = (num - current) / nValues;
		}
		result->push_back(combination);
	}
	return LCError();
}

LCError LCMathHelper::computeCornerCombinations(int nDims, int nValues,
	std::vector<Eigen::VectorXi> *result)
{
	LCError err;
	std::vector<Eigen::VectorXi> fullCombiantions;
	LCErrorReturn(LCMathHelper::computeCombinations(nDims, nValues, &fullCombiantions));
	for (auto fullCombiation : fullCombiantions)
	{
		bool isCorner = 0;
		for (int i = 0; i < nDims; i++)
		{
			if ((fullCombiation(i) == 0) || (fullCombiation(i) == nValues - 1))
			{
				isCorner = true;
			}
		}
		if (isCorner)
		{
			result->push_back(fullCombiation);
		}
	}
	return err;
}


void LCMathHelper::eigen2StdVector(const Eigen::VectorXd &eigenVec, std::vector<double> *stdVec)
{
	stdVec->clear();
	stdVec->reserve(eigenVec.size());
	for (int i = 0; i < eigenVec.size(); i++)
	{
		stdVec->push_back(eigenVec(i));
	}
}

void LCMathHelper::std2EigenVector(const std::vector<double> &stdVec, Eigen::VectorXd *eigenVec)
{
	Eigen::VectorXd vec(stdVec.size());
	for (int i = 0; i < stdVec.size(); i++)
	{
		vec(i) = stdVec[i];
	}
	*eigenVec = vec;
}


Eigen::VectorXd LCMathHelper::removeDimension(const Eigen::VectorXd &orig, int i)
{
	int projectedSize = orig.size() - 1;
	Eigen::VectorXd projected(projectedSize);
	projected.head(i) = orig.head(i);
	projected.tail(projectedSize - i) = orig.tail(projectedSize - i);
	return projected;
}

const double LCMathHelper::EPSILON = 0.000001;