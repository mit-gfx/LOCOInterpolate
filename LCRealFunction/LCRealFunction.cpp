#include <iostream>
#include "LCRealFunction.h"
#include "LCError.h"
#include "LCRealFunctionValue.h"


LCRealFunction::LCRealFunction(int nParams, bool linear)
{
	nParams_ = nParams;
	linear_ = linear;
	for (int i = 0; i < nParams; i++)
	{
		ranges_.push_back(0);
		ranges_.push_back(1);
	}
	if (!linear)
	{
		srand(5);
		nSinSums_ = 10; // the larger the more complex the function 
		amplitudeValues_ = Eigen::MatrixXd::Random(nParams_, nSinSums_);
		phaseValues_ = Eigen::MatrixXd::Random(nParams_, nSinSums_);
	}
}
int LCRealFunction::getNParams() const
{
	return nParams_;
}
double LCRealFunction::getMinRange(int iParam) const
{
	if (iParam >= nParams_)
	{
		std::cout << "parameter out of range" << std::endl;
	}
	return ranges_[iParam * 2];
}

double LCRealFunction::getMaxRange(int iParam) const
{
	if (iParam >= nParams_)
	{
		std::cout << "parameter out of range" << std::endl;
	}
	return ranges_[iParam * 2 + 1];
}

LCError LCRealFunction::evalDeriv(const std::vector<double> &params, int direction, LCFunctionDeriv **result)
{
	return LCError("not implemented");
}

LCError LCRealFunction::evalShapeInfo(const std::vector<double> &params, LCFunctionValue **shapeInfo)
{
	LCError err;
	//LCErrorReturn(validateParameters(params));
	double val = 1;
	for (int i = 0; i < nParams_; i++)
	{
		double fi = 0;
		for (int j = 0; j < nSinSums_; j++)
		{
			fi += amplitudeValues_(i,j)*sin(3*j* params[i] + phaseValues_(i,j));
		}
		val += fi;
	}
	if (linear_)
	{
		val = 1;
		for (int i = 0; i < nParams_; i++)
		{
			val += (i + 1)*params[i];
		}
	}

	*shapeInfo = new LCRealFunctionValue(val);

	return err;

}