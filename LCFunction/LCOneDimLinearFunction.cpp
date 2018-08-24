#include <iostream>
#include "LCOneDimLinearFunction.h"
#include "LCError.h"
#include "LCRealFunctionValue.h"


LCOneDimLinearFunction::LCOneDimLinearFunction(int nParams)
{
	nParams_ = nParams;
	ranges_ = { 1, 5, 4, 6, .1, .6 };
}

int LCOneDimLinearFunction::getNParams() const
{
	return nParams_;
}

double LCOneDimLinearFunction::getMinRange(int iParam) const
{
	if (iParam >= nParams_)
	{
		std::cout << "parameter out of range" << std::endl;
	}
	return ranges_[iParam * 2];
}

double LCOneDimLinearFunction::getMaxRange(int iParam) const
{
	if (iParam >= nParams_)
	{
		std::cout << "parameter out of range" << std::endl;
	}
	return ranges_[iParam * 2 + 1];
}

LCError LCOneDimLinearFunction::evalDeriv(const std::vector<double> &params, int direction, LCFunctionDeriv **result)
{
	return LCError("not implemented");
}

LCError LCOneDimLinearFunction::evalShapeInfo(const std::vector<double> &params, LCFunctionValue **shapeInfo)
{
	LCError err;
	//LCErrorReturn(validateParameters(params));
	double val = 1;
	for (int i = 0; i < nParams_; i++)
	{
		val += (i + 1)*params[i];
	}

	*shapeInfo = new LCRealFunctionValue(val);

	return err;

}