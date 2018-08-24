#ifndef _OC_HARD_CODED_LINEAR_SHAPE_
#define _OC_HARD_CODED_LINEAR_SHAPE_

#include "LCFunction.h"

/* This is used for testing the interpolation scheme. It a simple function from
R^N -> R where N are the number of params. If you set linear = true it will be a linear
function which can be used to test linear precision and if not it will a sum on
sines and cosines. */


class LCOneDimLinearFunction : public LCFunction {
public:
	LCOneDimLinearFunction(int nParams);
	int getNParams() const;
	double getMinRange(int iParam) const;
	double getMaxRange(int iParam) const;
	LCError evalDeriv(const std::vector<double> &params, int direction, LCFunctionDeriv **result);
	LCError evalShapeInfo(const std::vector<double> &params, LCFunctionValue **result);

private:
	int nParams_;
	std::vector<double> ranges_;
};


#endif