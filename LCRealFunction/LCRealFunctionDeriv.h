#ifndef _OC_REAL_DERIV_INFO_
#define _OC_REAL_DERIV_INFO_
#include <Eigen/Dense>
#include <vector> 
#include "LCFunctionDeriv.h"

class LCFunction;
class LCError;
class LCFunctionValue;
class LCRealFunctionValue;

class LCRealFunctionDeriv : public LCFunctionDeriv
{
public:
	LCRealFunctionDeriv(double val);
	LCRealFunctionDeriv(LCRealFunctionDeriv* other, double weight);
	LCError add(double weight, LCFunctionValue *other);
	void log();

	double val;
};

#endif 
