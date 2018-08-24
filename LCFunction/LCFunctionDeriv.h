#ifndef _OC_FUNC_DERIV_INFO_
#define _OC_FUNC_DERIV_INFO_
#include <Eigen/Dense>
#include <vector> 

class LCFunction;
class LCError;
class LCFunctionValue;
class LCRealFunctionValue;

class LCFunctionDeriv{
public:
	LCFunctionDeriv();
	static LCError newFromWeightedSum(const std::vector<std::pair<double, LCFunctionValue*>> weightedSamples, LCFunctionDeriv **result);
	virtual LCError add(double weight, LCFunctionValue *other) = 0;
	virtual void log() = 0;

};


#endif 
