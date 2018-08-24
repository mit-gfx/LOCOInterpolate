#ifndef _OC_REAL_SHAPE_INFO_
#define _OC_REAL_SHAPE_INFO_
#include <Eigen/Dense>
#include <vector> 
#include <string>
#include "LCFunctionValue.h"


class LCParametricShape;
class LCError;

class LCRealFunctionValue : public LCFunctionValue 
{
public:
	LCRealFunctionValue(double val);
	~LCRealFunctionValue();
	LCError computeDifference(LCFunctionValue *other, double *result);
	LCError computeErrorVector(LCFunctionValue *other, Eigen::VectorXd *result);
	LCError add(double weight, LCFunctionValue *other);
	LCError getHomeophicShapeInfo(LCFunctionValue *source, LCFunctionValue **result);
	LCError addBoundaryConditions(bool useForceAngle);
	void addNewPrecomputedPhysics(std::string name, std::vector<double> &values);
	void clearAllPhysics(){}
	void log();
	double val;
};

#endif 
