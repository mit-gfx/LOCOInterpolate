#ifndef _OC_FUNC_SHAPE_INFO_
#define _OC_FUNC_SHAPE_INFO_
#include <Eigen/Dense>
#include <vector> 
#include <string>


class LCParametricShape;
class LCError;

class LCFunctionValue{
public:
	LCFunctionValue();
	virtual ~LCFunctionValue() {};
	static LCError newFromWeightedSum(const std::vector<std::pair<double, LCFunctionValue*>> weightedSamples, LCFunctionValue **result);
	virtual LCError computeDifference(LCFunctionValue *other, double *result) = 0;
	virtual LCError computeErrorVector(LCFunctionValue *other, Eigen::VectorXd *result) = 0;
	virtual LCError add(double weight, LCFunctionValue *other) = 0;
	virtual LCError getHomeophicShapeInfo(LCFunctionValue *source, LCFunctionValue **result) = 0;
	virtual LCError addBoundaryConditions(bool useForceAngle) = 0;
	virtual void addNewPrecomputedPhysics(std::string name, std::vector<double> &values) = 0;
	virtual void clearAllPhysics() = 0;
	virtual void log() = 0;
};


#endif 
