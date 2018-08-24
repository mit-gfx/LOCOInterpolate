#ifndef _OC_PRECOMPUTED_SAMPLE_H
#define _OC_PRECOMPUTED_SAMPLE_H
#include <Eigen/Dense>
#include <vector> 
#include <unordered_map>
#include "LCBasisFunction.h"

class LCFunctionValue;
class LCError;
class LCBasisFunction;
class LCAdaptiveGridCell;
class LCFunction;

class GroupFunctionInfo
{
public:
	GroupFunctionInfo(){
		initialized = false;
		weightSum = 0;
	}
	void add(LCSample* sample, double weight);
	Eigen::VectorXd getCenter()
	{
		return centerSum / weightSum;
	}
	double getWeight()
	{
		return weightSum;
	}
private:
	Eigen::VectorXd centerSum;
	double weightSum;
	bool initialized;

};

class LCSample
{
public:
	LCSample(const Eigen::VectorXd center, LCBasisFunction *basisFunction, LCFunctionValue *shapeInfo);
	LCSample(const Eigen::VectorXd center, LCFunctionValue *shapeInfo);
	LCSample(const Eigen::VectorXd center, std::vector<LCBasisFunction*> &basisFunctions, LCFunctionValue *shapeInfo);
	void removeBasisFucntion(LCBasisFunction* basisFunction);

	LCError getInterpolationWeight(const Eigen::VectorXd pos, double *result);
	LCError getDerivInterpolationWeight(const Eigen::VectorXd pos, int direction, double *result);
	Eigen::VectorXd& getCenter();
	LCFunctionValue * getShapeInfo();
	std::unordered_map<LCBasisFunctionKey, LCBasisFunction*> getBasisFunctions();
	LCError refineBasisFunctions(LCAdaptiveGridCell *cell, int dir, double desiredDirSupport,
		LCBasisFunction::LCBasisFunctionType basisType, int *nRefines);
	void addBasisFunction(LCBasisFunction *basisFunction);
	LCError extractGroupBasisFunctions(std::unordered_map<LCBasisFunctionKey, GroupFunctionInfo*> * groupedBasisFunctions);
	void log(); 
	LCError validadeIsCoveredByCells(const std::vector<LCAdaptiveGridCell*> &adjacentCells, const LCFunction* paramShape);

	LCError getBasisType(LCBasisFunction::LCBasisFunctionType *type);
	void setShapeInfo(LCFunctionValue* shapeInfo);

private:
	Eigen::VectorXd center_;
	LCFunctionValue* shapeInfo_; // the shape information stored in the sample
	std::unordered_map<LCBasisFunctionKey, LCBasisFunction*> basisFunctions_;

};

#endif
