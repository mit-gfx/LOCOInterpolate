#ifndef _OC_BASIS_FUNCTION_
#define _OC_BASIS_FUNCTION_

#include <Eigen/Dense>
#include <iostream>

#include <vector> 
#include "LCMathHelper.h"
class LCError;
class LCSample;
class LCFunction;
class LCAdaptiveGridCell;

struct LCBasisFunctionKey;

class LCBasisFunction
{
public:
	enum LCBasisFunctionType{LINEAR_BSPLINE, CUBIC_BSPLINE};
	LCBasisFunction();
	virtual LCBasisFunctionType getType() = 0;
	virtual LCError eval(const Eigen::VectorXd &pos, double *result) = 0;
	virtual LCError evalDeriv(const Eigen::VectorXd &pos, int direction, double *result) = 0;
	virtual LCError refine(int dir, std::vector<LCBasisFunction*> *basisFunctions) = 0;
	virtual void log() = 0;
	Eigen::VectorXd& getSupport();
	Eigen::VectorXd& getCenter();
	double getWeight();
	void setWeightToZero();
	void setWeight(double weight);
	void addWeight(double val);
	bool isEqual(LCBasisFunction* other);
	bool checkRegionOverlap(const std::vector<double> &uncoveredRegion);
	LCError refineAllDirections(std::vector<LCBasisFunction*> *basisFunctions);
	LCError checkIsCoveredByCells(const std::vector<LCAdaptiveGridCell*> &adjacentCells, const LCFunction* paramShape);

	LCError computeReplacingWeights(const std::vector<std::vector<LCSample*>> &replacingSamples,
		int splitDirection, const Eigen::VectorXd &cellSize, int *replacingSamplesIndex, std::vector<double> *replacingWeights);

	static LCError newFromCopyWithWeights(LCBasisFunction* basisFunction, double weight, LCBasisFunction** result);
	static LCError newFromCopyKey(const LCBasisFunctionKey & basisFunctionKey, double weight, LCBasisFunction** result);

protected:
	Eigen::VectorXd center_;
	Eigen::VectorXd support_;
	double weight_;

};

class LCLinearBSpline : public LCBasisFunction
{
public:
	LCLinearBSpline(const Eigen::VectorXd &center, const Eigen::VectorXd &support, double weigth);
	LCBasisFunctionType getType() { return LINEAR_BSPLINE; }
	LCError eval(const Eigen::VectorXd &pos, double *result);
	LCError evalDeriv(const Eigen::VectorXd &pos, int direction, double *result);
	LCError refine(int dir, std::vector<LCBasisFunction*> *basisFunctions);
	void log();

};

class LCCubicBSpline : public LCBasisFunction
{
public:
	LCCubicBSpline(const Eigen::VectorXd &center, const Eigen::VectorXd &support, double weigth);
	LCBasisFunctionType getType() { return CUBIC_BSPLINE; }
	LCError eval(const Eigen::VectorXd &pos, double *result);
	LCError evalDeriv(const Eigen::VectorXd &pos, int direction, double *result);
	LCError refine(int dir, std::vector<LCBasisFunction*> *basisFunctions);
	void log();

};


struct LCBasisFunctionKey
{
	LCBasisFunctionKey(LCBasisFunction * bf)
	{
		center = bf->getCenter();
		support = bf->getSupport();
		type = bf->getType();
	}
	Eigen::VectorXd center;
	Eigen::VectorXd support;
	LCBasisFunction::LCBasisFunctionType type;

	bool operator==(const LCBasisFunctionKey &other) const
	{
		return (((center - other.center).norm() < LCMathHelper::EPSILON)
			&& ((support - other.support).norm() < LCMathHelper::EPSILON));
	}
	void log() const
	{
		std::cout << "center = " << center.transpose() << " support = " << support << std::endl;
	}
	// bool operator<(const edge &e) { // less operator
	bool operator<(const LCBasisFunctionKey &bf) const { // const-correct less operator
		return center.norm() < bf.center.norm();
	}
};

namespace std {
	template <>
	struct hash<LCBasisFunctionKey>
	{
		std::size_t operator()(const LCBasisFunctionKey& k) const
		{
			std::size_t hashValue = 0;

			for (int i = 0; i < k.center.size(); i++)
			{
				//hashValue ^= std::hash<double>()(round(100.0 * k.center(i))) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
			}
			//hashValue ^= std::hash<double>()(k.support(0)) + 0x9e3779b9 + (hashValue << 6) + (hashValue >> 2);
			return hashValue;

		}
	};
}

#endif 
