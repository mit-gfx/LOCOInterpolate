#include <iostream>
#include "LCSample.h"
#include "LCError.h"
#include "LCBasisFunction.h"
#include "LCAdaptiveGridCell.h"
#include "LCMathHelper.h"


void GroupFunctionInfo::add(LCSample* sample, double weight)
{
	if (!initialized)
	{
		centerSum = sample->getCenter() *weight;
		initialized = true;
	}
	else
	{
		centerSum += sample->getCenter() *weight;
	}
	weightSum += weight;
}

LCSample::LCSample(const Eigen::VectorXd center, 
		LCBasisFunction *basisFunction, LCFunctionValue *shapeInfo)
{
	center_ = center;
	addBasisFunction(basisFunction);
	shapeInfo_ = shapeInfo;
}

LCSample::LCSample(const Eigen::VectorXd center, LCFunctionValue *shapeInfo)
{
	center_ = center;
	shapeInfo_ = shapeInfo;
}

void LCSample::removeBasisFucntion(LCBasisFunction* basisFunction)
{

	LCBasisFunctionKey inputKey(basisFunction);
	auto bfMap = basisFunctions_.find(inputKey);
	if (bfMap != basisFunctions_.end())
	{
		basisFunctions_.erase(bfMap);
	}
	else
	{
		std::cout << "could not find to remove " << std::endl;
		basisFunction->log();
		inputKey.log();
		for (auto bf : basisFunctions_)
		{
			std::cout << "---- ";
			bf.first.log();
			std::cout << "-------- ";
			bf.second->log();
		}
		system("pause");
	}
	
}


LCSample::LCSample(const Eigen::VectorXd center, 
	std::vector<LCBasisFunction*> &basisFunctions, LCFunctionValue *shapeInfo)
{
	center_ = center;
	for (auto basisFunction : basisFunctions)
	{
		addBasisFunction(basisFunction);
	}
	shapeInfo_ = shapeInfo;
}


LCError LCSample::getInterpolationWeight(const Eigen::VectorXd pos, double *result)
{
	double combination = 0;
	for (auto basisFunction : basisFunctions_)
	{
		double singleEval = 0; 
		LCErrorReturn(basisFunction.second->eval(pos, &singleEval));
		combination += singleEval;
	}

	*result = combination;

	return LCError();
}


LCError LCSample::getDerivInterpolationWeight(const Eigen::VectorXd pos, int direction, double *result)
{
	double combination = 0;
	for (auto basisFunction : basisFunctions_)
	{
		double singleEval = 0;
		LCErrorReturn(basisFunction.second->evalDeriv(pos, direction, &singleEval));
		combination += singleEval;
	}

	*result = combination;

	return LCError();
}

Eigen::VectorXd & LCSample::getCenter()
{
	return center_;
}

LCFunctionValue * LCSample::getShapeInfo()
{
	return shapeInfo_;
}

void LCSample::setShapeInfo(LCFunctionValue * newShapeInfo)
{
	if (shapeInfo_ != nullptr)
	{
		delete shapeInfo_;
	}
	shapeInfo_ = newShapeInfo;
}


std::unordered_map<LCBasisFunctionKey, LCBasisFunction*> LCSample::getBasisFunctions()
{
	return basisFunctions_;
}


void LCSample::addBasisFunction(LCBasisFunction *basisFunction)
{
	//todo check for replications
	//basisFunction->log();
	LCBasisFunctionKey inputKey(basisFunction);
	auto bfMap = basisFunctions_.find(inputKey);
	if (bfMap == basisFunctions_.end())
	{
		basisFunctions_[inputKey] = basisFunction;
	}
	else
	{
		bfMap->second->addWeight(basisFunction->getWeight());
	}
}


LCError LCSample::refineBasisFunctions(LCAdaptiveGridCell *cell, int dir, double desiredDirSupport,
	LCBasisFunction::LCBasisFunctionType basisType, int *nRefines)
{
	LCError err;
	// first check if this sample needs refinement and the region this sample can affect
	// if it lies on the split region it does not need refinement
	// if not, we check if it's on the right or left side of the split and get the region that this sample can no longer cover
	//std::vector<std::vector<double>> uncoveredRegions;
	//LCErrorReturn(cell->getAffectedRangeWhenSplitting(getCenter(), dir, basisType, &uncoveredRegions));
	//refine all the basis functions inside the cell

	std::vector<LCBasisFunction*> newBasisFunctions;
	std::vector<LCBasisFunctionKey> basisFunctionToRemove;

	for (auto basisFunction : basisFunctions_)
	{
		if (!basisFunction.second->checkRegionOverlap(cell->ranges_))
		{
			continue;
		}
		Eigen::VectorXd currentSupport = basisFunction.second->getSupport();
		if ((currentSupport(dir) - desiredDirSupport) < 0.00001)
		{
			continue; // is already the corrent or lower resolution
		}
		if (std::abs(currentSupport(dir) - 2 * desiredDirSupport) > 0.000001)
		{
			std::cout << "error, the refinement is not half" << std::endl;
			return LCError("error, the refinement is not half");
		}
		(*nRefines)++;
		basisFunctionToRemove.push_back(basisFunction.first);
		basisFunction.second->refineAllDirections(&newBasisFunctions);
		newBasisFunctions.push_back(basisFunction.second);
		//basisFunction->refine(dir, &newBasisFunctions);
	}
	for (auto newBasisFunction : basisFunctionToRemove)
	{
		//basisFunctions_.find(newBasisFunction)->second->log();
		basisFunctions_.erase(newBasisFunction);
	}
	for (auto newBasisFunction : newBasisFunctions)
	{
		//newBasisFunction->log();
		this->addBasisFunction(newBasisFunction);
	}

	// for each of the basis functions that are in the uncoveredRegion
	// remove and add weight to the replacingSamples
	return err;

}


LCError LCSample::extractGroupBasisFunctions(std::unordered_map<LCBasisFunctionKey, GroupFunctionInfo*> * groupedBasisFunctions)
{
	LCError err;

	std::vector<LCBasisFunction*> remBasisFunctions;
	for (auto vis : *groupedBasisFunctions)
	{
		auto bfMap = basisFunctions_.find(vis.first);
		if (bfMap != basisFunctions_.end())
		{
			LCBasisFunction * basisFunc = bfMap->second;
			vis.second->add(this, basisFunc->getWeight());
			remBasisFunctions.push_back(basisFunc);
		}
	}

	for (auto basisFunction : remBasisFunctions)
	{
		this->removeBasisFucntion(basisFunction);
	}
	return err;

}
void LCSample::log()
{
	std::cout << "Precomputed sample:" << std::endl;
	std::cout << "center : " << center_.transpose() << std::endl;
	for (auto basisFunction : basisFunctions_)
	{
		basisFunction.second->log();
	}
}


LCError LCSample::validadeIsCoveredByCells(const std::vector<LCAdaptiveGridCell*> &adjacentCells, 
	const LCFunction* paramShape)
{
	LCError err;
	for (auto basisFunction : basisFunctions_)
	{
		err = basisFunction.second->checkIsCoveredByCells(adjacentCells, paramShape);
		if (!err.isOK())
		{
			std::cout << "error in locality validation" << std::endl;
			basisFunction.second->log();
			return LCError("Locality valiation failed");
		}
	}
	return err;
}


LCError LCSample::getBasisType(LCBasisFunction::LCBasisFunctionType *type)
{
	if (basisFunctions_.size() > 0)
	{
		*type = basisFunctions_.begin()->second->getType();
		return LCError();
	}
	return LCError("sample with no basis functions");
}

