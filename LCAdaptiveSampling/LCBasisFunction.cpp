#include <iostream>
#include "LCBasisFunction.h"
#include "LCSample.h"
#include "LCError.h"
#include "LCMathHelper.h"
#include "LCAdaptiveGridCell.h"
#include "LCFunction.h"

LCBasisFunction::LCBasisFunction()
{

}


bool LCBasisFunction::isEqual(LCBasisFunction* other)
{
	if (getType() != other->getType())
	{
		return false;
	}
	if ((getCenter() - other->getCenter()).norm() > LCMathHelper::EPSILON)
	{
		return false;
	}
	if ((getSupport() - other->getSupport()).norm() > LCMathHelper::EPSILON)
	{
		return false;
	}
	return true;
}

Eigen::VectorXd& LCBasisFunction::getSupport()
{
	return support_;
}

Eigen::VectorXd& LCBasisFunction::getCenter()
{
	return center_;
}

double LCBasisFunction::getWeight()
{
	return weight_;
}
void LCBasisFunction::setWeightToZero()
{
	weight_ = 0;
}

void LCBasisFunction::setWeight(double weight)
{
	weight_ = weight;
}

void LCBasisFunction::addWeight(double val)
{
	weight_ += val;
}


bool LCBasisFunction::checkRegionOverlap(const std::vector<double> &uncoveredRegion)
{

	int nParams = center_.size();
	for (int i = 0; i < nParams; i++)
	{
		double minBasisFunction = center_[i] - support_[i];
		double maxBasisFunction = center_[i] + support_[i];
		double minRegion = uncoveredRegion[2 * i];
		double maxRegion = uncoveredRegion[2 * i + 1];
		if ((maxRegion - minRegion) < LCMathHelper::EPSILON)
		{
			//return false;
		}
		double allowedError = support_[i] / 100.0;
		if ((maxBasisFunction < (minRegion + allowedError)) || (minBasisFunction >(maxRegion - allowedError)))
		{
			return false;
		}
	}
	return true;
}


LCError LCBasisFunction::computeReplacingWeights(const std::vector<std::vector<LCSample*>> &replacingSamplesList,
	int splitDirection, const Eigen::VectorXd &cellSize, int *replacingSamplesIndex, std::vector<double> *replacingWeights)
{
	int projectedSize = center_.size() - 1;
	if (projectedSize == 0)
	{
		int index = 0;
		for (auto replacingSamples : replacingSamplesList)
		{
			if (replacingSamples.size() != 1)
			{
				return LCError("one dim and the number of replacing samples is not one");
			}
			if ((center_ - replacingSamples[0]->getCenter()).norm() < LCMathHelper::EPSILON)
			{
				replacingWeights->push_back(1);
				*replacingSamplesIndex = index;
				return LCError();
			}
			index++;
		}
		return LCError("the center of the basis function is not the replacing sample on one dim");
	}
	int index = 0;
	for (auto replacingSamples : replacingSamplesList)
	{
		if ((int)replacingSamples.size() != pow(2, projectedSize))
		{
			return LCError("the number or replacing samples is not 2^(K-1)");
		}
		if (abs(replacingSamples[0]->getCenter()[splitDirection] - center_[splitDirection]) < LCMathHelper::EPSILON)
		{
			//the center of the basis function on the split direction needs to be the smale
			Eigen::VectorXd projectedMidVal = LCMathHelper::removeDimension(replacingSamples[0]->getCenter(), splitDirection);
			Eigen::VectorXd projectedCellSize = LCMathHelper::removeDimension(cellSize, splitDirection);
			Eigen::VectorXd projectedPos = LCMathHelper::removeDimension(center_, splitDirection);
			Eigen::VectorXd alpha = (projectedPos - projectedMidVal).cwiseQuotient(projectedCellSize);

			std::vector<Eigen::VectorXi> combinations;
			LCMathHelper::computeCombinations(projectedSize, 2, &combinations);
			double weightSum = 0;
			for (int i = 0; i < replacingSamples.size(); i++)
			{
				if (abs(replacingSamples[i]->getCenter()[splitDirection] - center_[splitDirection]) > LCMathHelper::EPSILON)
				{
					std::cout << "basis function center : " << center_.transpose() << std::endl;
					std::cout << "basis function support : " << support_.transpose() << std::endl;
					std::cout << "replacing sample : " << replacingSamples[i]->getCenter().transpose() << std::endl;
					return LCError("the center of the basis function is not the replacing sample restricted to one dim");
				}
				double midWeight = 1;
				for (int j = 0; j < projectedSize; j++)
				{
					if (combinations[i](j) == 1)
					{
						midWeight *= alpha[j];
					}
					else
					{
						midWeight *= (1 - alpha[j]);
					}
				}
				weightSum += midWeight;
				replacingWeights->push_back(midWeight);
			}
			if (abs(weightSum - 1) > LCMathHelper::EPSILON)
			{
				return LCError("sum of weights not zero");
			}
			*replacingSamplesIndex = index;
			return LCError();
		}
		index++;
	}
	std::cout << "basis function center : " << center_.transpose() << std::endl;
	std::cout << "basis function support : " << support_.transpose() << std::endl;
	std::cout << "replacing sampleList 0 : " << replacingSamplesList[0][0]->getCenter().transpose() << std::endl;
	std::cout << "replacing sampleList 1 : " << replacingSamplesList[1][0]->getCenter().transpose() << std::endl;
	std::cout << "replacing sampleList 2 : " << replacingSamplesList[2][0]->getCenter().transpose() << std::endl;
	return LCError("the center of the basis in not on the replacing sample");
}


LCError LCBasisFunction::newFromCopyWithWeights(LCBasisFunction* basisFunction, double weight, LCBasisFunction** result)
{
	LCError err;

	if (dynamic_cast<LCLinearBSpline*>(basisFunction))
	{
		*result = new LCLinearBSpline(basisFunction->getCenter(), basisFunction->getSupport(), basisFunction->getWeight() * weight);
		return err;
	}

	if (dynamic_cast<LCCubicBSpline*>(basisFunction))
	{
		*result = new LCCubicBSpline(basisFunction->getCenter(), basisFunction->getSupport(), basisFunction->getWeight() * weight);
		return err;
	}
	return LCError("Basis function type not specified");

}

LCError LCBasisFunction::newFromCopyKey(const LCBasisFunctionKey &basisFunction, double weight, LCBasisFunction** result)
{
	LCError err;

	Eigen::VectorXd center = basisFunction.center;
	Eigen::VectorXd support = basisFunction.support;
	if (basisFunction.type == LCBasisFunction::LCBasisFunctionType::LINEAR_BSPLINE)
	{
		*result = new LCLinearBSpline(center, support, weight);
		return err;
	}

	if (basisFunction.type == LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE)
	{
		*result = new LCCubicBSpline(center, support, weight);
		return err;
	}
	return LCError("Basis function type not specified");

}

LCError LCBasisFunction::refineAllDirections(std::vector<LCBasisFunction*> *basisFunctions)
{
	//std::cout << "original basis function" << std::endl;
	//this->log();
	LCError err;
	std::vector<LCBasisFunction*> allNewBasisFunctions;
	for (int dir = 0; dir < (this->getCenter()).size(); dir++)
	{
		std::vector<LCBasisFunction*> newBasisFunctions;
		refine(dir, &newBasisFunctions);
		for (LCBasisFunction* basisFunction : allNewBasisFunctions)
		{
			basisFunction->refine(dir, &newBasisFunctions);
		}
		allNewBasisFunctions.insert(allNewBasisFunctions.end(), newBasisFunctions.begin(), newBasisFunctions.end());
	}

	//std::cout << "final basis funcitons: " << allNewBasisFunctions.size() + 1 << std::endl;
	//this->log();
	//for (LCBasisFunction* basisFunction : allNewBasisFunctions)
	//{
	//	basisFunction->log();
	//}


	basisFunctions->insert(basisFunctions->end(), allNewBasisFunctions.begin(), allNewBasisFunctions.end());


	return err;
}


LCError LCBasisFunction::checkIsCoveredByCells(const std::vector<LCAdaptiveGridCell*> &adjacentCells, const LCFunction* paramShape)

{
	LCError err;

	Eigen::VectorXd center = getCenter();
	Eigen::VectorXd support = getSupport();
	std::vector<Eigen::VectorXi> combinations;
	LCMathHelper::computeCombinations(center.size(), 2, &combinations);
	for (auto combination : combinations)
	{
		Eigen::VectorXd multi = Eigen::VectorXd::Ones(center.size()) - 2 * combination.cast<double>();
		Eigen::VectorXd endCorner = center + multi.cwiseProduct(support);
		std::vector<double> centerVec;
		LCMathHelper::eigen2StdVector(center, &centerVec);
		std::vector<double> endCornerVec;
		bool combinationInCells = false;
		LCMathHelper::eigen2StdVector(endCorner, &endCornerVec);
		for (auto adjacentCell : adjacentCells)
		{
			if (adjacentCell->cointainsPoint(endCorner) || !paramShape->validateParameters(endCornerVec).isOK()
				|| !paramShape->validateParameters(centerVec).isOK())
			{
				combinationInCells = true;
			}
		}
		if (!combinationInCells)
		{
			return LCError("samples is outside neighbors");
		}
	}
	return err;
}
//----------------------------------------------------------------------------------------------
LCLinearBSpline::LCLinearBSpline(const Eigen::VectorXd &center,
	const Eigen::VectorXd &support, double weight) : LCBasisFunction()
{
	center_ = center;
	support_ = support;
	weight_ = weight;
}

LCError LCLinearBSpline::eval(const Eigen::VectorXd &pos, double *result)
{
	if (pos.size() != center_.size())
	{
		return LCError("Input position of wrong size");
	}
	double value = 1;
	for (int i = 0; i < center_.size(); i++)
	{
		double dimVal = 0;
		double d = abs(pos(i) - center_(i)) / support_(i);
		if (d <= 1)
		{
			dimVal = 1 - d;
		}
		value *= dimVal;
	}
	*result = value*weight_;
	return LCError();

}

LCError LCLinearBSpline::evalDeriv(const Eigen::VectorXd &pos, int direction, double *result)
{
	if (pos.size() != center_.size())
	{
		return LCError("Input position of wrong size");
	}
	double value = 1;
	for (int i = 0; i < center_.size(); i++)
	{
		double dimVal = 0;
		double d = (pos(i) - center_(i)) / support_(i);
		if (i == direction)
		{
			if (d <= 1)
			{
				if (d < 0)
				{
					dimVal = 1.0 / support_(i);
				}
				if (d > 0)
				{
					dimVal = -1.0 / support_(i);
				}
			}
		}
		else
		{
			if (d <= 1)
			{
				dimVal = 1 - d;
			}
		}
		value *= dimVal;
	}
	*result = value*weight_;
	return LCError();

}

LCError LCLinearBSpline::refine(int dir, std::vector<LCBasisFunction*> *basisFunctions)
{
	LCError err;
	Eigen::VectorXd newSupport = support_;
	newSupport(dir) = 0.5 *support_(dir);
	Eigen::VectorXd newCenterA = center_;
	newCenterA(dir) -= 0.5 *support_(dir);
	Eigen::VectorXd newCenterB = center_;
	newCenterB(dir) += 0.5 *support_(dir);
	double newValA, newValB;
	LCErrorReturn(eval(newCenterA, &newValA));
	LCErrorReturn(eval(newCenterB, &newValB));
	if (newValA != weight_*0.5)
	{
		return LCError("wrong weight eval");
	}
	if (newValB != weight_*0.5)
	{
		return LCError("wrong weight eval");
	}
	basisFunctions->push_back(new LCLinearBSpline(newCenterA, newSupport, newValA));
	basisFunctions->push_back(new LCLinearBSpline(newCenterB, newSupport, newValB));
	support_ = newSupport;
	return err;
}


void LCLinearBSpline::log()
{
	std::cout << "Linear basis center: " << center_.transpose();
	std::cout << "    support: " << support_.transpose();
	std::cout << "     weight: " << weight_ << std::endl;
}




//----------------------------------------------------------------------------------------------
LCCubicBSpline::LCCubicBSpline(const Eigen::VectorXd &center,
	const Eigen::VectorXd &support, double weight) : LCBasisFunction()
{
	if (center(0) > 0 && center(0) < 0.005 || center(1) > 0 && center(1) < 0.005)
	{
		std::cout << "small weight constructor " << center.transpose() << std::endl;
		system("pause");
	}

	center_ = center;
	support_ = support;
	weight_ = weight;
}

LCError LCCubicBSpline::eval(const Eigen::VectorXd &pos, double *result)
{
	if (pos.size() != center_.size())
	{
		return LCError("Input position of wrong size");
	}
	double value = 1;
	for (int i = 0; i < center_.size(); i++)
	{
		double dimVal = 0;
		double d = abs(pos(i) - center_(i)) / support_(i);
		if (d <= 1)
		{
			double p = 2 * ((pos(i) - center_(i)) / support_(i) + 1);
			int k = floor(p);
			double t = 1 - (p - k);
			if (k == 0)
			{
				dimVal = (1 - t) * (1 - t) * (1 - t) / 6;
			}
			if (k == 1)
			{
				dimVal = (3 * t * t * t - 6 * t * t + 4) / 6;
			}
			if (k == 2)
			{
				dimVal = (-3 * t * t* t + 3 * t * t + 3 * t + 1) / 6;
			}

			if (k == 3)
			{
				dimVal = t * t* t / 6;
			}
		}
		value *= dimVal;
	}
	*result = value*weight_;
	return LCError();

}

LCError LCCubicBSpline::evalDeriv(const Eigen::VectorXd &pos, int direction, double *result)
{
	if (pos.size() != center_.size())
	{
		return LCError("Input position of wrong size");
	}
	double value = 1;
	for (int i = 0; i < center_.size(); i++)
	{
		double dimVal = 0;
		double d = abs(pos(i) - center_(i)) / support_(i);
		if (d <= 1)
		{
			double p = 2 * ((pos(i) - center_(i)) / support_(i) + 1);
			int k = floor(p);
			double t = 1 - (p - k);
			if (i == direction)
			{
				double dt = -2 / support_(i);
				if (k == 0)
				{
					dimVal = -(1 - t)*(1 - t)*dt / 2;
				}
				if (k == 1)
				{
					dimVal = (3 * t*t - 4 * t)*dt / 2;
				}
				if (k == 2)
				{
					dimVal = (-3 * t*t + 2 * t + 1)*dt / 2;
				}
				if (k == 3)
				{
					dimVal = t*t*dt / 2;
				}
			}
			else
			{
				if (k == 0)
				{
					dimVal = (1 - t) * (1 - t) * (1 - t) / 6;
				}
				if (k == 1)
				{
					dimVal = (3 * t * t * t - 6 * t * t + 4) / 6;
				}
				if (k == 2)
				{
					dimVal = (-3 * t * t* t + 3 * t * t + 3 * t + 1) / 6;
				}

				if (k == 3)
				{
					dimVal = t * t* t / 6;
				}
			}
		}
		value *= dimVal;
	}
	*result = value*weight_;
	return LCError();

}

LCError LCCubicBSpline::refine(int dir, std::vector<LCBasisFunction*> *basisFunctions)
{

	//TODOCUBIC
	LCError err;
	Eigen::VectorXd newSupport = support_;
	newSupport(dir) = 0.5 *support_(dir);

	Eigen::VectorXd newCenterAM = center_;
	newCenterAM(dir) -= 0.25 *support_(dir);
	Eigen::VectorXd newCenterBM = center_;
	newCenterBM(dir) -= 0.5 *support_(dir);
	Eigen::VectorXd newCenterAP = center_;
	newCenterAP(dir) += 0.25 *support_(dir);
	Eigen::VectorXd newCenterBP = center_;
	newCenterBP(dir) += 0.5 *support_(dir);
	basisFunctions->push_back(new LCCubicBSpline(newCenterAM, newSupport, weight_ * 1.0 / 2.0));
	basisFunctions->push_back(new LCCubicBSpline(newCenterBM, newSupport, weight_ * 1.0 / 8.0));
	basisFunctions->push_back(new LCCubicBSpline(newCenterAP, newSupport, weight_ * 1.0 / 2.0));
	basisFunctions->push_back(new LCCubicBSpline(newCenterBP, newSupport, weight_ * 1.0 / 8.0));
	support_ = newSupport;
	weight_ *= 6.0 / 8.0;
	if (center_(0) > 0 && center_(0) < 0.005 || center_(1) > 0 && center_(1) < 0.005)
	{
		std::cout << "small weight refine " << center_.transpose() << std::endl;
		system("pause");
	}
	return err;
}





void LCCubicBSpline::log()
{
	std::cout << "CUBIC basis center: " << center_.transpose();
	std::cout << "    support: " << support_.transpose();
	std::cout << "     weight: " << weight_ << std::endl;
}
