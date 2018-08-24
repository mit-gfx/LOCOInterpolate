#include <iostream>
#include "LCAdaptiveGrid.h"
#include "LCFunction.h"
#include "LCError.h"
#include "LCMathHelper.h"
#include "LCSample.h"
#include "LCBasisFunction.h"
#include "LCAdaptiveGridCell.h"
#include "LCSampledFunction.h"
#include "LCFunctionValue.h"
#include <Eigen/Dense>
#include <queue>
#include <set>

LCAdaptiveGrid::LCAdaptiveGrid(LCFunction *parametricShape, LCBasisFunction::LCBasisFunctionType basisType,
	LCAdaptiveSamplingParams* params, bool evalCorners)

{
	std::cout << "first params constructor " << std::endl;
	params->log();
	params_ = params;
	parametricShape_ = parametricShape;
	std::cout << "threshold " << params_->maxTreeDepth_ << std::endl;
	allowedError_ = params_->threshold_ * Eigen::VectorXd::Ones(1);
	basisType_ = basisType;
	maxLevel_ = params->maxTreeDepth_;
	evalCorners_ = evalCorners;

	// bootstrap step: uniform sampling 
	std::cout << "doing bootstrap part" << std::endl;

	LCError err = bootstrapSample(params_->bootstrapGridSize_);
	if (!err.isOK())
	{
		std::cout << err.internalDescription() << std::endl;
	}

	std::cout << "doing adaptive part" << std::endl;
	//remember to uncomment this:
	err = adaptiveSample();
	if (!err.isOK())
	{
		std::cout << err.internalDescription() << std::endl;
	}
}

LCAdaptiveGrid::LCAdaptiveGrid(LCFunction *parametricShape, LCBasisFunction::LCBasisFunctionType basisType,
	Eigen::VectorXd allowedError, int bootstrapGridSize, int maxLevel, bool evalCorners)

{
	std::cout << "first constructor " << std::endl;
	parametricShape_ = parametricShape;
	allowedError_ = allowedError;
	basisType_ = basisType;
	maxLevel_ = maxLevel;
	evalCorners_ = evalCorners;
	LCAdaptiveSamplingParams* params = new LCAdaptiveSamplingParams();
	params_ = params;
	params_->log();

	// bootstrap step: uniform sampling 

	LCError err = bootstrapSample(bootstrapGridSize);
	if (!err.isOK())
	{
		std::cout << err.internalDescription() << std::endl;
	}

	//std::cout << "doing adaptive part" << std::endl;
	//remember to uncomment this:
	err = adaptiveSample();
	if (!err.isOK())
	{
		std::cout << err.internalDescription() << std::endl;
	}
}
LCAdaptiveGrid::LCAdaptiveGrid(LCFunction *parametricShape, LCSampledFunction *coarseShape,
	Eigen::VectorXd allowedError, int maxLevel, bool evalCorners)
{
	root_ = coarseShape->getRoot();
	samples_ = coarseShape->getSamples();
	borderCells_ = coarseShape->getBorderCells();
	parametricShape_ = parametricShape;
	allowedError_ = allowedError;
	evalCorners_ = evalCorners;
	for (auto sample : samples_)
	{
		if ((sample->getBasisType(&basisType_)).isOK())
		{
			break;
		}
	}
	maxLevel_ = maxLevel + parametricShape_->getNParams();
	// bootstrap step: uniform sampling 
	//std::cout << "doing adaptive part" << std::endl;

	LCError err = adaptiveSample();
	if (!err.isOK())
	{
		std::cout << err.internalDescription() << std::endl;
	}
}


LCError LCAdaptiveGrid::evalShapeInfo(const std::vector<double> &shapeParams, LCFunctionValue **result)
{
	auto params = parametricShape_->mapToStandardHypercube(shapeParams);
	LCAdaptiveGridCell *containingCell = nullptr;
	LCErrorReturn(getContainingCell(params, &containingCell));
	return containingCell->evalShapeInfo(params, basisType_, result);
}

LCError LCAdaptiveGrid::getContainingCell(const std::vector<double> &params, LCAdaptiveGridCell **result)
{
	LCError err = root_->getCointainingLeaf(params, result);
	if (err.isOK())
	{
		return err;
	}

	for (auto borderCell : borderCells_)
	{
		if (borderCell->evalPametersAreInCell(params).isOK())
		{
			*result = borderCell;
			return LCError();
		}
	}
	return LCError("could not find containing cell");
}

LCError LCAdaptiveGrid::getSmallestCointainingLeaf(const std::vector<double> &params, LCAdaptiveGridCell **result)
{

	LCError err = root_->getSmallestCointainingLeaf(params, result);
	if (err.isOK())
	{
		return err;
	}

	for (auto borderCell : borderCells_)
	{
		if (borderCell->evalPametersAreInCell(params).isOK())
		{
			*result = borderCell;
			return LCError();
		}
	}
	return LCError("could not find containing cell");
}



LCError LCAdaptiveGrid::getAdjacentCells(const std::vector<double> &params, std::unordered_set<LCAdaptiveGridCell*> *result)
{

	LCAdaptiveGridCell* aux;

	/*std::cout << "basisCenter: ";
	for (int i = 0; i < params.size(); i++)
	{
	std::cout << params[i] << " ";
	};
	std::cout << std::endl;*/

	bool interestingBasis = false;
	/*if ((std::abs(params[0] - 0.937509) < 0.000001) && (std::abs(params[1] - 0.875003) < 0.000001))
	{
	std::cout << "at the interesting basis center " << std::endl;
	std::cout << "basisCenter: ";
	interestingBasis = true;
	for (int i = 0; i < params.size(); i++)
	{
	std::cout << params[i] << " ";
	};
	std::cout << std::endl;
	}*/

	LCError outsideSupport = root_->getAdjacentLeaves(params, result);
	if (interestingBasis)
	{
		std::cout << "interesting basis center outsideSupport.isOK " << outsideSupport.isOK() << std::endl;
	}

	bool isOnBorder = false;
	//std::cout << "border cell size " << borderCells_.size() << std::endl;
	for (auto borderCell : borderCells_)
	{
		if (borderCell->evalPametersAreInCell(params).isOK())
		{
			result->insert(borderCell);
			isOnBorder = true;
		}
	}

	if (interestingBasis)
	{
		std::cout << "interesting basis center on the border " << isOnBorder << std::endl;
	}

	if (!isOnBorder && !outsideSupport.isOK())
	{
		system("pause");
		LCAdaptiveGridCell *leafCell;
		LCError err2 = root_->getCointainingLeaf(params, &leafCell);
		std::cout << "not on border and not outsideSupport basisCenter: ";
		for (int i = 0; i < params.size(); i++)
		{
			std::cout << params[i] << " ";
		}
		LCError outsideSupport2 = root_->getAdjacentLeaves(params, result);
		return LCError("Not found on border");
	}

	return LCError();
}

LCError LCAdaptiveGrid::getDoubleAdjacentCells(const std::vector<double> &params, std::unordered_set<LCAdaptiveGridCell*> *result)
{
	LCError err;


	std::unordered_set<LCAdaptiveGridCell*> imediateAdjacenteLeaves;
	err = getAdjacentCells(params, &imediateAdjacenteLeaves);
	if (!err.isOK())
	{
		std::cout << err.internalDescription() << std::endl;
		return err;
	}

	result->insert(imediateAdjacenteLeaves.begin(), imediateAdjacenteLeaves.end());


	for (auto adjCell : imediateAdjacenteLeaves)
	{
		//std::cout << "immediate adj: ";  adjCell->logRanges(0);
		for (auto doubleNeighbor : adjCell->getNeighbors())
		{
		//	std::cout << "double adj: ";  doubleNeighbor->logRanges(0);
			result->insert(doubleNeighbor);
		}
	}

	//result->insert(result->end(), doubleAdjSet.begin(), doubleAdjSet.end());
	return err;

}

LCAdaptiveGridCell * LCAdaptiveGrid::getRoot()
{
	return root_;
}

std::vector<LCSample*> LCAdaptiveGrid::getSamples()
{
	return samples_;
}

void LCAdaptiveGrid::setLeafInformation(LCAdaptiveGridCell* leafCell, const std::vector<LCSample*> &samples)
{
	std::vector<LCSample*> leafSamples;
	for (auto sample : samples)
	{
		if (leafCell->cointainsPoint(sample->getCenter()))
		{
			leafSamples.push_back(sample);
		}
	}
	LCFunctionValue *centerShapeInfo = nullptr;
	std::vector<double> center;
	leafCell->getCenter(&center);
	parametricShape_->evalShapeInfo(parametricShape_->mapFromStandardHypercube(center), &centerShapeInfo);
	leafCell->createHomeomorphicMap(centerShapeInfo, leafSamples);

}


LCError LCAdaptiveGrid::bootstrapUniformSplit(LCAdaptiveGridCell* cell, int  direction, int levels)
{
	LCError err;
	if (levels > 0)
	{
		int nextDirection = ((direction + 1) % parametricShape_->getNParams());
		if (nextDirection == 0)
		{
			levels--;
		}
		LCErrorReturn(cell->split(direction));
		LCErrorReturn(bootstrapUniformSplit(cell->getChild1(), nextDirection, levels));
		LCErrorReturn(bootstrapUniformSplit(cell->getChild2(), nextDirection, levels));
	}
	return err;
}
std::vector<LCAdaptiveGridCell*> LCAdaptiveGrid::getBorderCells()
{
	return borderCells_;
}


LCError LCAdaptiveGrid::bootstrapSample(int nLevels)
{
	LCError err;
	if (nLevels < 0)
	{
		return LCError("number of levels must be positive");
	}
	int nSamples = pow(2, nLevels);
	std::vector<Eigen::VectorXi> sampleCombinations;
	int nParams = parametricShape_->getNParams();
	LCErrorReturn(LCMathHelper::computeCombinations(nParams, nSamples + 1, &sampleCombinations));
	Eigen::VectorXd cellSize = (2 / (double)nSamples) * Eigen::VectorXd::Ones(nParams);
	Eigen::VectorXd shapeMinRange = (-1) * Eigen::VectorXd::Ones(nParams);
	Eigen::VectorXd shapeMaxRange = Eigen::VectorXd::Ones(nParams);

	for (auto sampleCombination : sampleCombinations)
	{
		Eigen::VectorXd alphaVec = sampleCombination.cast<double>() / (double)nSamples;
		Eigen::VectorXd paramsVec = shapeMinRange.cwiseProduct(Eigen::VectorXd::Ones(nParams) - alphaVec) + shapeMaxRange.cwiseProduct(alphaVec);
		std::vector<double> params;
		LCMathHelper::eigen2StdVector(paramsVec, &params);
		LCFunctionValue* shapeInfo = nullptr;
		LCErrorReturn(parametricShape_->evalShapeInfo(parametricShape_->mapFromStandardHypercube(params), &shapeInfo));
		switch (basisType_)
		{
		case LCBasisFunction::LCBasisFunctionType::LINEAR_BSPLINE:
		{
			LCBasisFunction *basisFunction = new LCLinearBSpline(paramsVec, cellSize, 1.0);
			LCSample* precomputedSample = new LCSample(paramsVec, basisFunction, shapeInfo);
			samples_.push_back(precomputedSample);
		}
		break;
		case LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE:
		{
			LCBasisFunction * basisFunction = new LCCubicBSpline(paramsVec, 2 * cellSize, 1.0);
			LCSample* precomputedSample = new LCSample(paramsVec, basisFunction, shapeInfo);
			samples_.push_back(precomputedSample);
		}
		break;
		default:
			return LCError("basis type unknown");
			break;
		}
	}

	std::vector<double> ranges(2 * nParams);
	for (int i = 0; i < nParams; i++)
	{
		ranges[2 * i] = -1;
		ranges[2 * i + 1] = 1;
	}

	root_ = new LCAdaptiveGridCell(0, nullptr, ranges);

	LCErrorReturn(bootstrapUniformSplit(root_, 0, nLevels));

	std::vector<LCAdaptiveGridCell*> leafCells;
	root_->getAllLeafCells(&leafCells);
	for (int i = 0; i < leafCells.size(); i++)
	{
		setLeafInformation(leafCells[i], samples_);
	}


	if (basisType_ == LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE)
	{

		//create border samples
		std::vector<Eigen::VectorXi> cornerCombiantions;
		LCErrorReturn(LCMathHelper::computeCornerCombinations(nParams, nSamples + 3, &cornerCombiantions));
		for (auto corner : cornerCombiantions)
		{
			Eigen::VectorXd alphaVec = (corner.cast<double>() - Eigen::VectorXd::Ones(nParams)) / (double)nSamples;
			Eigen::VectorXd sampleCenter = shapeMinRange.cwiseProduct(Eigen::VectorXd::Ones(nParams) - alphaVec) + shapeMaxRange.cwiseProduct(alphaVec);
			//std::cout << "corner sample center = " << sampleCenter.transpose() << std::endl;
			//std::cout << "corner alphaVec = " << alphaVec.transpose() << std::endl;
			Eigen::VectorXd projectedCorner = corner.cast<double>() - Eigen::VectorXd::Ones(nParams);
			for (int i = 0; i < nParams; i++)
			{
				if (projectedCorner(i) < 0)
				{
					projectedCorner(i) = 0;
				}
				if (projectedCorner(i) > nSamples)
				{
					projectedCorner(i) = nSamples;
				}
			}
			alphaVec = projectedCorner.cast<double>() / (double)nSamples;
			Eigen::VectorXd sampleRefereceCenter = shapeMinRange.cwiseProduct(Eigen::VectorXd::Ones(nParams) - alphaVec) + shapeMaxRange.cwiseProduct(alphaVec);
			LCBasisFunction * basisFunction = new LCCubicBSpline(sampleCenter, 2 * cellSize, 1.0);
			LCSample* sampleReference;
			LCErrorReturn(this->getSample(sampleRefereceCenter, &sampleReference));
			LCFunctionValue* refereceShapeInfo = sampleReference->getShapeInfo();
			//std::cout << "this is evalCorners_ " << evalCorners_ << std::endl;
			if (evalCorners_)
			{
				//std::cout << "this is the reference shape info before the if block " << refereceShapeInfo << std::endl;
				std::vector<double> sampleCenterVec;
				LCMathHelper::eigen2StdVector(sampleCenter, &sampleCenterVec);

				LCErrorReturn(parametricShape_->evalShapeInfo(parametricShape_->mapFromStandardHypercube(sampleCenterVec), &refereceShapeInfo));

				//std::cout << "this is the reference shape info after the if block " << refereceShapeInfo << std::endl;
			}
			LCSample* precomputedSample = new LCSample(sampleCenter, basisFunction, refereceShapeInfo);
			samples_.push_back(precomputedSample);
		}

		//create border cells
		std::vector<Eigen::VectorXi> cornerCellCombiantions;
		LCErrorReturn(LCMathHelper::computeCornerCombinations(nParams, nSamples + 2, &cornerCellCombiantions));
		std::unordered_map<int, LCAdaptiveGridCell*> cornerCellsFromCombination;
		for (auto corner : cornerCellCombiantions)
		{
			Eigen::VectorXd alphaVec = (corner.cast<double>() - Eigen::VectorXd::Ones(nParams)) / (double)nSamples;
			Eigen::VectorXd cellMin = shapeMinRange.cwiseProduct(Eigen::VectorXd::Ones(nParams) - alphaVec) + shapeMaxRange.cwiseProduct(alphaVec);
			Eigen::VectorXd cellMax = cellMin + (shapeMaxRange - shapeMinRange) / (double)nSamples;
			std::vector<double> cellRanges(2 * nParams);
			for (int i = 0; i < nParams; i++)
			{
				cellRanges[2 * i] = cellMin(i);
				cellRanges[2 * i + 1] = cellMax(i);
			}
			//std::cout << "cellMin = " << cellMin.transpose() << " cell max = " << cellMax.transpose() << std::endl;
			LCAdaptiveGridCell * borderCell = new LCAdaptiveGridCell(-1, nullptr, cellRanges);
			//std::cout << "new corner cell : "; borderCell->logRanges(0);
			Eigen::VectorXd cellCenter = 0.5*(cellMin + cellMax);
			std::vector<double> cellCenterVec;
			LCMathHelper::eigen2StdVector(cellCenter, &cellCenterVec);
			//add samples and neighbours
			std::unordered_set<LCAdaptiveGridCell*> neighbours;
			std::vector<LCSample*> cellSamples;
			for (auto sample : samples_)
			{
				if (borderCell->cointainsPoint(sample->getCenter()))
				{
					cellSamples.push_back(sample);
					std::vector<double> sampleCenter;
					LCMathHelper::eigen2StdVector(sample->getCenter(), &sampleCenter);
					root_->getAdjacentLeaves(sampleCenter, &neighbours);
					//std::cout << sample->getCenter() << std::endl;
				}
			}
			std::vector<double> borderCellCenter;
			borderCell->getCenter(&borderCellCenter);
			for (int i = 0; i < nParams; i++)
			{
				borderCellCenter[i] = std::max(ranges[2 * i], borderCellCenter[i]);
				borderCellCenter[i] = std::min(ranges[2 * i + 1], borderCellCenter[i]);
			}
			LCFunctionValue * borderCellShapeInfo;
			parametricShape_->evalShapeInfo(parametricShape_->mapFromStandardHypercube(borderCellCenter), &borderCellShapeInfo);
			LCErrorReturn(borderCell->createHomeomorphicMap(borderCellShapeInfo, cellSamples));


			for (auto neighbor : neighbours)
			{
				neighbor->addNeighbor(borderCell);
				borderCell->addNeighbor(neighbor);
			}
			borderCells_.push_back(borderCell);
			int combinationID = 0;
			for (int i = 0; i < nParams; i++)
			{
				combinationID += corner(i)*std::pow(nSamples + 2, i);
			}
			if (cornerCellsFromCombination.find(combinationID) != cornerCellsFromCombination.end())
			{
				std::cout << "error getting combiantionID " << std::endl;
				system("pause");
			}
			cornerCellsFromCombination[combinationID] = borderCell;
		}

		std::vector<Eigen::VectorXi> cornerNeighborCombinations;
		LCErrorReturn(LCMathHelper::computeCombinations(nParams, 3, &cornerNeighborCombinations));


		for (int i = 0; i < borderCells_.size(); i++)
		{
			for (auto cornerNeighbor : cornerNeighborCombinations)
			{
				Eigen::VectorXi cnVec = (cornerNeighbor - Eigen::VectorXi::Ones(nParams));
				if (cnVec.norm() == 0)
				{
					continue;
				}
				cnVec += cornerCellCombiantions[i];
				int combinationID = 0;
				if (cnVec.minCoeff() < 0)
				{
					continue;
				}
				for (int i = 0; i < nParams; i++)
				{
					combinationID += cnVec(i)*std::pow(nSamples + 2, i);
				}
				auto neighbor = cornerCellsFromCombination.find(combinationID);
				if (neighbor != cornerCellsFromCombination.end())
				{
					borderCells_[i]->addNeighbor(neighbor->second);
				}
			}
		}
	}


	return err;
}


LCError LCAdaptiveGrid::getSample(const Eigen::VectorXd &center, LCSample** result)
{
	for (auto sample : samples_)
	{
		if ((sample->getCenter() - center).norm() < LCMathHelper::EPSILON)
		{
			*result = sample;
			return LCError();
		}
	}
	std::stringstream error;
	error << "sample " << center.transpose() << " not found";
	return LCError(error.str());
}

void LCAdaptiveGrid::log()
{

	std::cout << "these are the adaptive grid nodes:" << std::endl;
	root_->log(0);
}

LCError LCAdaptiveGrid::adaptiveSample()
{
	LCError err;
	std::vector<LCAdaptiveGridCell*> leafCells;
	root_->getAllLeafCells(&leafCells);
	std::unordered_map<LCAdaptiveGridCell*, bool> hasBeenProcessed;
	auto cmp = [](LCAdaptiveGridCell* left, LCAdaptiveGridCell* right) {
		return left->size() < right->size();
	};
	std::priority_queue<LCAdaptiveGridCell*, std::vector<LCAdaptiveGridCell*>, decltype(cmp)> cellProcessingQueue(cmp);

	for (auto leafCell : leafCells)
	{
		hasBeenProcessed[leafCell] = false;
		cellProcessingQueue.push(leafCell);
	}

	while (!cellProcessingQueue.empty())
	{
		LCAdaptiveGridCell* cell = cellProcessingQueue.top();
		std::cout << "checking cell ";  cell->logRanges(0);
		cellProcessingQueue.pop();
		auto processingStatus = hasBeenProcessed.find(cell);
		if (processingStatus != hasBeenProcessed.end())
		{
			if (processingStatus->second)
			{
				continue;
			}
		}

		bool doSplit;
		LCErrorReturn(decideSplit(cell, &doSplit));
		//std::cout << "deciding to split a cell that points in direction " << cell->getSplitDirection() << std::endl;
		if (doSplit)
		{
			int nextDirection;
			LCErrorReturn(decideSplitDirection(cell, &nextDirection));
			std::vector<LCAdaptiveGridCell*> affectedCells;
			err = refineCell(cell, nextDirection, &affectedCells);
			cell->deleteAllSamples();
			if (!err.isOK())
			{
				return err;
			}

			//std::cout << "pushing a cell that points in direction " << cell->getChild1()->getSplitDirection() << std::endl;

			cellProcessingQueue.push(cell->getChild1());
			hasBeenProcessed[cell->getChild1()] = false;
			cellProcessingQueue.push(cell->getChild2());
			hasBeenProcessed[cell->getChild2()] = false;
			for (auto affectedCell : affectedCells)
			{
				if (affectedCell->getLevel() > 0)
				{
					cellProcessingQueue.push(affectedCell);
					hasBeenProcessed[affectedCell] = false;
				}

			}
		}


		hasBeenProcessed[cell] = true;
	}

	return err;
}

LCError LCAdaptiveGrid::decideSplitDirection(LCAdaptiveGridCell *cell, int *direction)
{
	switch (params_->splitPolicy_)
	{
	case LCAdaptiveSamplingParams::LCSplitPolicy::CYCLIC:
		// an alternative is to always split in a direction orthogonal to the previous split
		*direction = (cell->getParent()->getSplitDirection() + 1) % 2;
	//	std::cout << "decide split direction cyclic" << params_->splitPolicy_ << std::endl;
		break;
	case LCAdaptiveSamplingParams::LCSplitPolicy::ERROR_BASED:
		LCErrorReturn(cell->getOptimalSplitDirection(direction));
//		std::cout << "decide split direction errorBased" << params_->splitPolicy_ << std::endl;
		break;
	default:
		//std::cout << "decide split direction default" << params_->splitPolicy_ << std::endl;
		break;
	}
	return LCError();
}

LCError getJitterPoint(LCAdaptiveGridCell *cell, std::vector<double> *result)
{
	std::vector<double> jitterPoint;
	for (int i = 0; i < cell->ranges_.size() / 2; i++)
	{
		float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		double pos = cell->ranges_[2 * i] + r*(cell->ranges_[2 * i + 1] - cell->ranges_[2 * i]);
		jitterPoint.push_back(pos);
	}
	*result = jitterPoint;
	return LCError();
}


LCError LCAdaptiveGrid::decideSplit(LCAdaptiveGridCell* cell, bool *splitDecision)
{
	LCError err;

	//decide if split is necessary
	bool doSplit = cell->getLevel() < maxLevel_;
	if (doSplit)
	{
		std::vector<double> center;
		cell->getCenter(&center);
		Eigen::VectorXd  diff;
		LCFunctionValue* approxShapeInfo = nullptr;
		LCFunctionValue* actualShapeInfo = cell->getCenterShapeInfo();

		switch (params_->errorMetric_)
		{
		case LCAdaptiveSamplingParams::LCErrorMetric::CENTER_BASED:
			LCErrorReturn(cell->evalShapeInfo(center, basisType_, &approxShapeInfo));
			LCErrorReturn(approxShapeInfo->computeErrorVector(actualShapeInfo, &diff));
			if (diff.size() != allowedError_.size())
			{
				return LCError("the size of the error vector and allowed dif are not the same");
			}
			doSplit = !((diff - allowedError_).maxCoeff() <= 0);
			break;
		case LCAdaptiveSamplingParams::LCErrorMetric::RANDOM_SAMPLES:
			for (int i = 0; i < params_->nErrorSamples_; i++)
			{
				std::vector<double> otherPoint;
				getJitterPoint(cell, &otherPoint); //get random point in cell
				LCFunctionValue* otherActualShapeInfo;
				std::vector<double> mappedOtherPoint = parametricShape_->mapFromStandardHypercube(otherPoint);
				parametricShape_->evalShapeInfo(mappedOtherPoint, &otherActualShapeInfo); //get actual shape info of jitter point
				LCErrorReturn(cell->evalShapeInfo(otherPoint, basisType_, &approxShapeInfo));//get approximation 
				LCErrorReturn(approxShapeInfo->computeErrorVector(otherActualShapeInfo, &diff));
				if (diff.size() != allowedError_.size())
				{
					return LCError("the size of the error vector and allowed dif are not the same");
				}
				doSplit = doSplit || !((diff - allowedError_).maxCoeff() <= 0);
			}
			break;
		default:
			break;
		}
	}
	*splitDecision = doSplit;
	return err;
}


LCError LCAdaptiveGrid::refineCell(LCAdaptiveGridCell* cell, int splitDirection, std::vector<LCAdaptiveGridCell*> *affectedCells)
{
	//std::cout << "before split" << std::endl;
	//checkViolatesLocality();

	///czw debug
	std::cout << "going to split cell " << std::endl;
	cell->logRanges(0);
	//checkViolatesLocality();

	/*std::vector<double> ranges;
	ranges = cell->ranges_;
	if (ranges[0] == 0.5 && ranges[1] == 1 && ranges[2] == 0.875 && ranges[3] == 1)
	{
	std::cout << "refining the troublesome cell " << std::endl;
	cell->log(0);
	for (auto sample : samples_)
	{
	Eigen::VectorXd samCenter = sample->getCenter();
	if (samCenter(0) == 0.5 && samCenter(1) == 0.75)
	{
	std::cout << "found the troublesome sample" << std::endl;
	std::cout << "troublesome sample basis function size: " << sample->getBasisFunctions().size() << std::endl;
	sample->log();
	}
	}
	}*/
	////czw debug

	LCError err;
	//--- step1: add new samples
	std::vector<Eigen::VectorXd> midValues;
	std::vector<LCAdaptiveGridCell*> adjacentCells;
	cell->getMidValues(splitDirection, &midValues);

	std::vector<LCSample*> midSamples;
	std::vector<LCSample*> newMidSamples;
	std::vector<LCSample*> extraMidSamples;

	for (auto midValue : midValues)
	{
		//std::cout << "mid sample at: " << midValue.transpose() << std::endl;
		LCSample* midSample;
		//cell->logSamples();
		if (!cell->getSample(midValue, &midSample).isOK())
		{
			LCFunctionValue *shapeInfo = nullptr;
			std::vector<double> midValueVec;
			LCMathHelper::eigen2StdVector(midValue, &midValueVec);
			parametricShape_->evalShapeInfo(parametricShape_->mapFromStandardHypercube(midValueVec), &shapeInfo);
			midSample = new LCSample(midValue, shapeInfo);
			samples_.push_back(midSample);
			newMidSamples.push_back(midSample);
			//std::cout << "----added mid sample" << std::endl;
		}
		midSamples.push_back(midSample);
	}

	//--- step2: refine the basis function of the samples affecting the cell distributing	 weights to the new samples
	Eigen::VectorXd cellSize;
	cell->getCellSize(&cellSize);
	Eigen::VectorXd desiredSupport = cellSize;
	int projectedSize = cellSize.size() - 1;
	std::vector<LCSample*> samplesAffectingCell;
	switch (basisType_)
	{
	case LCBasisFunction::LCBasisFunctionType::LINEAR_BSPLINE:
		desiredSupport(splitDirection) = 0.5 * desiredSupport(splitDirection);
		cell->getSamples(&samplesAffectingCell);
		break;
	case LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE:
	{
		desiredSupport(splitDirection) = 0.5 * desiredSupport(splitDirection);
		desiredSupport = 2 * desiredSupport;
		cell->getDoubleAdjacentSamples(&samplesAffectingCell);
		break;
	}
	break;
	default:
		return LCError("basis type unknown");
		break;
	}

	//refine all the basis functions
	std::vector<std::pair<LCBasisFunction*, LCSample*>> removedFunctions;
	//std::cout << " desiredSupport(splitDirection)" << desiredSupport(splitDirection) << std::endl;
	int nRefines = 0;
	for (auto sample : samplesAffectingCell)
	{
		//std::cout << "sample affecting cell " << sample->getCenter().transpose() << std::endl;
		LCErrorReturn(sample->refineBasisFunctions(cell, splitDirection, desiredSupport(splitDirection), basisType_, &nRefines));
	}

	//--- step3: split the cell and create leaf information (homeomorphicMap on it)
	cell->split(splitDirection);
	std::vector<LCSample*> samples = samplesAffectingCell;
	samples.reserve(samples.size() + newMidSamples.size());
	samples.insert(samples.end(), newMidSamples.begin(), newMidSamples.end());
	setLeafInformation(cell->getChild1(), samples);
	setLeafInformation(cell->getChild2(), samples);

	//--- step4: get the leaf nodes that are affected by each new mid samples and add the sample
	for (auto neighbor : cell->getNeighbors())
	{
		for (auto midSample : newMidSamples)
		{
			if (neighbor->cointainsPoint(midSample->getCenter()))
			{
				neighbor->addSample(midSample);
			}
		}
	}


	// remove all the basis functions that do not get us locality 

	/*
	std::set<LCAdaptiveGridCell*> totalSupport;
	for (auto adj : cell->getChild1()->getNeighbors())
	{
	totalSupport.insert(adj);
	}
	for (auto adj : cell->getChild2()->getNeighbors())
	{
	totalSupport.insert(adj);
	}
	totalSupport.insert(cell->getChild1());
	totalSupport.insert(cell->getChild2());
	*/

	//czw debug
	std::set<LCAdaptiveGridCell*> totalSupport;
	this->getAllCellsInTheWorld(root_, &totalSupport);
	std::vector<LCAdaptiveGridCell*> padding = getBorderCells();
	for (auto padCell : padding)
	{
		totalSupport.insert(padCell);
	}
	//czw debug

	if (basisType_ == LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE)
	{
		std::vector<LCAdaptiveGridCell*> auxSupport(totalSupport.begin(), totalSupport.end());
		for (auto aux : auxSupport)
		{
			for (auto adj : aux->getNeighbors())
			{
				totalSupport.insert(adj);
			}
		}
	}

	int nBreakLocalityn = 0;
	//for each sample that is affecting the cell being refined, check if that sample is violating locality. 
	//Take the basis functions associated with that sample and remove them from the sample. Add them to a list. 
	for (auto sample : samplesAffectingCell)
	{
		/*if (sample->getCenter()(0) == -0.125 && sample->getCenter()(1) == 0)
		{
			std::cout << "check 1" << std::endl;
			sample->log();
			system("pause");
		}*/
	}
	for (auto sample : samplesAffectingCell)
	{

		//find all of the cells that are on the neighborhood of the sample after this split 
		std::unordered_set<LCAdaptiveGridCell*> adjacentCells;
		std::vector<double>params;
		LCMathHelper::eigen2StdVector(sample->getCenter(), &params);
		if (basisType_ == LCBasisFunction::LCBasisFunctionType::LINEAR_BSPLINE)
		{
			LCErrorReturn(getAdjacentCells(params, &adjacentCells));
		}
		if (basisType_ == LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE)
		{
			LCErrorReturn(getDoubleAdjacentCells(params, &adjacentCells));
		}

		//check if the basis function in non zero in any cell that is not in the adjacentCells
		//So that we don't have to check the whole world, we assume that locality was preserverd
		//in the previous iteration. And therefore we only need to check the totalSupprt and not the whole world
		std::vector<LCBasisFunction*> remBasisFunctions;
		for (auto basisFunction : sample->getBasisFunctions())
		{
			Eigen::VectorXd center = basisFunction.second->getCenter();
			Eigen::VectorXd support = basisFunction.second->getSupport();
			Eigen::VectorXd sampleCenter = sample->getCenter();
		
			bool localityIsBroken = false;
			for (auto supportCell : totalSupport)
			{
				//if a cell from this "reduced world" is not an adjacent cell
				if (center(0) == -0.25 && center(1) == 0 && support(0) == 0.125 && support(1) == 0.125 
					&& supportCell->ranges_[0] == -0.375 && supportCell->ranges_[1] == -0.3125 && supportCell->ranges_[2] == 0 && supportCell->ranges_[3] == 0.5)
				{
				//	std::cout << "here" << std::endl;
				//	std::cout << sampleCenter.transpose() << std::endl;
				//	basisFunction.second->log();
				//	supportCell->log(0);
				//	system("pause");
				}
				auto adj = adjacentCells.find(supportCell);
				if (adj == adjacentCells.end())
				{
					//check if the basis function overlaps this cell
					if (basisFunction.second->checkRegionOverlap(supportCell->ranges_))
					{
						localityIsBroken = true;
						break;
					}
				}
			}

			if (localityIsBroken)
			{
				nBreakLocalityn++;
				std::unordered_set<LCAdaptiveGridCell*> aux;
				std::vector<double> basisCenter;
				LCMathHelper::eigen2StdVector(basisFunction.second->getCenter(), &basisCenter);
				if (!this->getAdjacentCells(basisCenter, &aux).isOK())
				{
					std::cout << "locality is broken for somthing that it shouldn't be" << std::endl;
					basisFunction.second->log();
					std::cout << "support cells" << std::endl;
					for (auto supportCell : totalSupport)
					{
						supportCell->log(0);
					}
					std::cout << "adj cell sfor sample " << sample->getCenter().transpose() << std::endl;
					for (auto supportCell : adjacentCells)
					{
						supportCell->log(0);
					}
					for (auto supportCell : totalSupport)
					{
						auto adj = adjacentCells.find(supportCell);
						if (adj == adjacentCells.end())
						{
							if (basisFunction.second->checkRegionOverlap(supportCell->ranges_))
							{
								std::cout << "------locality is broken for cell"; supportCell->logRanges(0);
							}
						}
					}
					std::unordered_set<LCAdaptiveGridCell*> adjacentCells2;
					getDoubleAdjacentCells(params, &adjacentCells2);
					system("pause");
				}
				remBasisFunctions.push_back(basisFunction.second);
			}
		}

		for (auto basisFunction : remBasisFunctions)
		{
			removedFunctions.push_back(std::make_pair(basisFunction, sample));
			sample->removeBasisFucntion(basisFunction);
		}
	}

	//for (auto sample : samplesAffectingCell)
	//{
	//	if (sample->getCenter()(0) == -0.125 && sample->getCenter()(1) == 0)
	//	{
	//		std::cout << "check 2" << std::endl;
	//		sample->log();
	//		system("pause");
	//	}
	//}

	//Group all of the removed basis functions for reallocation 
	std::unordered_map<LCBasisFunctionKey, GroupFunctionInfo*> groupedBasisFunctions;
	for (auto basisFunction : removedFunctions)
	{

		LCBasisFunctionKey bfKey(basisFunction.first);
		GroupFunctionInfo* groupedBasisFunction = nullptr;
		auto bfGroupMap = groupedBasisFunctions.find(bfKey);
		if (bfGroupMap != groupedBasisFunctions.end())
		{
			groupedBasisFunction = bfGroupMap->second;
		}
		else{
			LCBasisFunction* newBasisFunction;
			err = LCBasisFunction::newFromCopyWithWeights(basisFunction.first, 0, &newBasisFunction);
			if (!err.isOK())
			{
				std::cout << err.internalDescription() << std::endl;
				return err;
			}
			groupedBasisFunction = new GroupFunctionInfo();
			groupedBasisFunctions[bfKey] = groupedBasisFunction;
		}

	/*	if (basisFunction.second->getCenter()(1) < 0.005 && basisFunction.second->getCenter()(1) > 0)
		{
			std::cout << "found a small grouped function while removing" << std::endl;
			std::cout << basisFunction.second->getCenter().transpose() << std::endl;
			system("pause");
		}*/

		groupedBasisFunction->add(basisFunction.second, basisFunction.first->getWeight());

		delete basisFunction.first;
	}
	//for (auto sample : samplesAffectingCell)
	//{
	//	if (sample->getCenter()(0) == -0.125 && sample->getCenter()(1) == 0)
	//	{
	//		std::cout << "check 3" << std::endl;
	//		sample->log();
	//		system("pause");
	//	}
	//}

	//Looping through all the samples and adding it to the groupedBasisFunctions
	for (auto sample : samples_)//remember to change this back
	{
		err = (sample->extractGroupBasisFunctions(&groupedBasisFunctions));
		if (!err.isOK())
		{
			std::cout << err.internalDescription() << std::endl;
			return err;
		}

	}

	//for (auto sample : samplesAffectingCell)
	//{
	//	if (sample->getCenter()(0) == -0.125 && sample->getCenter()(1) == 0)
	//	{
	//		std::cout << "check 4" << std::endl;
	//		sample->log();
	//		system("pause");
	//	}
	//}

	for (auto basisFunction : groupedBasisFunctions)
	{
		std::vector<double> basisCenter;
		/*if (basisFunction.second->getCenter()(1) < 0.005 && basisFunction.second->getCenter()(1) > 0)
		{
			std::cout << "found a small grouped function" << std::endl;
			std::cout << basisFunction.second->getCenter().transpose() << std::endl;
			system("pause");
		}*/
		LCMathHelper::eigen2StdVector(basisFunction.second->getCenter(), &basisCenter);
		std::vector<std::pair<LCSample*, double>> replacingWeights;
	/*	if (basisCenter[0] == -0.25 && basisCenter[1] == 0)
		{
			std::cout << " basis of interest " << std::endl;
			system("pause");
		}*/
		err = (getReplacingSamples(basisCenter, basisFunction.second->getCenter(), &extraMidSamples, &replacingWeights));
		if (!err.isOK())
		{
			std::cout << err.internalDescription() << std::endl;
			return err;
		}

		for (auto replacingWeight : replacingWeights)
		{
			if (replacingWeight.second > 0)
			{
				/*if (replacingWeight.first->getCenter()(0) == -0.125 && replacingWeight.first->getCenter()(1) == 0)
				{
					std::cout << "here in the replacing weihgt" << std::endl;
					replacingWeight.first->log();
					system("pause");
				}*/

				LCBasisFunction* newBasisFunction;
				err = (LCBasisFunction::newFromCopyKey(basisFunction.first, basisFunction.second->getWeight()* replacingWeight.second, &newBasisFunction));
				if (!err.isOK())
				{
					std::cout << err.internalDescription() << std::endl;
					return err;
				}
				replacingWeight.first->addBasisFunction(newBasisFunction);
			}
		}
		delete basisFunction.second;
	}

	//for (auto sample : samplesAffectingCell)
	//{
	//	if (sample->getCenter()(0) == -0.125 && sample->getCenter()(1) == 0)
	//	{
	//		std::cout << "check 5" << std::endl;
	//		sample->log();
	//		system("pause");
	//	}
	//}

	//--- step 5: add the extra mid samples
	for (auto midSample : extraMidSamples)
	{
		cell->getChild1()->addSample(midSample);
		cell->getChild2()->addSample(midSample);

	}
	for (auto neighbor : cell->getNeighbors())
	{
		for (auto midSample : extraMidSamples)
		{
			if (neighbor->cointainsPoint(midSample->getCenter()))
			{
				neighbor->addSample(midSample);
			}
		}
	}

	std::cout << "finish split" << std::endl;
	//cell->getChild1()->log(0);
	//cell->getChild2()->log(0);
	//checkViolatesLocality();

	return err;
}


LCError LCAdaptiveGrid::getAllCellsInTheWorld(LCAdaptiveGridCell *root, std::set<LCAdaptiveGridCell*> *result)
{
	if (root->isLeaf())
	{
		(*result).insert(root);
	}
	if (!root->isLeaf())
	{
		getAllCellsInTheWorld(root->getChild1(), result);
		getAllCellsInTheWorld(root->getChild2(), result);
	}
	return LCError();
}

LCError logAncestors(LCAdaptiveGridCell *cell)
{
	std::cout << "cell " << std::endl;
	cell->logRanges(0);
	std::cout << "cell split direction " << cell->getSplitDirection() << std::endl;
	if (cell->getParent() != NULL)
	{
		logAncestors(cell->getParent());
	}
	return LCError();
}

LCError LCAdaptiveGrid::checkViolatesLocality()
{
	LCError err;
	std::set<LCAdaptiveGridCell*> totalSupport;
	this->getAllCellsInTheWorld(root_, &totalSupport);
	std::vector<LCAdaptiveGridCell*> padding = getBorderCells();
	for (auto padCell : padding)
	{
		totalSupport.insert(padCell);
	}

	std::cout << "all cells in world " << totalSupport.size() << std::endl;
	bool localityIsBroken = false;
	for (auto sample : this->samples_)
	{
		std::unordered_set<LCAdaptiveGridCell*> adjacentCells;
		std::vector<double>params;
		LCMathHelper::eigen2StdVector(sample->getCenter(), &params);
		if (basisType_ == LCBasisFunction::LCBasisFunctionType::LINEAR_BSPLINE)
		{
			LCErrorReturn(getAdjacentCells(params, &adjacentCells));
		}
		if (basisType_ == LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE)
		{
			LCErrorReturn(getDoubleAdjacentCells(params, &adjacentCells));
		}

		std::vector<LCBasisFunction*> remBasisFunctions;
		for (auto basisFunction : sample->getBasisFunctions())
		{
			for (auto supportCell : totalSupport)
			{
				auto adj = adjacentCells.find(supportCell);
				if (adj == adjacentCells.end())
				{
					if (basisFunction.second->checkRegionOverlap(supportCell->ranges_))
					{
						localityIsBroken = true;
						std::cout << "locality violated for sample: " << sample->getCenter().transpose() << std::endl;
						sample->log();
						//system("pause");
						std::cout << "cell which is influenced by the sample, but shouldn't be" << std::endl;
						supportCell->log(0);
					/*	std::cout << "sample's adjacent cells" << std::endl;
						for (auto adjCell : adjacentCells)
						{
							adjCell->logRanges(0);
						}*/
						//basisFunction.second->checkRegionOverlap(supportCell->ranges_);
						//logAncestors(supportCell);
						system("pause");
						basisFunction.second->log();
						//root_->log(0);
						err = LCError("locality is broken");
						//break;
					}
				}
			}
		}
	}
	if (!localityIsBroken)
	{
		std::cout << "locality not violated " << std::endl;
	}
	return err;
}


LCError LCAdaptiveGrid::getReplacingSamples(std::vector<double> basisCenter, Eigen::VectorXd replacementCenter,
	std::vector<LCSample*> * extraMidSamples, std::vector<std::pair<LCSample*, double>> *replacingWeights)
{
	LCError err;
	std::unordered_set<LCAdaptiveGridCell*> adjacenteCells;


	err = this->getAdjacentCells(basisCenter, &adjacenteCells);

	if (adjacenteCells.size() < 1)
	{
		std::cout << "cell has no adjacent cells " << std::endl;
		std::cout << basisCenter.size() << std::endl;
		std::cout << basisCenter[0] << " " << basisCenter[1] << std::endl;
		return LCError("adjacent cells is empty");
	}
	auto itr = adjacenteCells.begin();
	std::vector<double> interRanges = (*itr)->ranges_;
	++itr;
	for (; itr != adjacenteCells.end(); ++itr)
	{
		std::vector<double> auxIterRanges;
		LCErrorReturn(LCAdaptiveGridCell::getIntersection(interRanges, (*itr)->ranges_, &auxIterRanges));
		interRanges = auxIterRanges;
	}
	int nParams = root_->ranges_.size() / 2;
	Eigen::VectorXd minRange(nParams);
	Eigen::VectorXd maxRange(nParams);
	for (int i = 0; i < nParams; i++)
	{
		minRange(i) = interRanges[2 * i];
		maxRange(i) = interRanges[2 * i + 1];
	}

	//check if replacementCenter is inside the region
	bool canDoLinearPrecision = true;
	for (int i = 0; i < nParams; i++)
	{
		if ((replacementCenter(i) < (minRange(i) - LCMathHelper::EPSILON)) ||
			(replacementCenter(i) > (maxRange(i) + LCMathHelper::EPSILON)))
		{
			std::cout << "cannot do linear precision" << std::endl;
			std::cout << "inter ranges ";
			for (int i = 0; i < interRanges.size(); i++)
			{
				std::cout << interRanges[i] << " ";
			}
			std::cout << std::endl;
			std::cout << "basisCenter: ";
			for (int i = 0; i < basisCenter.size(); i++)
			{
				std::cout << basisCenter[i] << " ";
			}
			std::cout << std::endl;
			std::cout << "replacement center = " << replacementCenter.transpose() << std::endl;
			//return LCError("cannot do linear precsion");
			system("pause");
			canDoLinearPrecision = false;
			break;
		}
	}

	if (!canDoLinearPrecision)
	{
		LCMathHelper::std2EigenVector(basisCenter, &replacementCenter);
	}


	double sumWeights = 0;
	std::vector<Eigen::VectorXi> combinations;
	LCMathHelper::computeCombinations(nParams, 2, &combinations);
	for (auto combination : combinations)
	{
		double weight = 1;
		for (int i = 0; i < nParams; i++)
		{
			double alphaParam = 1;
			double area = interRanges[2 * i + 1] - interRanges[2 * i];
			if (abs(area) > LCMathHelper::EPSILON)
			{
				alphaParam = (replacementCenter[i] - interRanges[2 * i]) / area;
			}
			if (combination(i) == 1)
			{
				weight *= alphaParam;
			}
			else
			{
				weight *= (1 - alphaParam);
			}
			/*if (weight < 0.001 && weight > 0)
			{
				std::cout << "\nsmall weight "<< weight << std::endl;
				std::cout << "area " << area << " alpha " << alphaParam << std::endl;
				std::cout << "i " << i << std::endl;
				double product = (replacementCenter[i] - interRanges[2 * i]) / area;
				std::cout << "product " << product << std::endl;
				std::cout << "combination " << combination.transpose() << std::endl;
				Eigen::VectorXd interRangesVec;
				LCMathHelper::std2EigenVector(interRanges, &interRangesVec);
				std::cout << "interRanges " << interRangesVec.transpose() << std::endl;
				std::cout << "replacementCenter " << replacementCenter.transpose() << std::endl;
				Eigen::VectorXd sampleCenter = combination.cast<double>().cwiseProduct(maxRange)
					+ (Eigen::VectorXd::Ones(nParams) - combination.cast<double>()).cwiseProduct(minRange);
				std::cout << "sample center " << sampleCenter.transpose() << std::endl;
				std::cout << "basis center " << basisCenter[0] << " " << basisCenter[1] << std::endl;
				std::cout << "adjacent cells " << std::endl;
				for (auto cell : adjacenteCells)
				{
					cell->log(0);
				}
				system("pause");
			}*/
		}
		if (weight > 0)
		{
			sumWeights += weight;
			Eigen::VectorXd sampleCenter = combination.cast<double>().cwiseProduct(maxRange)
				+ (Eigen::VectorXd::Ones(nParams) - combination.cast<double>()).cwiseProduct(minRange);
			LCSample *sample;
			//if (sampleCenter(0) == -0.125 && sampleCenter(1) == 0)
			//{
			//	std::cout << "sample center " << sampleCenter.transpose() << std::endl;
			//	std::cout << "combination " << sampleCenter.transpose() << std::endl;
			//	std::cout << "basis center" << basisCenter[0] << " " << basisCenter[1] << std::endl;
			//	std::cout << "weight " << weight << std::endl;
			//	std::cout << "adjacent cells " << std::endl;
			//	for (auto cell : adjacenteCells)
			//	{
			//		cell->log(0);
			//	}
			//	system("pause");
			//}
			err = (getSample(sampleCenter, &sample));
			if (!err.isOK())
			{
				// adding new mid sample
				LCFunctionValue *shapeInfo = nullptr;
				std::vector<double> midValueVec;
				LCMathHelper::eigen2StdVector(sampleCenter, &midValueVec);
				parametricShape_->evalShapeInfo(parametricShape_->mapFromStandardHypercube(midValueVec), &shapeInfo);
				sample = new LCSample(sampleCenter, shapeInfo);
				samples_.push_back(sample);//bookmark
				extraMidSamples->push_back(sample);
				//std::cout << " added mid sample at " << sampleCenter.transpose() << std::endl;
			}
			replacingWeights->push_back(std::make_pair(sample, weight));
		}
	}
	if (abs(sumWeights - 1) > LCMathHelper::EPSILON)
	{
		std::cout << "weights do not sum to one in get replacing sample : " << sumWeights << std::endl;

		std::cout << std::endl;
		system("pause");
		return LCError("weights do not sum to one: " + std::to_string(sumWeights));
	}
	return LCError();
}
LCError LCAdaptiveGrid::validateLocality()
{
	LCError err;

	//TODOCUBIC
	if (basisType_ == LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE)
	{
		return err;
	}

	for (int i = 0; i < samples_.size(); i++)
	{
		std::unordered_set<LCAdaptiveGridCell*> adjacentCells;
		std::vector<double>params;
		LCMathHelper::eigen2StdVector(samples_[i]->getCenter(), &params);
		LCErrorReturn(getAdjacentCells(params, &adjacentCells));
		std::vector<LCAdaptiveGridCell*> adjacentCellsVec;
		adjacentCellsVec.insert(adjacentCellsVec.end(), adjacentCells.begin(), adjacentCells.end());
		LCErrorReturn(samples_[i]->validadeIsCoveredByCells(adjacentCellsVec, parametricShape_));
	}
	return err;
}


LCError LCAdaptiveGrid::evalShapeInfoUsingAllSamples(const std::vector<double> &params,
	LCBasisFunction::LCBasisFunctionType basisType, LCFunctionValue **result)
{


	LCError err;
	Eigen::VectorXd pos;
	LCMathHelper::std2EigenVector(params, &pos);
	std::vector<LCFunctionValue*> shapesToDelete;
	std::vector<std::pair<LCSample*, LCFunctionValue*>> influencingSamples;

	for (auto samp : samples_)
	{
		influencingSamples.push_back(std::make_pair(samp, samp->getShapeInfo()));
	}
	//for all samples add to inclueng samples; 


	if (!err.isOK())
	{
		std::cout << err.internalDescription();
		return err;
	}

	std::vector<std::pair<double, LCFunctionValue*>> weightedSamples;
	double weightSum = 0;
	int count = 0;
	for (auto sample : influencingSamples)
	{
		double weight;
		err = sample.first->getInterpolationWeight(pos, &weight);
		if (!err.isOK())
		{
			std::cout << err.internalDescription();
			return err;
		}
		/*if (weight > 0)
		{

		std::cout << "sample center = " << sample.first->getCenter().transpose() << std::endl;
		}*/
		weightSum += weight;
		weightedSamples.push_back(std::make_pair(weight, sample.second));
	}



	if (abs(weightSum - 1) > LCMathHelper::EPSILON)
	{
		std::cout << "params " << params[0] << "' " << params[1] << std::endl;
		std::cout << "weights do not sum to one: " << weightSum << std::endl;
		system("pause");
	}
	err = (LCFunctionValue::newFromWeightedSum(weightedSamples, result));
	if (!err.isOK())
	{
		std::cout << err.internalDescription();
		return err;
	}

	for (auto shape : shapesToDelete)
	{
		delete shape;
	}

	return err;

}
