#include "LCSampledFunction.h"
#include "LCAdaptiveGridCell.h"
#include "LCSample.h"
#include "LCFunctionValue.h"
#include "LCError.h"
#include "LCMathHelper.h"
#include <set>

LCSampledFunction::LCSampledFunction(LCAdaptiveGridCell *root, 
	std::vector<double> ranges,
	std::vector<LCSample*> samples, std::vector<LCAdaptiveGridCell*> borderCells)
{
	root_ = root;
	samples_ = samples;
	basisType_ = (*samples_[0]->getBasisFunctions().begin()->second).getType();
	borderCells_ = borderCells;
	ranges_ = ranges;
}

std::vector<LCAdaptiveGridCell*> LCSampledFunction::getBorderCells()
{
	return borderCells_;
}
LCError LCSampledFunction::addBorderCells()
{
	LCError err;

	for (auto borderCell : borderCells_)
	{
		//add samples and neighbours
		std::unordered_set<LCAdaptiveGridCell*> neighbours;
		for (auto sample : borderCell->getPrecomputedSamples())
		{
			std::vector<double> sampleCenter;
			LCMathHelper::eigen2StdVector(sample.first->getCenter(), &sampleCenter);
			root_->getAdjacentLeaves(sampleCenter, &neighbours);
		}
	
		//addNeighbours
		for (auto neighbor : neighbours)
		{
			neighbor->addNeighbor(borderCell);
			borderCell->addNeighbor(neighbor);
		}		
	}
	return err;
}
int LCSampledFunction::getNParams() const
{
	return ranges_.size() / 2;
}
double LCSampledFunction::getMinRange(int iParam) const
{
	return ranges_[2 * iParam];
}
double LCSampledFunction::getMaxRange(int iParam) const
{
	return ranges_[2 * iParam + 1];
}

LCError LCSampledFunction::evalDeriv(const std::vector<double> &shapeParams, int direction, LCFunctionDeriv **result)
{
	auto params = mapToStandardHypercube(shapeParams);
	LCAdaptiveGridCell *containingCell = nullptr;
	LCErrorReturn(root_->getCointainingLeaf(params, &containingCell));
	return containingCell->evalDeriv(params, basisType_, direction, result);
}

LCError LCSampledFunction::evalShapeInfo(const std::vector<double> &shapeParams, LCFunctionValue **result)
{
	auto params = mapToStandardHypercube(shapeParams);

	LCAdaptiveGridCell *containingCell = nullptr;
	LCErrorReturn(root_->getCointainingLeaf(params, &containingCell));//added by czw params -> shapeParams
	LCErrorReturn(containingCell->evalShapeInfo(params, basisType_, result));

	return LCError();
}

LCAdaptiveGridCell* LCSampledFunction::getRoot()
{
	return root_;
}
std::vector<LCSample*> LCSampledFunction::getSamples()
{
	return samples_;
}
