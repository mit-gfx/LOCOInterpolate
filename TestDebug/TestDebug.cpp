#include <sstream>
#include <iostream>
#include <fstream>
#include "LCMathHelper.h"
#include "LCError.h"
#include "LCRealFunction.h"
#include "LCAdaptiveGrid.h"
#include "LCRealFunctionValue.h"
#include "LCAdaptiveGridCell.h"
#include "LCProtoConverter.h"
#include "LCSampledFunction.h"
#include "LCSample.h"
#include "LCTimer.h"
#include "LCFunctionValueMap.h"
#include "LCAdaptiveSamplingParams.h"
#include "lodepng.h"

//Test whether linear precision is satisfied 
LCError linearPrecisionTestCubic()
{
	int nParams = 3;
	LCError err;
	LCRealFunction* paramShape = new LCRealFunction(nParams, true);
	Eigen::VectorXd allowedError = 0.1 * Eigen::VectorXd::Ones(1);
	LCAdaptiveGrid grid(paramShape, LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE, allowedError, 0, 0, true);

	/* random refinement*/
	srand(1);
	int refined = 0;
	std::vector<LCAdaptiveGridCell*> affectedCells;
	for (int i = 0; i < 12; i++)
	{
		std::vector<LCAdaptiveGridCell*> leafs;
		grid.getRoot()->getAllLeafCells(&leafs);

		int id = rand() % leafs.size();
		int dir = rand() % nParams;
		if (leafs[id]->getLevel() < 5)
		{
			LCErrorReturn(grid.refineCell(leafs[id], dir, &affectedCells));
			refined++;
		}
	}

	grid.log();

	std::vector<Eigen::VectorXi> sampleCombinations;
	int nSamples = 10;
	LCErrorReturn(LCMathHelper::computeCombinations(nParams, nSamples + 1, &sampleCombinations));
	Eigen::VectorXd shapeMinRange(nParams);
	Eigen::VectorXd shapeMaxRange(nParams);
	for (int i = 0; i < nParams; i++)
	{
		shapeMinRange(i) = paramShape->getMinRange(i);
		shapeMaxRange(i) = paramShape->getMaxRange(i);
	}
	bool linearPrecision = true;
	
	return err;
}

LCError getLeafByRange(std::vector<double> ranges, std::vector<LCAdaptiveGridCell*> leafs, LCAdaptiveGridCell **leaf)

{
	for (LCAdaptiveGridCell *cell : leafs)
	{
		std::vector<double> cellRanges = cell->ranges_;
		bool equal = std::equal(ranges.begin(), ranges.end(), cellRanges.begin());
		if (equal)
		{
			*leaf = cell;
		}
	}
	return LCError();
}

LCError selectiveRefinement(LCAdaptiveSamplingParams *params, LCBasisFunction::LCBasisFunctionType basisType)
{
	LCRealFunction* originFunction = new LCRealFunction(2, false);
	LCAdaptiveGrid grid(originFunction, basisType, params, true);

	std::vector<int> splitDirection = {0, 1, 0, 1, 0, 0, 0};
	std::vector<double> vec1 = { -1, 1, -1, 1 };
	std::vector<double> vec2 = { -1, 0, -1, 1 };
	std::vector<double> vec3 = { -1, 0, 0, 1 };
	std::vector<double> vec4 = { -0.5, 0, 0, 1 };
	std::vector<double> vec5 = { -0.5, 0, 0, 0.5 };
	std::vector<double> vec6 = { -0.5, -0.25, 0, 0.5 };
	std::vector<double> vec7 = { -0.375, -0.25, 0, 0.5 };

	std::vector<std::vector<double>> ranges = { vec1, vec2, vec3, vec4, vec5, vec6, vec7 };

	std::vector<LCAdaptiveGridCell*> affectedCells;
	for (int i = 0; i < splitDirection.size(); i++)
	{
		std::vector<LCAdaptiveGridCell*> leafs;
		grid.getRoot()->getAllLeafCells(&leafs);
		LCAdaptiveGridCell *leaf;
		getLeafByRange(ranges[i], leafs, &leaf);
		leaf->logRanges(0);
		int dir = splitDirection[i];
		LCErrorReturn(grid.refineCell(leaf, dir, &affectedCells));
	}
	return LCError();
}

LCError testBugSumToOne(LCAdaptiveSamplingParams *params, LCBasisFunction::LCBasisFunctionType basisType)
{
	LCRealFunction* originFunction = new LCRealFunction(2, false);
	LCAdaptiveGrid grid(originFunction, basisType, params, true);
	std::vector<double> ranges;
	originFunction->getRanges(&ranges);
	LCSampledFunction* approxFunction = new LCSampledFunction(grid.getRoot(), ranges, grid.getSamples(), grid.getBorderCells());

	std::cout << "number of samples " << grid.getSamples().size() << std::endl;

	std::cout << "checking for locality" << std::endl;
	LCError locError = grid.checkViolatesLocality();
	if (!locError.isOK())
	{
	std::cout << "locality is broken!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	}
	system("pause"); 

	int imageResolution = 200;
	double delta = 1.0 / imageResolution;
	int width = imageResolution + 1;
	for (double y = 0; y < width; y += 1)
	{
		for (double x = 0; x < width; x += 1)
		{
			double i = y*delta;
			double j = x*delta;
			std::vector<double> pos(2);
			pos[1] = ranges[0] + i*(ranges[1] - ranges[0]);
			pos[0] = ranges[2] + j*(ranges[3] - ranges[2]);
			LCFunctionValue * info;
			LCError err = approxFunction->evalShapeInfo(pos, &info);
			if (!err.isOK())
			{
				std::cout << "err not ok: " << err.internalDescription() << std::endl;
			}
			LCRealFunctionValue* testInfo = dynamic_cast<LCRealFunctionValue*>(info);
			double temp = testInfo->val;
			int val = (int)10 * testInfo->val;
		}
	}
	return LCError();
}

int main(int argc, char* argv[])
{
	double threshold = 0.1;
	double maxTreeDepth = 10;
	int bootstrap = 2;
	LCAdaptiveSamplingParams *params = new LCAdaptiveSamplingParams(maxTreeDepth, threshold, bootstrap, 0, LCAdaptiveSamplingParams::LCSplitPolicy::ERROR_BASED, LCAdaptiveSamplingParams::LCErrorMetric::CENTER_BASED);
	LCBasisFunction::LCBasisFunctionType basisType = LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE;
	std::cout << "running LOCOInterpolate with the following parameters. See LCAdaptiveSamplingParams.h for details." << std::endl;
	params->log();
	//selectiveRefinement(params, basisType);
	testBugSumToOne(params, basisType);
}
