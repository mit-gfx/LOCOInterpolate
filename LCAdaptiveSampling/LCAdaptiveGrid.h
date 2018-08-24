#ifndef _OC_ADAPTIVE_GRID_
#define _OC_ADATIVE_GRID_

#include <Eigen/Dense>
#include <vector> 
#include "LCBasisFunction.h"
#include <unordered_set>
#include <set>
#include "LCAdaptiveSamplingParams.h"

class LCFunction;
class LCSampledFunction;
class LCAdaptiveGridCell;
class LCSample;
class LCError;
class LCFunctionValue;

class LCAdaptiveGrid{
public:
	LCAdaptiveGrid(LCFunction *parametricShape, LCBasisFunction::LCBasisFunctionType basisType,
		Eigen::VectorXd allowedError, int bootstrapGridSize, int maxLevel, bool evalCorners = false);
	LCAdaptiveGrid(LCFunction *parametricShape, LCSampledFunction *coarseShpe,
		Eigen::VectorXd allowedError, int maxLevel, bool evalCorners = false);
	LCAdaptiveGrid(LCFunction *parametricShape, LCBasisFunction::LCBasisFunctionType basisType,
		LCAdaptiveSamplingParams* params, bool evalCorners = false);
	void log();
	LCError evalShapeInfo(const std::vector<double> &params, LCFunctionValue **result);
	LCError refineCell(LCAdaptiveGridCell* cell, int splitDirection, std::vector<LCAdaptiveGridCell*> *affectedCells);
	LCError getContainingCell(const std::vector<double> &params, LCAdaptiveGridCell **result);
	LCError getSmallestCointainingLeaf(const std::vector<double> &params, LCAdaptiveGridCell **leafCell);
	LCError getAdjacentCells(const std::vector<double> &params, std::unordered_set<LCAdaptiveGridCell*> *result);
	LCError getDoubleAdjacentCells(const std::vector<double> &params, std::unordered_set<LCAdaptiveGridCell*> *result);
	LCAdaptiveGridCell* getRoot();
	std::vector<LCSample*> getSamples();
	LCError getSample(const Eigen::VectorXd &center, LCSample** result);
	std::vector<LCAdaptiveGridCell*> getBorderCells();
	LCError getReplacingSamples(std::vector<double> basisCenter,
		Eigen::VectorXd replacementCenter,
		std::vector<LCSample*> * extraMidSamples,
		std::vector<std::pair<LCSample*, double>>* replacingWeights);
	LCError validateLocality();

	// debug function for locality: evaluate the point but using ALL samples instead of only the LOCAL samples
	LCError evalShapeInfoUsingAllSamples(const std::vector<double> &params,
		LCBasisFunction::LCBasisFunctionType basisType, LCFunctionValue **result);

	//check if ANY of the samples in the grid violate locality
	LCError checkViolatesLocality();

	//Return a list of all the leaf cells in the grid
	LCError getAllCellsInTheWorld(LCAdaptiveGridCell *root, std::set<LCAdaptiveGridCell*> *result);

	//determines the split direction of the cell to be refined
	LCError decideSplitDirection(LCAdaptiveGridCell *cell, int *direction);
private:
	LCAdaptiveSamplingParams *params_;
	LCError bootstrapSample(int nSamples);
	LCError bootstrapUniformSplit(LCAdaptiveGridCell* cell, int  direction, int levels);
	LCError adaptiveSample();
	LCError decideSplit(LCAdaptiveGridCell* cell, bool *splitDecision);
	void setLeafInformation(LCAdaptiveGridCell* leafCell, const std::vector<LCSample*> &samples);
	int maxLevel_;
	std::string shapeName_;
	LCFunction *parametricShape_;
	Eigen::VectorXd allowedError_; // the error threshHold for keeping a sample
	LCAdaptiveGridCell *root_;
	std::vector<LCSample*> samples_;
	LCBasisFunction::LCBasisFunctionType basisType_;
	std::vector<LCAdaptiveGridCell*> borderCells_;
	bool evalCorners_;
};


#endif 
