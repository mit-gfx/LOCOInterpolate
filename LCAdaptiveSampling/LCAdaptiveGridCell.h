#ifndef _OC_ADAPTIVE_GRID_CELL_
#define _OC_ADAPTIVE_GRID_CELL_

#include <Eigen/Dense>
#include <vector> 
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include "LCBasisFunction.h"
#include <list>

class LCSample;
class LCError;
class LCFunctionValue;
class LCFunctionDeriv;
class LCFunctionValueMap;

class LCAdaptiveGridCell
{
public:
	LCAdaptiveGridCell(int level, LCAdaptiveGridCell* parent, std::vector<double> &ranges);
	LCAdaptiveGridCell(int level, LCAdaptiveGridCell* parent, std::vector<double> &ranges,
		LCFunctionValue *centerShapeInfo, std::unordered_map<LCSample*, LCFunctionValue*> &samples);
	LCAdaptiveGridCell(int level, LCAdaptiveGridCell* parent, std::vector<double> &ranges,
		int splitDirection, double splitVal);
	LCAdaptiveGridCell* getRoot();
	void setChildren(LCAdaptiveGridCell *child1, LCAdaptiveGridCell* child2);
	int getLevel();
	bool isLeaf();
	std::list<LCAdaptiveGridCell*> getNeighbors() const;
	LCAdaptiveGridCell* getParent();
	LCError split(int splitDirection);
	bool cointainsPoint(Eigen::VectorXd &point); // check if point is inside or on the border of a cell
	LCAdaptiveGridCell* getChild1();
	LCAdaptiveGridCell* getChild2();
	int getSplitDirection();
	double getSplitValue(); 
	double size();
	void getCenter(std::vector<double> *center);
	double getDirCellSize(int i);
	void getCellSize(Eigen::VectorXd *cellSize);
	void getAllLeafCells(std::vector<LCAdaptiveGridCell*> *leafCells);
	LCError getCointainingLeaf(const std::vector<double> &params, LCAdaptiveGridCell **leafCell);
	LCError getSmallestCointainingLeaf(const std::vector<double> &params, LCAdaptiveGridCell **leafCell);
	LCError getAdjacentLeaves(const std::vector<double> &params, std::unordered_set<LCAdaptiveGridCell*> *result);
	LCError getDoubleAdjacentLeaves(const std::vector<double> &params, std::unordered_set<LCAdaptiveGridCell*> *result);
	LCError getNeighborsOrthogonalToSplit(int splitDirection, std::vector<LCAdaptiveGridCell*> *result);
	static LCError getIntersection(const std::vector<double> &ranges1, const std::vector<double> &ranges2, std::vector<double> * result);
	LCError createHomeomorphicMap(LCFunctionValue *centerShapeInfo, std::vector<LCSample*> samples);
	LCError addSample(LCSample* sample);
	LCError addSample(LCSample* sample, LCFunctionValue* shapeInfo);
	LCError evalPametersAreInCell(const std::vector<double> & params);
	LCError getShapeInfoForSampleCenter(Eigen::VectorXd &center, LCFunctionValue **shapeInfo);
	LCError getClosestShapeInfo(Eigen::VectorXd &center, LCFunctionValue **shapeInfo);
	LCError getShapeInfoForSample(LCSample* sample, LCFunctionValue **shapeInfo);
	LCError getNeightbourShapeInfoMappingMapping(LCAdaptiveGridCell* neighbor, Eigen::VectorXd center, LCFunctionValueMap **shapeInfoMapping);
	void logAllSamples(std::string name);
	LCError getAllInfluencingSamples(LCBasisFunction::LCBasisFunctionType basisType,
		std::vector<std::pair<LCSample*, LCFunctionValue*>> *influencingSamples,
		std::vector<LCFunctionValue*> *shapesToDelete);
	LCError evalShapeInfo(const std::vector<double> &params, LCBasisFunction::LCBasisFunctionType basisType, LCFunctionValue **result);
	LCError evalDeriv(const std::vector<double> &params, LCBasisFunction::LCBasisFunctionType basisType, int direction, LCFunctionDeriv **result);
	LCError getMidValues(int splitDirection, std::vector<Eigen::VectorXd> *midValues);
	LCError getSample(const Eigen::VectorXd &center, LCSample** result);
	void getSamples(std::vector<LCSample*> *samples);
	void getDoubleAdjacentSamples(std::vector<LCSample*> *samples);
	LCFunctionValue *getCenterShapeInfo();
	void setCenterShapeInfo(LCFunctionValue * newShapeInfo);
	std::unordered_map<LCSample*, LCFunctionValue*> getPrecomputedSamples();
	LCError getAffectedRangeWhenSplitting(const Eigen::VectorXd &center, int dir, 
		LCBasisFunction::LCBasisFunctionType basisType, std::vector<std::vector<double>> *uncoverdRegion);
	void log(int level);
	void logRanges(int level);
	void logNeightbors(int level);
	void logSamples();
	LCError addAllNeighboursTopDown(); // to be called after protpo read in
	LCError getOptimalSplitDirection(int *optimalDirection);
	LCError getCornerId(Eigen::VectorXd &point, int *cornerID);
	void addNeighbor(LCAdaptiveGridCell* neighbor);
	// getOptimalSplitFromInterpolation(int *optimalDirection);
	LCError computeReplacingWeights(const Eigen::VectorXd &basisCenter, std::vector<std::pair<LCSample*, double>> *replacingWeights);
	LCError getApproximationError(Eigen::VectorXd *approxError);
	LCBasisFunction::LCBasisFunctionType getBasisType();
	void deleteAllSamples();
	std::vector<double> ranges_;
private:
	bool isLeaf_;
	int level_;
	LCAdaptiveGridCell* parent_;
	std::list<LCAdaptiveGridCell*> neighbors_;

	//intermediate node information
	LCAdaptiveGridCell *child1_;
	LCAdaptiveGridCell *child2_;
	int splitDirection_;
	double splitVal_;

	//Leaf node information
	LCFunctionValue *centerShapeInfo_;
	std::unordered_map<LCSample*, LCFunctionValue*> samples_;
	std::vector<LCAdaptiveGridCell*> borderCells_;
	//neighbor stuff
	bool checkIsNeighbor(const LCAdaptiveGridCell* other);
	LCError replaceNeighborWithChildren(LCAdaptiveGridCell* neighbor, LCAdaptiveGridCell* child1, LCAdaptiveGridCell* child2);

	std::map<std::pair<LCFunctionValue*, LCFunctionValue*>, LCFunctionValueMap*> storedMaps;

	double const CONTAINING_EPSILON = 0.0000001;
};


#endif 
