#include <iostream>
#include <fstream>
#include "LCAdaptiveGridCell.h"
#include "LCSample.h"
#include "LCFunctionValue.h"
#include "LCFunctionDeriv.h"
#include "LCError.h"
#include "LCBasisFunction.h"
#include "LCMathHelper.h"
#include "LCFunctionValueMap.h"
#include <unordered_map>
#include <algorithm>  

LCAdaptiveGridCell::LCAdaptiveGridCell(int level, LCAdaptiveGridCell* parent, std::vector<double> &ranges)
{
	isLeaf_ = true;
	child1_ = nullptr;
	child2_ = nullptr;
	ranges_ = ranges;
	level_ = level;
	parent_ = parent;

	//set neighbors
	if (parent_ != nullptr)
	{
		for (auto neighbor : parent_->getNeighbors())
		{
			if (checkIsNeighbor(neighbor))
			{
				neighbors_.push_back(neighbor);
			}
		}
	}

}


LCAdaptiveGridCell::LCAdaptiveGridCell(int level, LCAdaptiveGridCell* parent, std::vector<double> &ranges,
	LCFunctionValue *centerShapeInfo, std::unordered_map<LCSample*, LCFunctionValue*> &samples)
{
	isLeaf_ = true;
	child1_ = nullptr;
	child2_ = nullptr;
	ranges_ = ranges;
	level_ = level;
	parent_ = parent;
	centerShapeInfo_ = centerShapeInfo;
	samples_ = samples;
}

LCAdaptiveGridCell::LCAdaptiveGridCell(int level, LCAdaptiveGridCell* parent, std::vector<double> &ranges,
	int splitDirection, double splitVal)
{
	isLeaf_ = false;
	child1_ = nullptr;
	child2_ = nullptr;
	ranges_ = ranges;
	level_ = level;
	parent_ = parent;
	splitDirection_ = splitDirection;
	splitVal_ = splitVal;

}


LCAdaptiveGridCell* LCAdaptiveGridCell::getRoot()
{
	if (parent_ == nullptr)
	{
		return this;
	}
	return parent_->getRoot();
}

void LCAdaptiveGridCell::setChildren(LCAdaptiveGridCell *child1, LCAdaptiveGridCell* child2)
{
	child1_ = child1;
	child2_ = child2;
}


int LCAdaptiveGridCell::getLevel()
{
	return level_;
}

bool LCAdaptiveGridCell::isLeaf()
{
	return isLeaf_;
}


std::list<LCAdaptiveGridCell*> LCAdaptiveGridCell::getNeighbors() const
{
	return neighbors_;
}

LCAdaptiveGridCell* LCAdaptiveGridCell::getParent()
{
	return parent_;
}

LCError LCAdaptiveGridCell::split(int splitDirection)
{
	LCError err;

	if (level_ < 0)
	{
		return LCError("cannot split border cells");
	}

	if (splitDirection >= ranges_.size() / 2)
	{
		return LCError("split direction out of range");
	}
	isLeaf_ = false;
	splitDirection_ = splitDirection;
	splitVal_ = 0.5*(ranges_[splitDirection * 2] + ranges_[splitDirection * 2 + 1]);
	std::vector<double> ranges1 = ranges_;
	ranges1[splitDirection * 2 + 1] = splitVal_;
	std::vector<double> ranges2 = ranges_;
	ranges2[splitDirection * 2] = splitVal_;
	child1_ = new LCAdaptiveGridCell(level_ + 1, this, ranges1);
	child2_ = new LCAdaptiveGridCell(level_ + 1, this, ranges2);


	for (auto neighbor : neighbors_)
	{
		LCErrorReturn(neighbor->replaceNeighborWithChildren(this, child1_, child2_));
	}
	child1_->addNeighbor(child2_);
	child2_->addNeighbor(child1_);

	//TODO: delete homeomorphic map info

	return err;
}

LCError LCAdaptiveGridCell::addAllNeighboursTopDown()
{
	LCError err;
	if (parent_ != nullptr)
	{
		for (auto neighbor : parent_->getNeighbors())
		{
			if (checkIsNeighbor(neighbor))
			{
				neighbors_.push_back(neighbor);
			}
		}
	}
	if (!isLeaf_)
	{
		for (auto neighbor : neighbors_)
		{
			LCErrorReturn(neighbor->replaceNeighborWithChildren(this, child1_, child2_));
		}
		child1_->addNeighbor(child2_);
		child2_->addNeighbor(child1_);
		child1_->addAllNeighboursTopDown();
		child2_->addAllNeighboursTopDown();
	}
	return err;
}



bool LCAdaptiveGridCell::cointainsPoint(Eigen::VectorXd &point)
{
	for (int i = 0; i < (int)ranges_.size() / 2; i++)
	{
		double minVal = ranges_[i * 2] - CONTAINING_EPSILON;
		double maxVal = ranges_[i * 2 + 1] + CONTAINING_EPSILON;
		if ((point[i] < minVal) || (point[i] > maxVal))
		{
			return false;
		}
	}
	return true;
}



LCAdaptiveGridCell* LCAdaptiveGridCell::getChild1()
{
	return child1_;
}
LCAdaptiveGridCell* LCAdaptiveGridCell::getChild2()
{
	return child2_;
}

int LCAdaptiveGridCell::getSplitDirection()
{
	return splitDirection_;
}

double LCAdaptiveGridCell::getSplitValue()
{
	return splitVal_;
}
double LCAdaptiveGridCell::size()
{
	double size = 1;
	for (int i = 0; i < ranges_.size() / 2; i++)
	{
		size *= ranges_[2 * i + 1] - ranges_[2 * i];
	}
	return size;
}

void LCAdaptiveGridCell::getCenter(std::vector<double> * center)
{
	for (int i = 0; i < ranges_.size() / 2; i++)
	{
		center->push_back(0.5 *(ranges_[2 * i + 1] + ranges_[2 * i]));
	}
}

void LCAdaptiveGridCell::getCellSize(Eigen::VectorXd *result)
{
	int nParams = ranges_.size() / 2;
	Eigen::VectorXd cellSize(nParams);
	for (int i = 0; i < nParams; i++)
	{
		cellSize[i] = ranges_[2 * i + 1] - ranges_[2 * i];
	}
	*result = cellSize;
}


double LCAdaptiveGridCell::getDirCellSize(int dir)
{
	int nParams = ranges_.size() / 2;
	double result = 1;
	for (int i = 0; i < nParams; i++)
	{
		if (i != dir)
		{
			result *= (ranges_[2 * i + 1] - ranges_[2 * i]);

		}
	}
	return result;
}

void LCAdaptiveGridCell::getAllLeafCells(std::vector<LCAdaptiveGridCell*> *leafCells)
{
	if (this->isLeaf_)
	{
		leafCells->push_back(this);
	}
	else
	{
		child1_->getAllLeafCells(leafCells);
		child2_->getAllLeafCells(leafCells);
	}
}

LCError LCAdaptiveGridCell::getCointainingLeaf(const std::vector<double> &params, LCAdaptiveGridCell **leafCell)
{
	/*std::cout << "here in the get containing leaf function" << std::endl;
	for (int i = 0; i < params.size(); i++)
	{
	std::cout << params[i] << " ";
	}
	std::cout << std::endl;*/
	if (this->isLeaf_)
	{
		*leafCell = this;
		return evalPametersAreInCell(params);
	}
	if (params[this->splitDirection_] < this->splitVal_)
	{
		return child1_->getCointainingLeaf(params, leafCell);
	}
	return child2_->getCointainingLeaf(params, leafCell);
}

LCError LCAdaptiveGridCell::getSmallestCointainingLeaf(const std::vector<double> &params, LCAdaptiveGridCell **leafCell)
{
	LCError err;
	if (this->isLeaf_)
	{
		*leafCell = this;
		return evalPametersAreInCell(params);
	}

	//check if they are in both sides of the split

	if (abs(params[this->splitDirection_] - this->splitVal_) < LCMathHelper::EPSILON)
	{
		LCAdaptiveGridCell *cell1, *cell2;
		LCErrorReturn(child1_->getSmallestCointainingLeaf(params, &cell1));
		LCErrorReturn(child2_->getSmallestCointainingLeaf(params, &cell2));
		double size1 = cell1->getDirCellSize(this->splitDirection_);
		double size2 = cell2->getDirCellSize(this->splitDirection_);
		if (size1 < size2)
		{
			*leafCell = cell1;
			return err;
		}
		*leafCell = cell2;
		return err;
	}



	if (params[this->splitDirection_] < this->splitVal_)
	{
		return child1_->getSmallestCointainingLeaf(params, leafCell);
	}
	return child2_->getSmallestCointainingLeaf(params, leafCell);
}

LCError LCAdaptiveGridCell::getAdjacentLeaves(const std::vector<double> &params, std::unordered_set<LCAdaptiveGridCell*> *result)
{

	if ((std::abs(params[0] - 0.937509) < 0.000001) && (std::abs(params[1] - 0.875003) < 0.000001))
	{
		std::cout << "bad basis center get adjacent leaves " << std::endl;
		std::cout << "basisCenter: ";
		for (int i = 0; i < params.size(); i++)
		{
			std::cout << params[i] << " ";
		};
		std::cout << std::endl;
		this->log(0);
		std::cout << "this split direction " << this->splitDirection_ << std::endl;
		std::cout << "this leaf: " << this->isLeaf_ << std::endl;
		std::cout << "this split val " << this->splitVal_ << std::endl;
		std::cout << "params at the split direction " << params[this->splitDirection_] << std::endl;
		if (params[this->splitDirection_] < this->splitVal_ + CONTAINING_EPSILON)
		{
			std::cout << " params smaller than splitVal. Params: " << params[this->splitDirection_] << " splitVal " << this->splitVal_ + CONTAINING_EPSILON << std::endl;
		}
		if (params[this->splitDirection_] > this->splitVal_ - CONTAINING_EPSILON)
		{
			std::cout << " params larger than splitVal. Params: " << params[this->splitDirection_] << " splitVal " << this->splitVal_ - CONTAINING_EPSILON << std::endl;
		}
		//std::cout << "params in cell " << err.isOK() << std::endl;
	}

	LCError err;
	if (this->isLeaf_)
	{
		err = evalPametersAreInCell(params);

		if (err.isOK())
		{
			result->insert(this);
		}

		/*std::cout << "results ";
		for (auto *cell : *result)
		{
		cell->log(0);
		}
		std::cout << std::endl;*/
		return err;
	}
	if (params[this->splitDirection_] <= this->splitVal_ + CONTAINING_EPSILON)
	{
		//	std::cout << " params smaller than splitVal. Params: " << params[this->splitDirection_] << " splitVal " << this->splitVal_ + CONTAINING_EPSILON << std::endl;
		LCErrorReturn(child1_->getAdjacentLeaves(params, result));
	}
	if (params[this->splitDirection_] >= this->splitVal_ - CONTAINING_EPSILON)
	{
		//	std::cout << " params larger than splitVal. Params: " << params[this->splitDirection_] << " splitVal " << this->splitVal_ - CONTAINING_EPSILON << std::endl;
		LCErrorReturn(child2_->getAdjacentLeaves(params, result));
	}

	//std::cout << "no cases hit " << std::endl;
	return err;
}


LCError LCAdaptiveGridCell::getDoubleAdjacentLeaves(const std::vector<double> &params, std::unordered_set<LCAdaptiveGridCell*> *result)
{
	return LCError("not implemented");
}



LCError LCAdaptiveGridCell::getNeighborsOrthogonalToSplit(int splitDirection, std::vector<LCAdaptiveGridCell*> *result)
{
	double minSplit = ranges_[2 * splitDirection];
	double maxSplit = ranges_[2 * splitDirection + 1];

	for (auto neigbor : neighbors_)
	{
		std::vector<double> center;
		neigbor->getCenter(&center);
		double neighborCenterSplit = center[splitDirection];
		if (((neighborCenterSplit - minSplit) > -LCMathHelper::EPSILON) && ((neighborCenterSplit - maxSplit) < LCMathHelper::EPSILON))
		{
			result->push_back(neigbor);
		}
	}

	return LCError();
}





LCError LCAdaptiveGridCell::getIntersection(const std::vector<double> &ranges1, const std::vector<double> &ranges2, std::vector<double> * result)
{
	int nDims = ranges1.size() / 2;
	if (ranges2.size() / 2 != nDims)
	{
		return LCError("cannot compare ranges of different sizes");
	}
	for (int i = 0; i < nDims; i++)
	{
		double minDim = std::max(ranges1[2 * i], ranges2[2 * i]);
		double maxDim = std::min(ranges1[2 * i + 1], ranges2[2 * i + 1]);
		result->push_back(minDim);
		result->push_back(maxDim);
		if (maxDim < minDim - LCMathHelper::EPSILON)
		{
			return LCError("empty intersection");
		}
	}
	return LCError();
}




LCError LCAdaptiveGridCell::createHomeomorphicMap(LCFunctionValue *centerShapeInfo, std::vector<LCSample*> samples)
{
	LCError err;
	centerShapeInfo_ = centerShapeInfo;

	std::vector<double> center;
	getCenter(&center);
	Eigen::VectorXd centerVector;
	LCMathHelper::std2EigenVector(center, &centerVector);

	for (auto sample : samples)
	{
		//std::cout << "before mapping homeo" << std::endl;
		//system("pause");
		LCFunctionValue* shapeInfo = nullptr;

		//std::cout << "computing homeomorphic mapping for center shape at " << centerVector.transpose() << " to shape " << sample->getCenter().transpose() << std::endl;
		err = (sample->getShapeInfo()->getHomeophicShapeInfo(centerShapeInfo_, &shapeInfo));
		if (!err.isOK())
		{
			std::cout << err.internalDescription() << std::endl;
		}

		//map the physics - 7-11-18 removed by czw because contained references to tetCAD	
		/*OCTetCADShapeInfo* sampleShapeInfo = dynamic_cast<OCTetCADShapeInfo*>(sample->getShapeInfo());
		if (sampleShapeInfo->isTetCAD())
		{
		std::cout << "cast to tet CAD in create homeomorphic map "<<std::endl;

		if (sampleShapeInfo->nPhysics() > 0)
		{
		std::cout << "physics to compute " << std::endl;
		LCFunctionValueMap * map;
		LCFunctionValueMap::newMapFromShapePair(sample->getShapeInfo(), shapeInfo, &map);
		map->mapPrecomputedPhysics(sample->getShapeInfo(), shapeInfo);
		delete map;
		}
		}*/

		samples_[sample] = shapeInfo;

	}


	return err;
}


LCError LCAdaptiveGridCell::addSample(LCSample* sample)
{
	LCError err;
	LCFunctionValue* shapeInfo = nullptr;
	LCErrorReturn(sample->getShapeInfo()->getHomeophicShapeInfo(centerShapeInfo_, &shapeInfo));
	samples_[sample] = shapeInfo;

	//map the physics - removed by czw 7-11-18 because contained references to physics
	/*OCTetCADShapeInfo* sampleShapeInfo = dynamic_cast<OCTetCADShapeInfo*>(sample->getShapeInfo());//no more dynamic cast
	if (sampleShapeInfo)
	{
	if (sampleShapeInfo->nPhysics() > 0)
	{
	LCFunctionValueMap * map;
	LCFunctionValueMap::newMapFromShapePair(sample->getShapeInfo(), shapeInfo, &map);
	map->mapPrecomputedPhysics(sample->getShapeInfo(), shapeInfo);
	delete map;
	}
	}*/
	return err;
}

LCError LCAdaptiveGridCell::addSample(LCSample* sample, LCFunctionValue* shapeInfo)
{
	LCError err;
	if (samples_.find(sample) != samples_.end())
	{
		return LCError("cant add aleady existing sample");
	}
	samples_[sample] = shapeInfo;
	return err;
}

LCError LCAdaptiveGridCell::evalPametersAreInCell(const std::vector<double> & params)
{
	/*std::cout << "here in evalParametrsAreInCell and these are the parameters ";
	for (int i = 0; i < params.size(); i++)
	{
	std::cout << params[i] << " ";
	}
	std::cout << std::endl;*/

	if (params.size() != ranges_.size() / 2)
	{
		return LCError("invalid number of parameters");
	}
	for (int i = 0; i < ranges_.size() / 2; i++)
	{
		if ((params[i] < (ranges_[2 * i]) - LCMathHelper::EPSILON)
			|| (params[i] > (ranges_[2 * i + 1]) + LCMathHelper::EPSILON))
		{
			/*
			if ((params[i] < (ranges_[2 * i]) - LCMathHelper::EPSILON)){
			std::cout << "(params[i] < (ranges_[2 * i]) - LCMathHelper::EPSILON) " << std::endl;
			}
			if ((params[i] > (ranges_[2 * i + 1]) + LCMathHelper::EPSILON)){
			std::cout << "(params[i] > (ranges_[2 * i + 1]) + LCMathHelper::EPSILON) " << std::endl;
			}
			std::cout << "params[i] " << params[i] << " " << i << std::endl;

			std::cout << "here in evalParametrsAreInCell and found bad parameters ";
			for (int i = 0; i < ranges_.size() / 2; i++)
			{
			std::cout << params[i] << " ";
			}
			std::cout << std::endl;
			std::cout << "here in evalParametrsAreInCell and these are the ranges ";
			for (int i = 0; i < ranges_.size(); i++)
			{
			std::cout << ranges_[i] << " ";
			}
			std::cout << std::endl;
			*/
			return LCError("parameter ouside cell");
		}
	}
	return LCError();
}


LCError LCAdaptiveGridCell::getShapeInfoForSample(LCSample* sample, LCFunctionValue **shapeInfo)
{
	auto mapSample = samples_.find(sample);
	if (mapSample != samples_.end())
	{
		*shapeInfo = mapSample->second;
		return LCError();
	}

	return LCError("sample not found");
}


LCError LCAdaptiveGridCell::getShapeInfoForSampleCenter(Eigen::VectorXd &center, LCFunctionValue **shapeInfo)
{
	for (auto sample : samples_)
	{
		if ((sample.first->getCenter() - center).norm() < LCMathHelper::EPSILON)
		{
			*shapeInfo = sample.second;
			return LCError();
		}
	}

	return LCError("sample with center not found");
}

LCError LCAdaptiveGridCell::getClosestShapeInfo(Eigen::VectorXd &center, LCFunctionValue **shapeInfo)
{
	double foundMin = false;
	double minDist = std::numeric_limits<double>::max();
	for (auto sample : samples_)
	{
		double dist = (sample.first->getCenter() - center).norm();
		if (dist < minDist)
		{
			*shapeInfo = sample.second;
			foundMin = true;
			minDist = dist;
		}
	}

	if (!foundMin)
	{
		return LCError("could not find closest sample");
	}
	return LCError();

}

//czw 7-11-18 - got rid of logging function because it contains references to mesh
/*void LCAdaptiveGridCell::logAllSamples(std::string name)
{
std::cout << "outputting samples for cell at ";
logRanges(0);
for (auto s : samples_)
{
{
std::stringstream filename;
filename << "..//..//..//data//meshes//sampleslist_" << name << s.first->getCenter().transpose();
TriangularMesh tri_meshSource;
tri_meshSource.getFromTetMesh((s.first->getShapeInfo())->getMesh());
OCMorphingUtils::VisualizationUtils::writeMeshToPLYFile(tri_meshSource, (filename.str() + "first.ply").c_str());
TriangularMesh tri_meshOrig;
tri_meshOrig.getFromTetMesh((s.first->getShapeInfo())->getMesh());
OCMorphingUtils::VisualizationUtils::writeMeshToPLYFile(tri_meshOrig, (filename.str() + "second.ply").c_str());
std::cout << "finished mapping" << std::endl;
}
}
}*/

LCError LCAdaptiveGridCell::getNeightbourShapeInfoMappingMapping(LCAdaptiveGridCell* neighbor, Eigen::VectorXd center,
	LCFunctionValueMap **result)
{
	//get the point of contact between the neightbours
	LCError err;
	double foundMin = false;
	double minDist = std::numeric_limits<double>::max();
	auto neighborSamples = neighbor->getPrecomputedSamples();
	LCFunctionValue *neighborShapeInfo = nullptr;
	LCFunctionValue *thisShapeInfo = nullptr;
	//LCFunctionValue *orgShapeInfo = nullptr;
	//LCFunctionValue *orgShapeInf2o = nullptr;
	//std::cout << "neighbor = "; neighbor->logRanges(0);
	//std::cout << "center = " << center.transpose() << std::endl;
	//std::cout << "neighbor centerd at " << center.transpose() << std::endl;
	for (auto thisSample : samples_)
	{
		auto neightborSample = neighborSamples.find(thisSample.first);
		if (neightborSample != neighborSamples.end())
		{
			double dist = (thisSample.first->getCenter() - center).norm();
			if (dist < minDist)
			{
				//std::cout << "closest sample centerd at " << thisSample.first->getCenter().transpose() << std::endl;
				//std::cout << "neighbor sample centerd at " << neightborSample->first->getCenter().transpose() << std::endl;
				thisShapeInfo = thisSample.second;
				neighborShapeInfo = neightborSample->second;
				//orgShapeInfo = thisSample.first->getShapeInfo();
				//orgShapeInf2o = neightborSample->first->getShapeInfo();
				foundMin = true;
				minDist = dist;
			}
		}
	}

	if (!foundMin)
	{
		std::cout << "not closest sample centerd at" << std::endl;
		system("pause");
		return LCError("could not find closest sample");
	}


	auto storedMap = storedMaps.find(std::make_pair(thisShapeInfo, neighborShapeInfo));
	if (storedMap != storedMaps.end())
	{
		*result = storedMap->second;
	}
	else{
		LCFunctionValueMap *shapeInfoMapping;
		LCErrorReturn(LCFunctionValueMap::newMapFromShapePair(neighborShapeInfo, thisShapeInfo, &shapeInfoMapping));
		storedMaps[std::make_pair(thisShapeInfo, neighborShapeInfo)] = shapeInfoMapping;
		*result = shapeInfoMapping;
	}


	return err;

}

LCError LCAdaptiveGridCell::getAllInfluencingSamples(LCBasisFunction::LCBasisFunctionType basisType,
	std::vector<std::pair<LCSample*, LCFunctionValue*>> *influencingSamples,
	std::vector<LCFunctionValue*> *shapesToDelete)
{

	LCError err;
	std::unordered_map<LCSample*, bool> processedSamples;
	//std::cout << " h 0" << std::endl;
	for (auto sample : samples_)
	{
		influencingSamples->push_back(sample);
		processedSamples[sample.first] = true;
	}
	//std::cout << " h 1" << std::endl;

	switch (basisType)
	{
	case LCBasisFunction::LCBasisFunctionType::LINEAR_BSPLINE:
		break;
	case LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE:
	{
		for (auto neighbor : neighbors_)
		{
			std::vector<LCSample*> neighborSamples;
			neighbor->getSamples(&neighborSamples);
			for (auto neighborSample : neighborSamples)
			{
				if (processedSamples.find(neighborSample) == processedSamples.end())
				{
					LCFunctionValueMap* shappeInfoMap;
					LCErrorReturn(getNeightbourShapeInfoMappingMapping(neighbor, neighborSample->getCenter(), &shappeInfoMap));
					processedSamples[neighborSample] = true;
					LCFunctionValue* neightborShapeInfo;
					LCErrorReturn(neighbor->getShapeInfoForSample(neighborSample, &neightborShapeInfo));
					LCFunctionValue* mappedShapeInfo; //needs to be deleted later
					LCErrorReturn(shappeInfoMap->newShapeFromMap(neightborShapeInfo, &mappedShapeInfo));
					shapesToDelete->push_back(mappedShapeInfo);
					influencingSamples->push_back(std::make_pair(neighborSample, mappedShapeInfo));
					//delete shappeInfoMap;
				}
			}
		}
	}
	break;
	default:
		std::cout << "type unknown" << std::endl;
		return LCError("basis type unknown");
		break;
	}

	return err;
}



LCError LCAdaptiveGridCell::evalShapeInfo(const std::vector<double> &params,
	LCBasisFunction::LCBasisFunctionType basisType, LCFunctionValue **result)
{


	LCError err;
	err = (evalPametersAreInCell(params));
	if (!err.isOK())
	{
		std::cout << err.internalDescription();
		return err;
	}
	Eigen::VectorXd pos;
	LCMathHelper::std2EigenVector(params, &pos);
	std::vector<LCFunctionValue*> shapesToDelete;
	std::vector<std::pair<LCSample*, LCFunctionValue*>> influencingSamples;


	err = (getAllInfluencingSamples(basisType, &influencingSamples, &shapesToDelete));

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
		this->logRanges(0);
		std::cout << "params " << params[0] << ", " << params[1] << std::endl;
		std::cout << "weights do not sum to one: " << weightSum << std::endl;
		//system("pause"); 
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

LCError LCAdaptiveGridCell::evalDeriv(const std::vector<double> &params,
	LCBasisFunction::LCBasisFunctionType basisType, int direction, LCFunctionDeriv **result)
{

	LCError err;
	LCErrorReturn(evalPametersAreInCell(params));
	Eigen::VectorXd pos;
	LCMathHelper::std2EigenVector(params, &pos);
	std::vector<LCFunctionValue*> shapesToDelete;
	std::vector<std::pair<LCSample*, LCFunctionValue*>> influencingSamples;

	LCErrorReturn(getAllInfluencingSamples(basisType, &influencingSamples, &shapesToDelete));
	std::vector<std::pair<double, LCFunctionValue*>> weightedSamples;
	double weightSum = 0;
	for (auto sample : influencingSamples)
	{
		double weight;
		LCErrorReturn(sample.first->getDerivInterpolationWeight(pos, direction, &weight));
		weightSum += weight;
		weightedSamples.push_back(std::make_pair(weight, sample.second));
	}

	if (abs(weightSum) > LCMathHelper::EPSILON)
	{
		//return LCError("deriv weights do not sum to zero");
	}
	LCErrorReturn(LCFunctionDeriv::newFromWeightedSum(weightedSamples, result));
	for (auto shape : shapesToDelete)
	{
		delete shape;
	}
	return err;

}



LCError LCAdaptiveGridCell::getMidValues(int splitDirection, std::vector<Eigen::VectorXd> *midValues)
{
	LCError err;
	int nParams = ranges_.size() / 2;

	if (nParams == 1)
	{
		Eigen::VectorXd midValue(1);
		midValue(0) = 0.5*(ranges_[0] + ranges_[1]);
		(*midValues).push_back(midValue);
		return err;
	}

	double midPoint = 0.5*(ranges_[2 * splitDirection] + ranges_[2 * splitDirection + 1]);
	std::vector<Eigen::VectorXi> combinations;
	LCErrorReturn(LCMathHelper::computeCombinations(nParams - 1, 2, &combinations));
	for (auto combination : combinations)
	{
		Eigen::VectorXd midValue(nParams);
		int count = 0;
		for (int i = 0; i < nParams; i++)
		{
			if (i == splitDirection)
			{
				midValue[i] = midPoint;
			}
			else
			{
				midValue[i] = (1 - combination[count]) * ranges_[2 * i] + combination[count] * ranges_[2 * i + 1];
				count++;
			}
		}
		midValues->push_back(midValue);
	}

	return err;
}


LCError LCAdaptiveGridCell::getSample(const Eigen::VectorXd &center, LCSample** result)
{
	for (auto sample : samples_)
	{
		//std::cout << (sample.first->getCenter() - center).norm() << std::endl;
		if ((sample.first->getCenter() - center).norm() < CONTAINING_EPSILON)
		{
			*result = sample.first;
			return LCError();
		}
	}
	return LCError("sample not founnd");
}


void LCAdaptiveGridCell::getSamples(std::vector<LCSample*> *samples)
{
	for (auto sample : samples_)
	{
		samples->push_back(sample.first);
	}
}


void LCAdaptiveGridCell::getDoubleAdjacentSamples(std::vector<LCSample*> *samples)
{
	std::vector<LCSample*> allSamples;
	getSamples(&allSamples);
	for (auto neighbor : neighbors_)
	{
		std::vector<double> nCenter;
		neighbor->getCenter(&nCenter);
		neighbor->getSamples(&allSamples);

	}

	//remove duplicates
	std::set<LCSample*> uniqueSamples;
	for (unsigned i = 0; i < allSamples.size(); ++i)
	{
		uniqueSamples.insert(allSamples[i]);
	}
	samples->assign(uniqueSamples.begin(), uniqueSamples.end());

}

LCFunctionValue* LCAdaptiveGridCell::getCenterShapeInfo()
{
	return centerShapeInfo_;
}


void LCAdaptiveGridCell::setCenterShapeInfo(LCFunctionValue * newShapeInfo)
{
	if (centerShapeInfo_ != nullptr)
	{
		delete centerShapeInfo_;
	}
	centerShapeInfo_ = newShapeInfo;
}

std::unordered_map<LCSample*, LCFunctionValue*> LCAdaptiveGridCell::getPrecomputedSamples()
{
	return samples_;
}


LCError LCAdaptiveGridCell::getAffectedRangeWhenSplitting(const Eigen::VectorXd &center, int dir,
	LCBasisFunction::LCBasisFunctionType basisType, std::vector<std::vector<double>> *uncoverdRegion)
{

	if (ranges_.size() / 2 < dir)
	{
		return LCError("paramters out of range");
	}

	double cellCenter = 0.5*  (ranges_[2 * dir] + ranges_[2 * dir + 1]);
	double cellMin = ranges_[2 * dir];
	double cellMax = ranges_[2 * dir + 1];
	double sampleCenter = center(dir);
	std::vector<double> halfRanges = ranges_;
	if (cellCenter > (sampleCenter - LCMathHelper::EPSILON))
	{
		halfRanges[2 * dir] = cellCenter;
	}
	if (cellCenter < (sampleCenter + LCMathHelper::EPSILON))
	{
		halfRanges[2 * dir + 1] = cellCenter;
	}
	if (basisType == LCBasisFunction::LCBasisFunctionType::LINEAR_BSPLINE)
	{
		uncoverdRegion->push_back(halfRanges);
		return LCError();
	}

	if (basisType == LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE)
	{
		if ((sampleCenter < cellMin - LCMathHelper::EPSILON) || (sampleCenter > cellMax + LCMathHelper::EPSILON))
		{
			uncoverdRegion->push_back(halfRanges);
			return LCError();
		}

		/*		std::vector<double> rootRanges = getRoot()->getRanges();
		if ( (sampleCenter > cellCenter) && (abs(cellMin - rootRanges[2 * dir]) < LCMathHelper::EPSILON))
		{
		std::vector<double> borderRanges = ranges_;
		borderRanges[2 * dir] = cellMin - 10;
		borderRanges[2 * dir + 1] = cellMin;
		uncoverdRegion->push_back(borderRanges);
		return LCError();
		}
		if ((sampleCenter < cellCenter) && (abs(cellMax - rootRanges[2 * dir + 1]) < LCMathHelper::EPSILON))
		{
		std::vector<double> borderRanges = ranges_;
		borderRanges[2 * dir] = cellMax;
		borderRanges[2 * dir + 1] = cellMax + 10;
		uncoverdRegion->push_back(borderRanges);
		return LCError();
		}
		*/
		for (auto neighbor : neighbors_)
		{
			std::vector<double> nCenter;
			neighbor->getCenter(&nCenter);
			double neighborCenter = nCenter[dir];
			if ((neighborCenter < cellMin) && (sampleCenter > cellCenter))
			{
				uncoverdRegion->push_back(neighbor->ranges_);
			}
			if ((neighborCenter > cellMax) && (sampleCenter < cellCenter))
			{
				uncoverdRegion->push_back(neighbor->ranges_);
			}
		}
		return LCError();
	}


	return LCError("unknown basis type");
}

void LCAdaptiveGridCell::log(int level)
{
	logRanges(level);
	if (!isLeaf_)
	{
		child1_->log(level + 1);
		child2_->log(level + 1);
	}
}

void LCAdaptiveGridCell::logRanges(int level)
{
	for (int i = 0; i < level; i++)
	{
		std::cout << "-";
	}
	std::cout << "> [";
	for (int i = 0; i < ranges_.size(); i++)
	{
		std::cout << ranges_[i] << " ";
	}
	std::cout << "]" << std::endl;

}

void LCAdaptiveGridCell::logSamples()
{
	std::cout << "total samples : " << samples_.size() << std::endl;
	for (auto sample : samples_)
	{
		std::cout << "sample = " << sample.first->getCenter().transpose() << std::endl;
	}
}

bool LCAdaptiveGridCell::checkIsNeighbor(const LCAdaptiveGridCell * other)
{
	//std::vector<double> otherRanges = other->ranges_;
	int nParams = ranges_.size() / 2;
	for (int i = 0; i < nParams; i++)
	{
		double minA = ranges_[2 * i];
		double maxA = ranges_[2 * i + 1];
		double minB = other->ranges_[2 * i];
		double maxB = other->ranges_[2 * i + 1];
		if (((minB - maxA) > CONTAINING_EPSILON) || ((minA - maxB) > CONTAINING_EPSILON))
		{
			return false;
		}
	}
	return true;
}

LCError LCAdaptiveGridCell::replaceNeighborWithChildren(LCAdaptiveGridCell* neighbor,
	LCAdaptiveGridCell* child1, LCAdaptiveGridCell* child2)
{
	/*if (child1->ranges_[0] == -0.375 && child1->ranges_[1] == -0.3125 && child1->ranges_[2] == 0 && child1->ranges_[3] == 0.5)
	{
		std::cout << "found the child " << std::endl;
		neighbor->log(0);
		this->log(0);
		system("pause");
	}*/

	int origSize = neighbors_.size();
	neighbors_.remove(neighbor);
	if (checkIsNeighbor(child1))
	{
		neighbors_.push_back(child1);
		/*if (child1->ranges_[0] == -0.375 && child1->ranges_[1] == -0.3125 && child1->ranges_[2] == 0 && child1->ranges_[3] == 0.5)
		{
			std::cout << "found the child and added it" << std::endl;
		}*/
	}
	if (checkIsNeighbor(child2))
	{
		neighbors_.push_back(child2);
	}
	if (neighbors_.size() - origSize < 0)
	{
		return LCError("replacing neighbor with children did not work properly");
	}

	return LCError();
}
void LCAdaptiveGridCell::addNeighbor(LCAdaptiveGridCell* neighbor)
{
	//if (neighbor->ranges_[0] == -0.375 && neighbor->ranges_[1] == -0.3125 && neighbor->ranges_[2] == 0 && neighbor->ranges_[3] == 0.5)
	//{
	//	std::cout << "found the neighbor " << std::endl;
	//	neighbor->log(0);
	//	this->log(0);
	//	system("pause");
	//}
	neighbors_.push_back(neighbor);
}

LCError getCornerId(int *optimalDirection);


LCError LCAdaptiveGridCell::getCornerId(Eigen::VectorXd &point, int *cornerId)
{
	*cornerId = 0;
	for (int i = 0; i < (int)ranges_.size() / 2; i++)
	{
		double minVal = ranges_[i * 2] - CONTAINING_EPSILON;
		double maxVal = ranges_[i * 2 + 1] + CONTAINING_EPSILON;
		if ((point[i] > ranges_[i * 2 + 1] - CONTAINING_EPSILON) && (point[i] < ranges_[i * 2 + 1] + CONTAINING_EPSILON))
		{
			*cornerId += pow(2, i);
		}
		else
		{
			if ((point[i] < ranges_[i * 2] - CONTAINING_EPSILON) || (point[i] > ranges_[i * 2] + CONTAINING_EPSILON))
			{
				return LCError("not a corner");
			}
		}
	}
	return LCError();
}

LCError LCAdaptiveGridCell::getOptimalSplitDirection(int *optimalDirection)
{
	//get all corner shapes
	std::unordered_map<int, LCFunctionValue*> cornerShapes;
	for (auto sample : samples_)
	{
		int cornerID;
		if (getCornerId(sample.first->getCenter(), &cornerID).isOK())
		{
			std::vector<double> corner;
			LCMathHelper::eigen2StdVector(sample.first->getCenter(), &corner);
			cornerShapes[cornerID] = sample.second;
		}
	}

	//compute the error in each direction
	LCError err;
	std::vector<Eigen::VectorXi> combinations;
	double nDirections = ranges_.size() / 2;
	LCMathHelper::computeCombinations(nDirections - 1, 2, &combinations);
	*optimalDirection = 0;
	double maxDirectionError = 0;
	for (int i = 0; i < nDirections; i++)
	{
		double directionError = 0;
		for (auto combination : combinations)
		{
			int pos1 = 0;
			int pos2 = 0;
			int c = 0;
			for (int j = 0; j < nDirections; j++)
			{
				if (j == i)
				{
					pos1 += 0;
					pos2 += pow(2, j);
				}
				else
				{
					pos1 += combination(c) * pow(2, j);
					pos2 += combination(c) * pow(2, j);
					c++;
				}
			}

			auto shapeRef1 = cornerShapes.find(pos1);
			auto shapeRef2 = cornerShapes.find(pos2);
			if ((shapeRef1 == cornerShapes.end()) || (shapeRef2 == cornerShapes.end()))
			{
				return LCError("could not find corner shape");
			}
			double diff;
			LCErrorReturn(shapeRef1->second->computeDifference(shapeRef2->second, &diff));
			directionError += diff;
		}
		if (directionError > maxDirectionError)
		{
			maxDirectionError = directionError;
			*optimalDirection = i;
		}
	}

	return err;

}


LCError LCAdaptiveGridCell::computeReplacingWeights(const Eigen::VectorXd &basisCenter,
	std::vector<std::pair<LCSample*, double>> *replacingWeights)
{
	LCError err;
	int nParams = ranges_.size() / 2;
	Eigen::VectorXd minRange(nParams);
	Eigen::VectorXd maxRange(nParams);
	for (int i = 0; i < nParams; i++)
	{
		minRange(i) = ranges_[2 * i];
		maxRange(i) = ranges_[2 * i + 1];
	}

	Eigen::VectorXd alphaVec = (basisCenter - minRange).cwiseQuotient(maxRange - minRange);
	std::vector<Eigen::VectorXi> combinations;
	LCMathHelper::computeCombinations(nParams, 2, &combinations);
	for (auto combination : combinations)
	{
		double weight = 1;
		for (int j = 0; j < nParams; j++)
		{
			if (combination(j) == 1)
			{
				weight *= alphaVec(j);
			}
			else
			{
				weight *= (1 - alphaVec(j));
			}
		}
		if (weight > 0)
		{

			Eigen::VectorXd sampleCenter = combination.cast<double>().cwiseProduct(maxRange)
				+ (Eigen::VectorXd::Ones(nParams) - combination.cast<double>()).cwiseProduct(minRange);
			LCSample *sample;
			LCErrorReturn(getSample(sampleCenter, &sample));
			replacingWeights->push_back(std::make_pair(sample, weight));
		}
	}
	return err;
}


LCError LCAdaptiveGridCell::getApproximationError(Eigen::VectorXd *approxError)
{

	LCError err;
	std::vector<double> params;
	getCenter(&params);
	LCFunctionValue * approximation;
	LCErrorReturn(evalShapeInfo(params, getBasisType(), &approximation));
	LCErrorReturn(approximation->computeErrorVector(centerShapeInfo_, approxError));
	delete approximation;

	return err;
}




LCBasisFunction::LCBasisFunctionType LCAdaptiveGridCell::getBasisType()
{
	LCBasisFunction::LCBasisFunctionType type;
	for (auto sample : samples_)
	{
		if ((sample.first->getBasisType(&type)).isOK())
		{
			return type;
		}
	}
	std::cout << "this cell has no basis function - this is an error";
	this->logRanges(0);
	return type;
}

void LCAdaptiveGridCell::deleteAllSamples()
{
	for (auto s : samples_)
	{
		delete s.second;
	}
	samples_.clear();
}

