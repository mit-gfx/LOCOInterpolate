#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <iostream>
#include "LCProtoConverter.h"
#include "LCAdaptiveGridCell.h"
#include "LCSampledFunction.h"
#include "LCSample.h"
#include "LCRealFunctionValue.h"
#include "LCFunctionValue.h"
#include "LCBasisFunction.h"
#include "LCError.h"

LCError PrecomputedParametricShapeConverter::outputToProto(std::string filename, bool useAscii, LCSampledFunction* shape)
{
	LCError err;
	proto::PrecomputedParametricShape  shapeProto;
	size_t lastindex = filename.find_last_of(".");
	std::string filepath = filename.substr(0, lastindex);
	std::cout << "the filepath is " << filepath << std::endl;
	std::cout << "the filename is " << filename << std::endl;

	LCErrorReturn(convertToProto(shape, filepath, &shapeProto));
	LCErrorReturn(saveToFile(shapeProto, filename.c_str(), useAscii));
	return err; 
}

LCError PrecomputedParametricShapeConverter::loadFromProto(std::string filename, bool useAscii, LCSampledFunction** shape)
{
	LCError err;
	proto::PrecomputedParametricShape shapeProto;
	LCErrorReturn(loadFromFile(filename.c_str(), useAscii, &shapeProto));
	LCErrorReturn(convertToCplusplus(shapeProto, shape));
	return err; 
}

LCError convertBasisFuntionToProto(LCBasisFunction *basisFuntion, proto::BasisFunction* result)
{
	LCLinearBSpline* linearBSpline = dynamic_cast<LCLinearBSpline*>(basisFuntion);
	if (linearBSpline != nullptr)
	{
		proto::LinearBSpline* protoBasisFunction = new proto::LinearBSpline();
		Eigen::VectorXd center = linearBSpline->getCenter();
		for (int i = 0; i < center.size(); i++)
		{
			protoBasisFunction->add_center(center(i));
		}
		Eigen::VectorXd support = linearBSpline->getSupport();
		for (int i = 0; i < support.size(); i++)
		{
			protoBasisFunction->add_support(support(i));
		}
		protoBasisFunction->set_weight(linearBSpline->getWeight());
		result->set_allocated_linearbspline(protoBasisFunction);
		return LCError();
	}
	LCCubicBSpline* cubicBSpline = dynamic_cast<LCCubicBSpline*>(basisFuntion);
	if (cubicBSpline != nullptr)
	{
		proto::CubicBSpline* protoBasisFunction = new proto::CubicBSpline();
		Eigen::VectorXd center = cubicBSpline->getCenter();
		for (int i = 0; i < center.size(); i++)
		{
			protoBasisFunction->add_center(center(i));
		}
		Eigen::VectorXd support = cubicBSpline->getSupport();
		for (int i = 0; i < support.size(); i++)
		{
			protoBasisFunction->add_support(support(i));
		}
		protoBasisFunction->set_weight(cubicBSpline->getWeight());
		result->set_allocated_cubicbspline(protoBasisFunction);
		return LCError();
	}
	return LCError("undefined type of basis function");
}


LCError convertShapeInfoToProto(LCFunctionValue * shapeInfo, std::string filepath, int* meshCounter, proto::ShapeInfo* result)
{
	LCRealFunctionValue* funcShapeInfo = dynamic_cast<LCRealFunctionValue*>(shapeInfo);
	if (funcShapeInfo != nullptr)
	{
		proto::FunctionTestShapeInfo* protoShapeInfo = new proto::FunctionTestShapeInfo();
		protoShapeInfo->set_val(funcShapeInfo->val);
		result->set_allocated_functionshapeinfo(protoShapeInfo);
		return LCError();
	}
	/* removed by czw 7-11-2018
	OCTetCADShapeInfo* tetShapeInfo = dynamic_cast<OCTetCADShapeInfo*>(shapeInfo);
	if (tetShapeInfo != nullptr)
	{
		proto::TetMesh* protoTetMesh = new proto::TetMesh();
		std::stringstream filename;
		filename << filepath << "_mesh" << *meshCounter << ".proto";
		*meshCounter = *meshCounter + 1; 
		protoTetMesh->set_filename(filename.str());
		LCErrorReturn(OCTetCADShapeInfoConverter::outputToProto(filename.str(), tetShapeInfo));
		result->set_allocated_tetmesh(protoTetMesh);
		return LCError();
	}
	*/
	return LCError("undefined type of shape info");
}

LCError convertetPrecomputedSampleToProto(LCSample *sample, int id, std::string filepath, int * meshCounter, proto::PrecomputedSample* result)
{
	LCError err;
	result->set_id(id);
	Eigen::VectorXd center = sample->getCenter();
	for (int i = 0; i < center.size(); i++)
	{
		result->add_center(center(i));
	}
	for (auto basisFunction : sample->getBasisFunctions())
	{
		convertBasisFuntionToProto(basisFunction.second, result->mutable_basisfuntions()->Add());
	}
	LCErrorReturn(convertShapeInfoToProto(sample->getShapeInfo(),filepath, meshCounter, result->mutable_shapeinfo()));
	return err;
}

LCError convertetAdaptiveGridCellToProto(LCAdaptiveGridCell * cell, std::string filepath, int * meshCounter, std::unordered_map<LCSample*, int> &samplesToIds,
	proto::AdaptiveGridCell* result)
{
	LCError err;
	for (double range : cell->ranges_)
	{
		result->add_ranges(range); 
	}
	if (cell->isLeaf())
	{
		proto::AdaptiveGridLeaf* leaf = new proto::AdaptiveGridLeaf();
		LCErrorReturn(convertShapeInfoToProto(cell->getCenterShapeInfo(), filepath, meshCounter, leaf->mutable_centershapeinfo()));
		for (auto sample : cell->getPrecomputedSamples())
		{
			int sampleId = samplesToIds[sample.first];
			proto::HomeomorphicSample *protoSample = new proto::HomeomorphicSample();
			LCErrorReturn(convertShapeInfoToProto(sample.second, filepath, meshCounter, protoSample->mutable_shapeinfo()));
			protoSample->set_precomputedsampleid(sampleId);
			leaf->mutable_homeomorphicsamples()->AddAllocated(protoSample);
		}
		result->set_allocated_leaf(leaf);
		return err;
	}
	proto::AdaptiveGridInterNode* interNode = new proto::AdaptiveGridInterNode();
	LCErrorReturn(convertetAdaptiveGridCellToProto(cell->getChild1(), filepath, meshCounter, samplesToIds, interNode->mutable_child1()));
	LCErrorReturn(convertetAdaptiveGridCellToProto(cell->getChild2(), filepath, meshCounter, samplesToIds, interNode->mutable_child2()));
	interNode->set_spliddirection(cell->getSplitDirection());
	interNode->set_splitval(cell->getSplitValue());
	result->set_allocated_internode(interNode);

	return err;

}
LCError PrecomputedParametricShapeConverter::convertToProto(LCSampledFunction *paramShape, std::string filepath,
	proto::PrecomputedParametricShape* result) {
	int meshCounter = 0;

	std::vector<LCSample*> samples = paramShape->getSamples();
	std::unordered_map<LCSample*, int> samplesToIds;
	for (int i = 0; i < samples.size(); i++)
	{
		proto::PrecomputedSample *protoSample = new proto::PrecomputedSample();
		LCErrorReturn(convertetPrecomputedSampleToProto(samples[i], i, filepath, &meshCounter, protoSample));
		result->mutable_samples()->AddAllocated(protoSample);
		samplesToIds[samples[i]] = i;
	}
	for (auto borderCell : paramShape->getBorderCells())
	{
		convertetAdaptiveGridCellToProto(borderCell, filepath,  &meshCounter, samplesToIds, result->add_bordercells());
	}
	LCErrorReturn(convertetAdaptiveGridCellToProto(paramShape->getRoot(),  filepath, &meshCounter, samplesToIds, result->mutable_root()));

	std::vector<double> ranges;
	paramShape->getRanges(&ranges);
	for (int i = 0; i < ranges.size(); i++)
	{
		result->add_range(ranges[i]);
	}

	return LCError();
}

//----------------------------------------------------------------------------------------------------------

LCError convertBasisFuntionToCpp(const proto::BasisFunction &basisFunction, LCBasisFunction **result)
{
	if (basisFunction.has_linearbspline())
	{
		auto linearBSpline = basisFunction.linearbspline(); 
		Eigen::VectorXd center(linearBSpline.center_size());
		for (int i = 0; i < linearBSpline.center_size(); i++)
		{
			center(i) = linearBSpline.center(i);
		}
		Eigen::VectorXd support(linearBSpline.support_size());
		for (int i = 0; i < linearBSpline.support_size(); i++)
		{
			support(i) = linearBSpline.support(i);
		}
		double weight = linearBSpline.weight();
		*result = new LCLinearBSpline(center, support, weight);
		return LCError();
	}
	if (basisFunction.has_cubicbspline())
	{
		auto cubicBSpline = basisFunction.cubicbspline();
		Eigen::VectorXd center(cubicBSpline.center_size());
		for (int i = 0; i < cubicBSpline.center_size(); i++)
		{
			center(i) = cubicBSpline.center(i);
		}
		Eigen::VectorXd support(cubicBSpline.support_size());
		for (int i = 0; i < cubicBSpline.support_size(); i++)
		{
			support(i) = cubicBSpline.support(i);
		}
		double weight = cubicBSpline.weight();
		*result = new LCCubicBSpline(center, support, weight);
		return LCError();
	}
	return LCError("type of basis function not recognized");
}

/* removed by czw 7-11-2018
LCError convertTetMeshToCpp(const proto::TetMesh &tetMesh, OCTetCADShapeInfo * result)
{
	for (int i = 0; i < tetMesh.vertices_size(); i++)
	{
		proto::Vector3d protoVer = tetMesh.vertices(i);
		Eigen::Vector3d vertex(protoVer.x(), protoVer.y(), protoVer.z());
		result->tetMesh.vertices.push_back(vertex);
	}
	for (int i = 0; i < tetMesh.tets_size(); i++)
	{
		proto::Vector4i protoVer = tetMesh.tets(i);
		Eigen::Vector4i vertex(protoVer.x(), protoVer.y(), protoVer.z(), protoVer.w());
		result->tetMesh.tets.push_back(vertex);
	}

	for (int i = 0; i < tetMesh.boundaryconditions_size(); i++)
	{
		OCBoundaryCondition * boundaryCondition = new OCBoundaryCondition();
		auto protoBoundary = tetMesh.boundaryconditions(i);
		proto::Vector3d force = protoBoundary.force();
		boundaryCondition->value.force = Eigen::Vector3d(force.x(), force.y(), force.z());

		if (protoBoundary.isfixed())
		{
			boundaryCondition->value.type = OCBoundaryValue::OCBoundaryType::FIXED_BOUNDARY;
		}
		else
		{
			boundaryCondition->value.type = OCBoundaryValue::OCBoundaryType::BOUNDARY_FORCE;
		}
		for (int j = 0;  j < protoBoundary.vertices_size(); j++)
		{
			boundaryCondition->vertices.push_back(protoBoundary.vertices(j));
		}
		result->boundaryConditions.push_back(boundaryCondition);	
	}
	for (int i = 0; i < tetMesh.precomputedphysics_size(); i++)
	{
		OCPrecomputedPhysics *precomputedPhysics = new OCPrecomputedPhysics();
		auto protoPhysics = tetMesh.precomputedphysics(i);

		precomputedPhysics->name = protoPhysics.name();
		for (int j = 0; j < protoPhysics.values_size(); j++)
		{
			precomputedPhysics->values.push_back(protoPhysics.values(j));
		}

		result->precomputedPhysics.push_back(precomputedPhysics);
	}
	return LCError();
}
*/

LCError convertShapeInfoToCpp(const proto::ShapeInfo &shapeInfo, LCFunctionValue**result)
{
	if (shapeInfo.has_functionshapeinfo())
	{
		if (!shapeInfo.functionshapeinfo().has_val())
		{
			return LCError("value for functio shape info not found");
		}
		*result = new LCRealFunctionValue(shapeInfo.functionshapeinfo().val());
		return LCError();
	}
	/*removed by czw 7-11-2018
	if (shapeInfo.has_tetmesh())
	{
		if (shapeInfo.tetmesh().has_filename())
		{
			OCTetCADShapeInfo* tetShapeInfo = new OCTetCADShapeInfo(shapeInfo.tetmesh().filename());
			*result = tetShapeInfo;
		}
		else
		{
			OCTetCADShapeInfo* tetShapeInfo = new OCTetCADShapeInfo();
			LCErrorReturn(convertTetMeshToCpp(shapeInfo.tetmesh(), tetShapeInfo));
			*result = tetShapeInfo;
		}
		return LCError();
	}
	*/
	return LCError("type fo shape info not recognized");
}

LCError convertPrecomputedSampleToCpp(const proto::PrecomputedSample &sample, int* id, LCSample** result)
{
	LCError err;
	if (!sample.has_id())
	{
		return LCError("sample does not have id");
	}
	*id = sample.id();
	Eigen::VectorXd center(sample.center_size());
	for (int i = 0; i < sample.center_size(); i++)
	{
		center(i) = sample.center(i);
	}
	LCFunctionValue * shapeInfo = nullptr;
	if (sample.has_shapeinfo())
	{
		LCErrorReturn(convertShapeInfoToCpp(sample.shapeinfo(), &shapeInfo));
	}
	std::vector<LCBasisFunction*> basisFunctions;
	for (int i = 0; i < sample.basisfuntions_size(); i++)
	{
		LCBasisFunction* basisFunction;
		LCErrorReturn(convertBasisFuntionToCpp(sample.basisfuntions(i), &basisFunction));
		basisFunctions.push_back(basisFunction);
	}
	*result = new LCSample(center, basisFunctions, shapeInfo);
	return err; 
}

LCError convertetAdaptiveGridCellToCpp(const proto::AdaptiveGridCell &cell, 
	LCAdaptiveGridCell* parent, int level, std::unordered_map<int, LCSample*> &idsToSamples,
	LCAdaptiveGridCell **result)
{
	LCError err;
	std::vector<double> ranges;
	for (int i = 0; i < cell.ranges_size(); i++)
	{
		ranges.push_back(cell.ranges(i));
	}
	if (cell.has_leaf())
	{
		auto leaf = cell.leaf();
		LCFunctionValue * centerShapeInfo = nullptr;
		LCErrorReturn(convertShapeInfoToCpp(leaf.centershapeinfo(), &centerShapeInfo));
		std::unordered_map<LCSample*, LCFunctionValue*> samples;
		for (int i = 0; i < leaf.homeomorphicsamples_size(); i++)
		{
			int sampleID = leaf.homeomorphicsamples(i).precomputedsampleid();
			auto idToSample = idsToSamples.find(sampleID);
			if (idToSample == idsToSamples.end())
			{
				return LCError("could not map the sample from id");
			}
			LCSample* sample = idsToSamples.find(sampleID)->second;
			LCFunctionValue* shapeInfo = nullptr;
			LCErrorReturn(convertShapeInfoToCpp(leaf.homeomorphicsamples(i).shapeinfo(), &shapeInfo));
			samples[sample] = shapeInfo;
		}
		*result = new LCAdaptiveGridCell(level, parent, ranges, centerShapeInfo, samples);
		return err; 
	}
	if (cell.has_internode())
	{
		auto interNode = cell.internode(); 
		LCAdaptiveGridCell *gridCell = new LCAdaptiveGridCell(level, parent, ranges, interNode.spliddirection(), interNode.splitval());
		LCAdaptiveGridCell *child1 = nullptr;
		LCAdaptiveGridCell *child2 = nullptr;
		LCErrorReturn(convertetAdaptiveGridCellToCpp(interNode.child1(), gridCell, level + 1, idsToSamples, &child1));
		LCErrorReturn(convertetAdaptiveGridCellToCpp(interNode.child2(), gridCell, level + 1, idsToSamples, &child2));
		gridCell->setChildren(child1, child2);
		*result = gridCell;
		return err;
	}

	return  LCError("not an interal node or ");
}

LCError PrecomputedParametricShapeConverter::convertToCplusplus(const proto::PrecomputedParametricShape &protoParamShape,
	LCSampledFunction ** result)
{
	LCError err;
	//step 1: read all the precomputed samples
	std::vector<LCSample*> samples;
	std::unordered_map<int, LCSample*> idsToSamples;
	for (int i = 0; i < protoParamShape.samples_size(); i++)
	{
		int id;
		LCSample* sample = nullptr;
		LCErrorReturn(convertPrecomputedSampleToCpp(protoParamShape.samples(i), &id, &sample));
		idsToSamples[id] = sample;
		samples.push_back(sample);
	}
	LCAdaptiveGridCell* root = nullptr;
	LCErrorReturn(convertetAdaptiveGridCellToCpp(protoParamShape.root(), nullptr, 0, idsToSamples, &root));
	root->addAllNeighboursTopDown(); 
	int bootstrapNSamples = 0;
	if (protoParamShape.has_bootstrapnsamples())
	{
		bootstrapNSamples = protoParamShape.bootstrapnsamples();
	}
	std::vector<LCAdaptiveGridCell*> borderCells;
	for (int i = 0; i < protoParamShape.bordercells_size(); i++)
	{
		//int id;
		LCAdaptiveGridCell* borderCell = nullptr;
		LCErrorReturn(convertetAdaptiveGridCellToCpp(protoParamShape.bordercells(i), nullptr, 0, idsToSamples, &borderCell));
		borderCells.push_back(borderCell);
	}

	std::vector<double> ranges = root->ranges_; 
	if (protoParamShape.range_size() > 0)
	{
		ranges.clear();
		for (int i = 0; i < protoParamShape.range_size(); i++)
		{
			ranges.push_back(protoParamShape.range(i));
		}
	}


	*result = new LCSampledFunction(root, ranges, samples, borderCells);
	(*result)->addBorderCells(); 
	return err; 

}
