// TestPrecomputationNoMesh.cpp : Defines the entry point for the console application.
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

//Encode from raw pixels to disk with a single function call
//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height)
{
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

//Helper function to collect the keys from a map
//Input: map from LCBasisFunctionKey to LCBasisFunction
//Output: vector containing the LCBasisFunctionKey keys in the map
void extractKeys(std::unordered_map<LCBasisFunctionKey, LCBasisFunction*> map, std::vector<LCBasisFunctionKey> *keys)
{
	for (auto kv : map)
	{
		keys->push_back(kv.first);
	}
}

//Determines whether two vectors are identical up to a threshhold difference
//Input: Two vectors and a threshhold error
//Output: Returns a !OK error if there exists an i such that the two vectors differ by more than the threshhold at index i
LCError compareVectors(std::vector<double> vec1, std::vector<double> vec2, double threshhold)
{
	if (vec1.size() != vec2.size())
	{
		return LCError("Vectors must be same size for comparison");
	}
	double error = 0;
	double maxError = 0;

	for (int i = 1; i < vec1.size(); i++)
	{
		double cellError = std::abs(vec1[i] - vec2[i]);
		error += cellError;
		maxError = std::max({ maxError, cellError });
	}
	double averageError = error / vec1.size();
	std::cout << "max error " << maxError << " averageError " << averageError << std::endl;
	return LCError();
}

//helper function that transforms an array of floating point values into an array of values that can be written to produce a colored PNG
//input: values = array of length width*width containing floating point values
//output: image = array of length 4*width*width.
void colorData(std::vector<double> *values, std::vector<unsigned char> *image, int width, bool blackWhite = false){
	double min_val = *(std::min_element(std::begin(*values), std::end(*values)));
	double max_val = *(std::max_element(std::begin(*values), std::end(*values)));
	double range = max_val - min_val;
	for (double y = 0; y < width; y += 1)
	{
		for (double x = 0; x < width; x += 1)
		{
			double val = (*values)[y*width + x];
			double scaled_val = (val - min_val) / range * 255;
			int scaled_int = (int)scaled_val;
			(*image)[4 * width * y + 4 * x + 0] = scaled_int;
			(*image)[4 * width * y + 4 * x + 1] = (255 - scaled_int);
			(*image)[4 * width * y + 4 * x + 2] = scaled_int / 2;
			(*image)[4 * width * y + 4 * x + 3] = 255;
			if (blackWhite)
			{
				if (val > 0)
				{
					(*image)[4 * width * y + 4 * x + 0] = 0;
					(*image)[4 * width * y + 4 * x + 1] = 0;
					(*image)[4 * width * y + 4 * x + 2] = 0;
					(*image)[4 * width * y + 4 * x + 3] = 255;
				}
				else
				{
					(*image)[4 * width * y + 4 * x + 0] = 255;
					(*image)[4 * width * y + 4 * x + 1] = 255;
					(*image)[4 * width * y + 4 * x + 2] = 255;
					(*image)[4 * width * y + 4 * x + 3] = 0;
				}
			}
		}
	}
}


//Writes a png file representing the sampled parametric shape. Also returns a list of the samples values.
//Input: An 2-dimensional parametric shape 
//	     Filename of png to be written
//Output: Writes a png file that represents the parametric shape 
//		  Also returns a list of the parametric shape's sample values
LCError outputMatrix(LCFunction* shape, int imageResolution, std::string inputFileName, std::vector<double> * values){	
	std::vector<double> ranges;
	shape->getRanges(&ranges);//the upper and lower bounds of the parametric shape.

	if (ranges.size() != 4)
	{
		return (LCError("Matrix output need a 2D input"));
	}

	double delta = 1.0 / imageResolution;
	int width = imageResolution + 1;
	std::vector<unsigned char> image;
	image.resize(4 * width * width);

	for (double y = 0; y < width; y += 1)
	{
		for (double x = 0; x < width; x += 1)
		{
			double i = y*delta;
			double j = x*delta;
			std::vector<double> pos(2);//coordinates that specify a point in the parametric shape
			pos[1] = ranges[0] + i*(ranges[1] - ranges[0]);
			pos[0] = ranges[2] + j*(ranges[3] - ranges[2]);
			LCFunctionValue * info;
			LCErrorReturn(shape->evalShapeInfo(pos, &info));//evaluate the sample at these coordinates
			LCRealFunctionValue* testInfo = dynamic_cast<LCRealFunctionValue*>(info); //convert to a 1-dimensional sample
			double temp = testInfo->val;
			values->push_back(temp);

		}
	}
	colorData(values, &image, width);//"color" the data
	encodeOneStep(inputFileName.c_str(), image, width, width);//write the image to a png
	return LCError();
}

LCError outputMatrixSamples(std::vector<LCSample*> *samples, std::string inputFileName, int width = 5){
	std::cout << inputFileName << std::endl;

	std::vector<double> values(width * width);

	std::vector<unsigned char> image;
	image.resize(4 * width * width);

	for (int y = 0; y < width; y++)
	{
		for (int x = 0; x < width; x++)
		{
			values[y*width + x] = 0.0;
		}
	}

	double maxCoordX = 0; double maxCoordY = 0;
	double minCoordX = 0; double minCoordY = 0;
	for (int i = 0; i < samples->size(); i++)
	{
		LCSample* sample = (*samples)[i];
		Eigen::VectorXd center = sample->getCenter();
		maxCoordX = std::max({ maxCoordX, center[0] });
		minCoordX = std::min({ minCoordX, center[0] });
		maxCoordY = std::max({ maxCoordY, center[1] });
		minCoordY = std::min({ minCoordY, center[1] });
	}

	for (int i = 0; i < samples->size(); i++)
	{
		LCSample* sample = (*samples)[i];
		Eigen::VectorXd center = sample->getCenter();
		int x = (center[0] + -1.0*minCoordX) / (maxCoordX - minCoordX) * (width - 1);
		int y = (center[1] + -1.0*minCoordY) / (maxCoordY - minCoordY) * (width - 1);
		values[y*width + x] = 255;
	}
	std::cout << "minCoord X" << minCoordX << " minCoordY " << minCoordY << " maxCoordX " << maxCoordX << " maxCoordY " << maxCoordY << std::endl;

	colorData(&values, &image, width, true);
	encodeOneStep(inputFileName.c_str(), image, width, width);
	return LCError();
}

//Test whether the result of a linear interpolation can approximate a linear function
LCError testFunctionLinearApproximation(LCAdaptiveSamplingParams params, LCBasisFunction::LCBasisFunctionType basisType, std::string outName, int width = 5)
{

	LCRealFunction* originFunction = new LCRealFunction(2, false);//A function that is RxR -> R
	std::cout << "params error samples " << params.nErrorSamples_ << std::endl;
	params.log();
	LCAdaptiveGrid grid(originFunction, basisType, &params, true);//The kd-tree data structure used to store the sampled function
	std::vector<double> ranges;
	originFunction->getRanges(&ranges);
	LCSampledFunction* approxFunction = new LCSampledFunction(grid.getRoot(), ranges, grid.getSamples(), grid.getBorderCells());//The interpolated function

	std::cout << "number of samples " << grid.getSamples().size() << std::endl;
	std::vector<LCSample*> samples = grid.getSamples();

	for (int i = 0; i < samples.size(); i++)
	{
		LCSample* sample = samples[i];
	}

	std::vector<double> exactVals;
	std::vector<double> approxVals;
	outputMatrix(originFunction, 20, outName + "_origin.png", &exactVals);
	outputMatrix(approxFunction, 20, outName + "_approx.png", &approxVals);
	outputMatrixSamples(&samples, outName + "_samples.png", width);
	compareVectors(exactVals, approxVals, 0.1);//Compare the sampled values of the original function and the approximated function
	return LCError();
}

//Test whether parametric shapes can be read and written as proto messages
LCError protoTest()
{
	std::string outFileName = "oneDimensionalTestCmake.proto";
	int nParams = 2;
	LCError err;
	LCRealFunction* paramShape = new LCRealFunction(nParams, true);//function that is RxR -> R
	LCAdaptiveSamplingParams params(0, 0.1, 0, 0);
	LCAdaptiveGrid grid(paramShape, LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE, &params, true);//Parametric shape

	/*Randomly refine cells of the kd-tree*/
	srand(1);
	int refined = 0;
	std::vector<LCAdaptiveGridCell*> affectedCells;
	for (int i = 0; i < 12; i++)
	{
		std::vector<LCAdaptiveGridCell*> leafs;
		grid.getRoot()->getAllLeafCells(&leafs);//the leaves of the kd tree

		int id = rand() % leafs.size();
		int dir = rand() % nParams;
		if (leafs[id]->getLevel() < 5)
		{
			LCErrorReturn(grid.refineCell(leafs[id], dir, &affectedCells));//refine a randomly chosen leaf
			refined++;
		}
	}

	std::vector<double> ranges;
	paramShape->getRanges(&ranges);
	LCSampledFunction* shape = new LCSampledFunction(grid.getRoot(), ranges, grid.getSamples(), grid.getBorderCells());//The shape that results from the random refinement

	err = PrecomputedParametricShapeConverter::outputToProto(outFileName, false, shape);//output shape to proto message format
	LCSampledFunction* outputPrecomputedShape;
	LCError err2 = PrecomputedParametricShapeConverter::loadFromProto(outFileName, false, &outputPrecomputedShape);//Read the outputted shape into a c++ shape object

	if (!err2.isOK())
	{
		return LCError("Could not read from proto file");
	}

	std::vector<LCAdaptiveGridCell*> leafCells;
	outputPrecomputedShape->getRoot()->getAllLeafCells(&leafCells);

	//Compare the sample values of the original shape to the samples of the shape that were read from the proto message
	std::vector<double> protoReadValues;
	std::vector<double> protoWriteValues;
	outputMatrix(outputPrecomputedShape, 10, "readFromProto.png", &protoReadValues);
	outputMatrix(shape, 10, "writtenToProto.png", &protoWriteValues);
	std::cout << "Testing stored proto sample values" << std::endl;
	LCError comparison = compareVectors(protoReadValues, protoWriteValues, 0.1);//compare the written and read sample values
	if (comparison.isOK()){
		std::cout << "proto sample values are correct" << std::endl;
	}
	else
	{
		std::cout << "proto sample values are incorrect" << std::endl;
	}

	//Compare the basis functions of the original shape to the samples of the shape that were read from the proto message
	std::vector<LCSample*> outputSamples = outputPrecomputedShape->getSamples();
	std::vector<LCSample*> inputSamples = shape->getSamples();
	int numSamples = inputSamples.size();
	bool sampleBases = true;
	for (int i = 0; i < numSamples; i++)//For each sample, find the associated basis functions and compare them
	{
		LCSample* in = inputSamples[i];
		LCSample* out = outputSamples[i];

		std::vector<LCBasisFunctionKey> inputKeys;
		extractKeys(in->getBasisFunctions(), &inputKeys);//Get the keys of the basis functions associated with this sample. The keys contain all the relevant basis function information - i.e. center, support, etc.
		std::multiset<LCBasisFunctionKey> inputKeySet(inputKeys.begin(), inputKeys.end());
	
		std::vector<LCBasisFunctionKey> outputKeys;
		extractKeys(out->getBasisFunctions(), &outputKeys);
		std::multiset<LCBasisFunctionKey> outputKeySet(outputKeys.begin(), outputKeys.end());

		bool identicalBases = inputKeySet.size() == outputKeySet.size();
		for (auto key : inputKeySet)
		{
			identicalBases = identicalBases && outputKeySet.find(key) != outputKeySet.end();
		}
		sampleBases = sampleBases && identicalBases;
	}
	if (sampleBases)
	{
		std::cout << "proto basis functions are correct" << std::endl;
	}
	else
	{
		std::cout << "proto basis functions are incorrect" << std::endl;
	}

	return err2;
}

//Test whether linear precision is satisfied 
LCError linearPrecisionTestCubic()
{
	int nParams = 3;
	LCError err;
	LCRealFunction* paramShape = new LCRealFunction(nParams, true);
	LCAdaptiveSamplingParams params(0, 0.1, 0, 0);
	LCAdaptiveGrid grid(paramShape, LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE, &params, true);//Parametric shape
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
	for (auto sampleCombination : sampleCombinations)
	{

		Eigen::VectorXd alphaVec = sampleCombination.cast<double>() / (double)nSamples;
		Eigen::VectorXd paramsVec = shapeMinRange.cwiseProduct(Eigen::VectorXd::Ones(nParams) - alphaVec) + shapeMaxRange.cwiseProduct(alphaVec);
		std::vector<double> params;
		LCMathHelper::eigen2StdVector(paramsVec, &params);
		LCFunctionValue * actual, *approx;
		paramShape->evalShapeInfo(params, &actual);
		grid.evalShapeInfo(params, &approx);
		double error = std::abs((dynamic_cast<LCRealFunctionValue*>(actual))->val - (dynamic_cast<LCRealFunctionValue*>(approx))->val);
		if (error > 0.00001)
		{
			/* if the funciton is linear this should never be called! */
			linearPrecision = false;
		}
	}

	if (linearPrecision)
	{
		std::cout << "linear precison was verified " << std::endl;
	}
	else
	{
		std::cout << "linear precison was NOT verified " << std::endl;
	}

	return err;
}

LCError linearPrecisionTestLinear()
{
	int nParams = 1;
	LCError err;
	LCRealFunction* paramShape = new LCRealFunction(nParams, true);
	LCAdaptiveSamplingParams params(0, 0.1, 0, 0);
	LCAdaptiveGrid grid(paramShape, LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE, &params, true);//Parametric shape

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
	for (auto sampleCombination : sampleCombinations)
	{

		Eigen::VectorXd alphaVec = sampleCombination.cast<double>() / (double)nSamples;
		Eigen::VectorXd paramsVec = shapeMinRange.cwiseProduct(Eigen::VectorXd::Ones(nParams) - alphaVec) + shapeMaxRange.cwiseProduct(alphaVec);
		std::vector<double> params;
		LCMathHelper::eigen2StdVector(paramsVec, &params);
		LCFunctionValue * actual, *approx;
		paramShape->evalShapeInfo(params, &actual);
		grid.evalShapeInfo(params, &approx);
		double error = std::abs((dynamic_cast<LCRealFunctionValue*>(actual))->val - (dynamic_cast<LCRealFunctionValue*>(approx))->val);
		if (error > 0.00001)
		{
			/* if the funciton is linear this should never be called! */
			linearPrecision = false;
		}
	}
	if (linearPrecision)
	{
		std::cout << "linear precison was verified " << std::endl;
	}
	else
	{
		std::cout << "linear precison was NOT verified " << std::endl;
	}

	return err;
}


int main(int argc, char* argv[])
{

	protoTest();//Test whether parametric shapes can be read and written as proto messages
	std::vector<double> thresholds = { 0.1 };
	std::vector<int> maxDepths = { 3 };
	std::vector<int> bootstraps = { 2 };
	std::vector<int> errorSamples = { 0 };

	for (int i = 0; i < 1; i++)
	{
		double threshold = thresholds[i];
		int maxDepth = maxDepths[i];
		int bootstrap = bootstraps[i];	
		int errorSample = errorSamples[i];
		bool optimalDirection = true;

		LCAdaptiveSamplingParams params(maxDepth, threshold, bootstrap, errorSample);
		std::cout << "\nthreshold " << threshold << "\nmax depth " << maxDepth 
				  << "\nbootstrap " << bootstrap << "\nerror samples " << errorSample 
				  << "\noptimal direction " << optimalDirection << std::endl;

		std::string name = "cubic" + std::string("_thresh_") + std::to_string(threshold) + std::string("_maxDepth_") + std::to_string(maxDepth) + std::string("_bootstrap_") + std::to_string(bootstrap);
		testFunctionLinearApproximation(params, LCBasisFunction::LCBasisFunctionType::CUBIC_BSPLINE, name, 100);
	}

	linearPrecisionTestCubic(); //test whether linear precision is satisfied when using a cubic interpolation
	linearPrecisionTestLinear();//test whether linear precision is satisfied when using a linear interpolation

	system("pause");
}
