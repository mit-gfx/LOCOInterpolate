#ifndef _LC_ADAPTIVE_SAMPLING_PARAMS_
#define _LC_ADAPTIVE_SAMPLING_PARAMS_

#include <Eigen/Dense>
#include <vector> 
#include <unordered_set>

struct LCAdaptiveSamplingParams
{
	enum LCSplitPolicy{ ERROR_BASED, CYCLIC };//Determines how cells of the adaptive grid will be split. 
											  //ERROR_BASED splits in the direction of largest error. 
											  //CYCLIC splits the cell in a direction orthagonal to the direction the parent was split.
	enum LCErrorMetric{ CENTER_BASED, RANDOM_SAMPLES };//Determines when cells are split. 
													   //CENTER_BASED splits the cell whenever the error at the center of the cell is greater than the threshold
													   //RANDOM_SAMPLES splits the cell whenever the error at any of n random samples in the cell is greater than the threshold

	int maxTreeDepth_; //the maximum depth of the adaptive grid. In other words, the height of the k-d tree?
	double threshold_;//the tolerated error within a cell. See errorMetric_.
	int bootstrapGridSize_;//the number of uniform cell divisions that will occur BEFORE adaptive sampling begins.
	int nErrorSamples_;//the number of random samples that will be taken. See errorMetric_.
	LCSplitPolicy splitPolicy_;//determines the split direction when cell division occurs
	LCErrorMetric errorMetric_;//determines what triggers a cell division.

	LCAdaptiveSamplingParams(int maxTreeDepth = 0, double allowedError = 0.5, 
		int bootstrap = 0, int errorSamples = 0, LCSplitPolicy policy = LCSplitPolicy::ERROR_BASED, LCErrorMetric errorMetric=LCErrorMetric::CENTER_BASED )
	{
		std::cout << "lc adaptive sampling constructor" << std::endl;
		splitPolicy_ = policy;
		maxTreeDepth_ = maxTreeDepth;
		threshold_ = allowedError;
		bootstrapGridSize_ = bootstrap;
		nErrorSamples_ = errorSamples;
		errorMetric_ = errorMetric;
	}

	void log()
	{
		std::cout << "split policy " << splitPolicy_ << " maxTreeDepth " << maxTreeDepth_ << " threshold " << threshold_ << " bootstrapGridSize_ " << bootstrapGridSize_ << " nErrorSamples " << nErrorSamples_ << " errorMetric " << errorMetric_ << std::endl;
	}
};

#endif 
