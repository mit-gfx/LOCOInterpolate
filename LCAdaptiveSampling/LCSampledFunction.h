#ifndef _OC_PRECOMPUTED_PARAMETRIC_SHAPE_
#define _OC_PRECOMPUTED_PARAMETRIC_SHAPE_

#include <Eigen/Dense>
#include "LCFunction.h"
#include "LCBasisFunction.h"

class LCAdaptiveGridCell;
class LCSample;

class LCSampledFunction : public LCFunction {
public:

	LCSampledFunction(LCAdaptiveGridCell *root, std::vector<double> ranges, 
		std::vector<LCSample*> samples,
		std::vector<LCAdaptiveGridCell*> borderCells);
	LCError addBorderCells();
	int getNParams() const;
	double getMinRange(int iParam) const;
	double getMaxRange(int iParam) const;
	LCError evalDeriv(const std::vector<double> &params, int direction, LCFunctionDeriv **result);
	LCError evalShapeInfo(const std::vector<double> &params, LCFunctionValue **result);


	LCAdaptiveGridCell* getRoot();
	std::vector<LCSample*> getSamples();
	std::vector<LCAdaptiveGridCell*> getBorderCells();

private:
	LCAdaptiveGridCell *root_;
	std::vector<double> ranges_;
	std::vector<LCSample*> samples_;
	LCBasisFunction::LCBasisFunctionType basisType_;
	std::vector<LCAdaptiveGridCell*> borderCells_;

};

#endif

