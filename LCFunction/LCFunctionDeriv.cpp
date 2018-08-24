#include <iostream>
#include "LCFunctionDeriv.h"
#include "LCFunctionValue.h"
#include "LCRealFunctionValue.h"
#include "LCError.h"
#include "LCRealFunctionDeriv.h"

LCFunctionDeriv::LCFunctionDeriv()
{
}

LCError LCFunctionDeriv::newFromWeightedSum(const std::vector<std::pair<double, LCFunctionValue*>> weightedSamples, LCFunctionDeriv **result)
{

	LCError err;
	if (weightedSamples.size() == 0)
	{
		return LCError("cannot create new from empty list of weighted samples");
	}


	if (dynamic_cast<LCRealFunctionValue*>(weightedSamples[0].second))
	{
		LCRealFunctionDeriv* derivInfo = new LCRealFunctionDeriv(0);
		for (auto weightedSample : weightedSamples)
		{
			LCErrorReturn(derivInfo->add(weightedSample.first, weightedSample.second));
		}
		*result = derivInfo;
		return err;
	}


	//removed by czw 07-11-2018
	/*if (dynamic_cast<OCTetCADShapeInfo*>(weightedSamples[0].second))
	{
		auto tetShapeInfo = dynamic_cast<OCTetCADShapeInfo*>(weightedSamples[0].second);
		OCTetCADDerivInfo* derivInfo = new OCTetCADDerivInfo(tetShapeInfo, weightedSamples[0].first);
		for (int i = 1; i < weightedSamples.size(); i++)
		{
			tetShapeInfo = dynamic_cast<OCTetCADShapeInfo*>(weightedSamples[i].second);
			LCErrorReturn(derivInfo->add(weightedSamples[i].first, tetShapeInfo));
		}

		*result = derivInfo;
		return err;
	}*/

	return LCError("Shape info type not specified");

}
