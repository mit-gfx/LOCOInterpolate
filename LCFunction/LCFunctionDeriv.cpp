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

	return LCError("Shape info type not specified");

}
