#include <iostream>
#include "LCFunctionValue.h"
#include "LCError.h"
#include "time.h"
#include "LCRealFunctionValue.h"

LCFunctionValue::LCFunctionValue()
{
}

LCError LCFunctionValue::newFromWeightedSum(const std::vector<std::pair<double, LCFunctionValue*>> weightedSamples, LCFunctionValue **result)
{
	LCError err;
	if (weightedSamples.size() == 0)
	{
		return LCError("cannot create new from empty list of weighted samples");
	}


	if (dynamic_cast<LCRealFunctionValue*>(weightedSamples[0].second))
	{
		LCRealFunctionValue* shapeInfo = new LCRealFunctionValue(0);
		for (auto weightedSample : weightedSamples)
		{
			LCErrorReturn(shapeInfo->add(weightedSample.first, weightedSample.second));
		}
		*result = shapeInfo;
		return err;
	}

	return LCError("Shape info type not specified");
}

