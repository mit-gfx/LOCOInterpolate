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

	/* removed by czw - 7-11-2018 
	if (dynamic_cast<OCTetCADShapeInfo*>(weightedSamples[0].second))
	{
		auto tetShapeInfo = dynamic_cast<OCTetCADShapeInfo*>(weightedSamples[0].second);

		OCTetCADShapeInfo* shapeInfo = new OCTetCADShapeInfo(tetShapeInfo, weightedSamples[0].first);

		if (tetShapeInfo->precomputedPhysics.size() > 3) {
			std::cout << "size = " << tetShapeInfo->precomputedPhysics.size() << " ";
			for (int i = 0; i < tetShapeInfo->precomputedPhysics.size(); ++i)
				std::cout << tetShapeInfo->precomputedPhysics[i]->values[0] << " ";
			std::cout << std::endl;
		}
		for (int i = 1; i < weightedSamples.size(); i++)
		{
			tetShapeInfo = dynamic_cast<OCTetCADShapeInfo*>(weightedSamples[i].second);
			LCErrorReturn(shapeInfo->add(weightedSamples[i].first, tetShapeInfo));
		}
		*result = shapeInfo;
		return err;
	}
	*/

	return LCError("Shape info type not specified");
}

