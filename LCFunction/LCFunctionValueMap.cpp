#include <iostream>
#include "LCFunctionValueMap.h"
#include "LCRealFunctionValueMap.h"
#include "LCRealFunctionValue.h"
#include "LCError.h"
#include <ctime>

LCFunctionValueMap::LCFunctionValueMap()
{
}

LCError LCFunctionValueMap::newMapFromShapePair(LCFunctionValue *sourceShape, LCFunctionValue *targetShape, LCFunctionValueMap **result)
{
	LCError err;

	if (dynamic_cast<LCRealFunctionValue*>(sourceShape))
	{
		*result = new LCRealFunctionValueMap();
		return err;
	}
/*	auto sourceTet = dynamic_cast<OCTetCADShapeInfo*>(sourceShape);

	if (sourceTet)
	{
		auto targetTet = dynamic_cast<OCTetCADShapeInfo*>(targetShape);
		if (!targetTet)
		{
			return LCError("shape info type inconsistent");
		}
		OCTetCADShapeInfoMap* mapResult = new OCTetCADShapeInfoMap();
		mapResult->buildMap(sourceTet, targetTet);
		*result = mapResult;
		return err;
	}
*/
	return LCError("Shape info type not specified");

}
