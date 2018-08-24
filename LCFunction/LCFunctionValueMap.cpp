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

	return LCError("Shape info type not specified");

}
