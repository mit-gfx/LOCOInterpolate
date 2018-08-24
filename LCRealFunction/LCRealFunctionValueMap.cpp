#include <iostream>
#include "LCFunctionValueMap.h"
#include "LCRealFunctionValue.h"
#include "LCError.h"
#include <ctime>
#include "LCRealFunctionValueMap.h"

LCRealFunctionValueMap::LCRealFunctionValueMap()
{

}

LCError LCRealFunctionValueMap::newShapeFromMap(LCFunctionValue *sourceShape, LCFunctionValue **result)
{
	LCError err;
	LCRealFunctionValue* testSource = dynamic_cast<LCRealFunctionValue*>(sourceShape);
	if (testSource)
	{
		*result = new LCRealFunctionValue(testSource->val);
		return err;
	}

	return LCError("wrong shape info type");

}