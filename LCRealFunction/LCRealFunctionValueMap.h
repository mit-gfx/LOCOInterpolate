#ifndef _OC_REAL_VALUE_MAP_
#define _OC_REAL_VALUE_MAP_
#include <Eigen/Dense>
#include <vector> 

class OCParametricShape;
class LCError;
class LCFunctionValue;

class LCRealFunctionValueMap : public LCFunctionValueMap
{
public:
	LCRealFunctionValueMap();
	LCError newShapeFromMap(LCFunctionValue *sourceShape, LCFunctionValue **result);
	LCError mapPrecomputedPhysics(LCFunctionValue *sourceShape, LCFunctionValue *target);
};

#endif 
