#ifndef _OC_FUNC_SHAPE_INFO_MAP_
#define _OC_FUNC_SHAPE_INFO_MAP_
#include <Eigen/Dense>
#include <vector> 

class OCParametricShape;
class LCError;
class LCFunctionValue;
class LCRealFunctionValue;

/* maps each point in the sourceShape (on a neighbor) to the targetShape (on the current cell)*/
class LCFunctionValueMap{
public:
	LCFunctionValueMap();
	virtual ~LCFunctionValueMap() {}
	static LCError newMapFromShapePair(LCFunctionValue *sourceShape, LCFunctionValue *targetShape, LCFunctionValueMap **result);
	virtual LCError newShapeFromMap(LCFunctionValue *sourceShape, LCFunctionValue **result) = 0;
};

#endif 
