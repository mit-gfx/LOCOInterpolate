#include <iostream>
#include "LCFunctionValue.h"
#include "LCError.h"
#include "time.h"
#include "LCRealFunctionValue.h"

/*bool LCRealFunctionValue::isTetCAD()
{
	return false;
}*/


LCRealFunctionValue::LCRealFunctionValue(double val) : LCFunctionValue()
{
	this->val = val;
}


LCError LCRealFunctionValue::computeDifference(LCFunctionValue *other, double *result)
{
	LCRealFunctionValue* otherTest = dynamic_cast<LCRealFunctionValue*>(other);
	if (!otherTest)
	{
		return LCError("error in dynamic casting");
	}
	*result = std::abs(val - otherTest->val);

	return LCError();
}


LCError LCRealFunctionValue::computeErrorVector(LCFunctionValue *other, Eigen::VectorXd *result)
{
	LCRealFunctionValue* otherTest = dynamic_cast<LCRealFunctionValue*>(other);
	if (!otherTest)
	{
		return LCError("error in dynamic casting");
	}
	*result = std::abs(val - otherTest->val) * Eigen::VectorXd::Ones(1);

	return LCError();
}



LCError LCRealFunctionValue::add(double weight, LCFunctionValue *other)
{
	LCRealFunctionValue* otherTest = dynamic_cast<LCRealFunctionValue*>(other);
	if (!otherTest)
	{
		return LCError("error in dynamic casting");
	}
	val +=  weight*otherTest->val;

	return LCError();
}

LCError LCRealFunctionValue::getHomeophicShapeInfo(LCFunctionValue *source, LCFunctionValue **result)
{
	LCRealFunctionValue* sourceTest = dynamic_cast<LCRealFunctionValue*>(source);
	if (!sourceTest)
	{
		return LCError("error in dynamic casting");
	}
	*result = new LCRealFunctionValue(this->val);

	return LCError();
}

LCError LCRealFunctionValue::addBoundaryConditions(bool useForceAngle)
{
	return LCError();
}


void LCRealFunctionValue::addNewPrecomputedPhysics(std::string name, std::vector<double> &values)
{
}


void LCRealFunctionValue::log()
{
	std::cout << "val = " << val << std::endl;
}

