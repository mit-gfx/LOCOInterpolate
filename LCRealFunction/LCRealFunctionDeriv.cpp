#include <iostream>
#include "LCError.h"
#include "LCRealFunctionDeriv.h"
#include "LCRealFunctionValue.h"

LCRealFunctionDeriv::LCRealFunctionDeriv(double val) : LCFunctionDeriv()
{
	this->val = val;
}

LCError LCRealFunctionDeriv::add(double weight, LCFunctionValue *other)
{
	LCRealFunctionValue* otherTest = dynamic_cast<LCRealFunctionValue*>(other);
	if (!otherTest)
	{
		return LCError("error in dynamic casting");
	}
	val += weight*otherTest->val;

	return LCError();
}

void LCRealFunctionDeriv::log()
{
	std::cout << "val = " << val << std::endl;
}
