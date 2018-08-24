#include "LCFunction.h"
#include "LCError.h"
#include <iostream>

LCFunction::LCFunction()
{
}

LCError LCFunction::getMidPoint(std::vector<double> *result)
{
	for (int i = 0; i < getNParams(); i++)
	{
		double minVal = getMinRange(i);
		double maxVal = getMaxRange(i);
		if (maxVal < minVal)
		{
			std::cout << "get MidPoint" << std::endl;
			return(LCError("ranges are invalid for get MidPoint"));
		}
		result->push_back(0.5 * (minVal + maxVal));
	}
	return LCError();
}

LCError LCFunction::validateParameters(const std::vector<double> &parameters) const
{
	if (parameters.size() != getNParams())
	{
		return LCError("number of parameters is incorrect");
	}

	for (size_t i = 0; i < parameters.size(); i++)
	{
		double minVal = getMinRange(i);
		double maxVal = getMaxRange(i);
		
		if ((parameters[i] < minVal) || (parameters[i] > maxVal))
		{
			std::cout << "minVal " << minVal << " maxVal " << maxVal << std::endl;
			std::cout << parameters[i] << std::endl;
			return(LCError("ranges are invalid for validateParameters"));
		}
	}
	return LCError();
}

LCError LCFunction::getRelativeParamVal(int i, double alpha, double *result )
{
	if ((i >= getNParams()) || (alpha < 0) || (alpha > 1))
	{
		return LCError("parameter out of range");
	}

	*result = getMinRange(i)*(1 - alpha) + getMaxRange(i)*alpha;
	return LCError();

}



std::vector<double> LCFunction::mapToStandardHypercube(const std::vector<double> & params)
{

	std::vector<double> result(params.size());
	for (int i = 0; i < params.size(); i++)
	{

		double alpha = (params[i] - getMinRange(i)) / (getMaxRange(i) - getMinRange(i));
		result[i] = alpha * 2 - 1;
	}

	return result;


}



std::vector<double> LCFunction::mapFromStandardHypercube(const std::vector<double> &unitParams)
{
	//std::cout << "mapping from standard hypercube " << std::endl;
	std::vector<double>result(unitParams.size());
	for (int i = 0; i < unitParams.size(); i++)
	{
		double alpha = (unitParams[i] + 1.) / 2.;
	//	std::cout << i << " min range " << getMinRange(i) << " " << getMaxRange(i) << std::endl;
		result[i] = getMinRange(i)*(1 - alpha) + getMaxRange(i)*alpha;
	}
	//std::cout << "mapping from standard hypercube done " << std::endl;

	return result;

}


void LCFunction::getRanges(std::vector<double> * ranges)
{
	for (int i = 0; i < getNParams(); i++)
	{
		ranges->push_back(getMinRange(i));
		ranges->push_back(getMaxRange(i));
	}
}
