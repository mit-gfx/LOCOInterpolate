#ifndef OC_TIMER
#define OC_TIMER

#include <vector>
#include <ctime>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
class LCError;

class LCTimer
{
public:
	LCTimer()
	{
		lastTime = clock();
	}
	void reset()
	{
		lastTime = clock();
	}
	void addMeasure(std::string tag)
	{
		times.push_back(clock() - lastTime);
		timeTags.push_back(tag);
		lastTime = clock();
	}
	void log()
	{
		for (int i = 0; i < times.size(); i++)
		{
			std::cout << "Time for " << timeTags[i] << ": " << times[i] / CLOCKS_PER_SEC << std::endl;
		}
	}
	void log(std::string tag)
	{
		times.push_back(clock() - lastTime);
		timeTags.push_back(tag);
		log();
	}
private:
	std::vector<int> times;
	std::vector<std::string> timeTags;
	double lastTime;
};

#endif
