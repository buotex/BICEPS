#include "Timer.h"
namespace Pepsplice{
Timer::Timer()
{
	std::time(&starttime);
	lasttime = starttime;
}

Timer::~Timer()
{
}

double Timer::timeSinceStart()
{
	std::time_t thistime;
	std::time(&thistime);
	return thistime - starttime;;
}

double Timer::startTime()
{
	return starttime;
}

double Timer::timeSinceLast()
{
	std::time_t thistime;
	std::time(&thistime);
	std::time_t timesincelast = thistime - lasttime;
	lasttime = thistime;
	return timesincelast;
}

double Timer::timeSinceLastWithoutReset()
{
	std::time_t thistime;
	std::time(&thistime);
	std::time_t timesincelast = thistime - lasttime;
	return timesincelast;
}

void Timer::reset()
{
	time(&starttime);
	lasttime = starttime;
}
}