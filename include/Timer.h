#pragma once
#include <chrono>
#include <ostream>

class Timer
{
public:
	Timer()
	{
		restart();
	}
	void restart()
	{
		m_time = std::chrono::high_resolution_clock::now();
	}
	double getDuration() const
	{
		auto tmp = std::chrono::high_resolution_clock::now();
		return std::chrono::duration<double>(tmp - m_time).count();
	}

private:
	std::chrono::high_resolution_clock::time_point m_time;
};

auto& operator<<(std::ostream& out, const Timer& timer)
{
	return out << timer.getDuration() << " sec";
}