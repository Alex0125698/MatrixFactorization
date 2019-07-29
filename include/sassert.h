#pragma once
#include <cassert>

void sassert(bool cond)
{
	if (!cond)
		assert(false);
}