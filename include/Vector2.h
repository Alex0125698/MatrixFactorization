#pragma once
#include <initializer_list>
#include <algorithm>
#include "sassert.h"
/*
class Vector2
{
public:
	class iterator
	{
	public:
		auto operator*()
		{
			assert(ptr != nullptr);
			return *ptr;
		}
		auto operator++()
		{
			++ptr;
		}
		auto operator+(std::size_t off)
		{
			return iterator{ ptr + off };
		}
		double* __restrict ptr{ nullptr };
		friend class ::Vector2;
	};
public:
	Vector2() = default;
	Vector2(const std::initializer_list<double> list)
	{
		it_begin.ptr = new double[list.size()*2];
		it_end = it_begin + list.size();
		it_end_alloc = it_begin + list.size()*2;
		std::copy(list.begin(), list.end(), begin());
	}
	Vector2(const Vector2& other)
	{
		data = new double[other.size()];
		std::copy(other.begin(), other.end(), data);
	}
	~Vector2()
	{
		if (it_begin.ptr) delete it_begin.ptr;
	}
	auto begin()
	{
		return it_begin;
	}
	auto end()
	{
		return it_end;
	}


	iterator it_begin;
	iterator it_end;
	iterator it_end_alloc;



};
*/