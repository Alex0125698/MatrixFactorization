#pragma once
#include <cstdint>
#include <cassert>
#include <vector>
#include <ostream>
#include <algorithm>
#include <random>
#include <initializer_list>
// not used
#include <bitset>
#include <valarray>

void sassert(bool cond)
{
	if (!cond)
		assert(false);
}

// 2nd/Column index sequential in memory
class Matrix
{
public:
	Matrix() = default;

	Matrix(const std::initializer_list<std::initializer_list<double>> list)
	{
		if (list.size() == 0) return;
		sassert(list.begin()->size() != 0);

		rows = list.size();
		columns = list.begin()->size();
		data.reserve(rows * columns);

		for (auto& i : list)
		{
			sassert(i.size() == columnCount());
			for (auto& j : i)
			{
				data.push_back(j);
			}
		}
	}

	Matrix(size_t rows, size_t columns)
		: rows{ rows }, columns{ columns }, data(rows* columns)
	{
		sassert((rows != 0 && columns != 0) || (rows == 0 && columns == 0));
	}

	bool empty() const
	{
		return data.empty();
	}

	struct msize
	{
		size_t rows;
		size_t columns;
	};

	msize size() const
	{
		return msize{ rows,columns };
	}

	size_t elementCount() const
	{
		return data.size();
	}

	size_t rowCount() const
	{
		return rows;
	}

	size_t columnCount() const
	{
		return columns;
	}

	auto& operator()(size_t row, size_t column)
	{
		sassert(!empty());
		sassert(row < rowCount());
		sassert(column < columnCount());
		return data[row * columns + column];
	}

	const auto& operator()(size_t row, size_t column) const
	{
		sassert(!empty());
		sassert(row < rowCount());
		sassert(column < columnCount());
		return data[row * columns + column];
	}

	struct range
	{
		range(const size_t first_index, const size_t last_index)
			:first_index{ first_index }, last_index{ last_index }
		{
			sassert(last_index >= first_index);
		}
		const size_t first_index{ 0 };
		const size_t last_index{ 0 };
		auto count() const
		{
			return last_index - first_index + 1;
		}
	};

	auto slice(const range rows, const range columns)
	{
		sassert(rows.last_index < rowCount());
		sassert(columns.last_index < columnCount());

		Matrix res(rows.count(), columns.count());

		for (size_t r = 0; r < res.rowCount(); ++r)
		{
			const auto in_begin = data.begin() + (r + rows.first_index) * columnCount() + columns.first_index;
			const auto in_end = data.begin() + (r + rows.first_index) * columnCount() + columns.last_index + 1;
			const auto out_begin = res.data.begin() + r * res.columnCount();

			std::copy(in_begin, in_end, out_begin);
		}

		return std::move(res);
	}

	auto operator()(const range rows, const range columns)
	{
		return slice(rows, columns);
	}

	std::vector<double> data;
	size_t rows{ 0 };
	size_t columns{ 0 };

	static auto identity(size_t size)
	{
		Matrix res(size, size);
		for (size_t i = 0; i < res.data.size(); i += res.rowCount() + 1)
			res.data[i] = 1;
		return res;
	}
};

auto transpose(const Matrix m)
{
	Matrix res(m.columnCount(), m.rowCount());
	size_t res_index = 0;
	size_t start_index = 0;

	for (size_t i = 0; i < m.data.size(); ++i)
	{
		res.data[res_index] = m.data[i];
		res_index += m.rowCount();

		if (res_index >= m.data.size())
		{
			++start_index;
			res_index = start_index;
		}
	}

	return res;
}

auto& operator<<(std::ostream & out, const Matrix & m)
{
	for (size_t r = 0; r < m.rowCount(); ++r)
	{
		out << '\n';
		for (size_t c = 0; c < m.columnCount(); ++c)
		{
			out << ' ' << m(r, c);
		}
	}
	out << '\n';
	return out;
}

bool operator==(const Matrix::msize s1, const Matrix::msize s2)
{
	return (s1.rows == s2.rows) && (s1.columns == s2.columns);
}

// Elementwise addition of two matrices
auto operator+(const Matrix m1, const Matrix m2)
{
	sassert(m1.size() == m2.size());
	sassert(m1.rowCount() > 0);
	Matrix ret = m1;
	for (size_t i = 0; i < ret.data.size(); ++i)
		ret.data[i] += m2.data[i];
	return std::move(ret);
}

// Elementwise subtraction of two matrices
auto operator-(const Matrix m1, const Matrix m2)
{
	sassert(m1.size() == m2.size());
	sassert(m1.rowCount() > 0);
	Matrix ret = m1;
	for (size_t i = 0; i < ret.data.size(); ++i)
		ret.data[i] -= m2.data[i];
	return std::move(ret);
}

// Elementwise multiplication of two matrices
auto operator*(const Matrix m1, const Matrix m2)
{
	sassert(m1.size() == m2.size());
	sassert(m1.rowCount() > 0);
	Matrix ret = m1;
	for (size_t i = 0; i < ret.data.size(); ++i)
		ret.data[i] *= m2.data[i];
	return std::move(ret);
}

// Elementwise division of two matrices
auto operator/(const Matrix m1, const Matrix m2)
{
	sassert(m1.size() == m2.size());
	sassert(m1.rowCount() > 0);
	Matrix ret = m1;
	for (size_t i = 0; i < ret.data.size(); ++i)
		ret.data[i] /= m2.data[i];
	return std::move(ret);
}

// Matrix-wise multiplication
// Complexity: O(elementCount(m1)*columnCount(m2)); ~10 Sec: ?; no multithread, no SSE, no cache localisation
auto matmul(const Matrix m1, const Matrix m2)
{
	sassert(m1.columnCount() == m2.rowCount());
	Matrix res(m1.rowCount(), m2.columnCount());

	for (size_t r = 0; r < res.rowCount(); ++r)
	{
		for (size_t c = 0; c < res.columnCount(); ++c)
		{
			for (size_t k = 0; k < m1.columnCount(); ++k)
			{
				res(r, c) += m1(r, k) * m2(k, c);
			}
		}
	}

	return std::move(res);
}

// Matrix-wise multiplication
// Complexity: O(elementCount(m1)*columnCount(m2)); ~10 Sec: ?; no multithread, no SSE, no cache localisation
auto matmul2(const Matrix m1, const Matrix m2)
{
	sassert(m1.columnCount() == m2.rowCount());
	Matrix res(m1.rowCount(), m2.columnCount());

	auto m2_trans = transpose(m2);

	for (size_t r = 0; r < res.rowCount(); ++r)
	{
		for (size_t c = 0; c < res.columnCount(); ++c)
		{
			for (size_t k = 0; k < m1.columnCount(); ++k)
			{
				res(r, c) += m1(r, k) * m2_trans(c, k);
			}
		}
	}

	return std::move(res);
}

// Raise a matrix to some integer power (this uses matrix-wise multiplication)
// Complexity: O(log(pow)*elementCount(m)^(3/2)); ~10 Sec: m(400,400) ^ 2000000000; no multithread, no SSE
auto pow2(const Matrix m, const int32_t pow)
{
	sassert(m.rowCount() == m.columnCount());
	sassert(pow > 0); // negative powers not implemented yet
	if (pow == 0) return Matrix::identity(m.rowCount());

	Matrix pow_of_2 = m;
	Matrix working_mat;

	if (pow & 1u)
		working_mat = m;
	else
		working_mat = Matrix::identity(m.rowCount());

	uint32_t i = 2;
	while (i <= pow)
	{
		pow_of_2 = matmul2(pow_of_2, pow_of_2);
		if (i & pow)
			working_mat = matmul2(pow_of_2, working_mat);
		i <<= 1u;
	}

	return std::move(working_mat);
}

auto pow3(const Matrix m, const int32_t pow)
{
	sassert(m.rowCount() == m.columnCount());
	sassert(pow > 0); // negative powers not implemented yet
	if (pow == 0) return Matrix::identity(m.rowCount());

	Matrix result = m;
	for (uint32_t i = 1; i < pow; ++i)
		result = matmul(result, m);

	return std::move(result);
}

auto operator^(const Matrix m, const int32_t pow)
{
	sassert(m.rowCount() == m.columnCount());
	sassert(pow > 0); // negative powers not implemented yet
	if (pow == 0) return Matrix::identity(m.rowCount());

	Matrix result = m;
	int32_t pow_tmp = pow;
	int32_t tmp = 1;

	while (pow_tmp - tmp > 0)
	{
		result = matmul(result, result);
		pow_tmp -= tmp;
		tmp *= 2;
	}

	while (pow_tmp > 0)
	{
		result = matmul(result, m);
		--pow_tmp;
	}

	return std::move(result);
}

constexpr uint64_t sign_mask     = 0b10000000'00000000'00000000'00000000'00000000'00000000'00000000'00000000;
constexpr uint64_t exponent_mask = 0b01111111'11110000'00000000'00000000'00000000'00000000'00000000'00000000;
constexpr uint64_t mantissa_mask = 0b00000000'00001111'11111111'11111111'11111111'11111111'11111111'11111111;

void randomize(Matrix& m, const double min, const double max)
{
	// high-quality random-number generator
	static std::mt19937 mt(1);

	for (auto& d : m.data)
	{
		uint64_t tmp = uint64_t(mt()) << 32ull | uint64_t(mt());
		union di
		{
			double d;
			uint64_t i;
		};

		di num;
		num.i = tmp & (mantissa_mask);
		num.i |= 0b00111111'11110000'00000000'00000000'00000000'00000000'00000000'00000000;
		//d = (num.d - 1.0) * (max - min) - min;
		d = 1.013523 * (num.d - 0.5) / m.rowCount();
	}
}