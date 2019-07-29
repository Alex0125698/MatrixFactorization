#include <iostream>
#include "Timer.h"
#include "Matrix.h"

// Factorize V into H, W s.t. W x H = V
// Using Lee and Seung's multiplicative update rule
void factorise(const Matrix& inV, Matrix& outH, Matrix& outW, const uint32_t innerDim)
{
	Matrix W(inV.rowCount(),innerDim);
	Matrix H(innerDim, inV.columnCount());
	randomize(W, 0, 1);
	randomize(H, 0, 1);

	for (uint32_t n = 1; n <= 100000; ++n)
	{
		H = H * matmul(transpose(W),inV) / matmul(transpose(W), matmul(W, H));
		W = W * matmul(inV, transpose(H)) / matmul(matmul(W,H),transpose(H));
	}

	outH = std::move(H);
	outW = std::move(W);
}

void tests()
{
	Matrix m3(60, 60);
	Matrix m4(60, 60);
	//std::cout << m3 << std::endl;
	randomize(m3, 1, 2);

	Timer timer;
	//std::cout << matmul(m3, m3)(0,0) << '\n';
	//std::cout << timer << "(matmul original)" << std::endl;
	//timer.restart();
	//std::cout << matmul2(m4, m4)(0,0) << '\n';
	//std::cout << timer << "(matmul cache)" << std::endl;
	//timer.restart();

/*	std::cout << (m3 ^ 20000)(0, 0) << std::endl;
	std::cout << "duration = " << timer << std::endl*/;
	//timer.restart();
	std::cout << (pow2(m3, 2000000000))(0, 0) << '\n';
	std::cout << "duration2 = " << timer << std::endl;
	timer.restart();
	//std::cout << (pow3(m3, 20000))(0, 0) << std::endl;
	//std::cout << "duration3 = " << timer << std::endl;

	Matrix m1 = { {0.1,0,0.7953,0.01} , {0,0.2,0,0.01}, {0.05,0.1,0,0.8} }; // 3x4
	Matrix m2 = { {1,2,3} , {4,5,6}, {7,8,9}, {0,10,11} }; // 4x3

	std::cout << "m2 slice = " << m2({ 1,2 }, { 1,2 }) << std::endl;
	std::cout << "m1 = " << m1 << std::endl;
	std::cout << "m2 = " << m2 << std::endl;
	std::cout << "m1 x m2 = " << (matmul(m1, m2) ^ 2) << std::endl;
	std::cout << "m1 x m2 = " << transpose(matmul(m1, m2) ^ 2) << std::endl;
	std::cout << Matrix::identity(4) << std::endl;
}

int main()
{
	std::ios_base::sync_with_stdio(false); // massive speedup

	Matrix W(9, 2);
	Matrix H(2,11);
	randomize(W, 0, 1);
	randomize(H, 0, 1);
	Matrix V = matmul(W, H);

	std::cout << "actual W = " << W << std::endl;
	std::cout << "actual H = " << H << std::endl;
	std::cout << "V = W x H = " << V << std::endl;

	Matrix estW, estH;
	Timer timer;
	factorise(V, estH, estW, 2);
	std::cout << timer << " (factorize run time)" << std::endl;

	std::cout << "W estimate = " << estW << std::endl;
	std::cout << "H estimate = " << estH << std::endl;
	std::cout << "V estimate = " << matmul(estW,estH) << std::endl;

}