#include "SimplexMethod.h"
#include <iostream>


int main()
{
	setlocale(LC_ALL, "russian");

	std::vector<std::vector<double>> v{
		std::vector<double>{5.0,  1.0, 9.0, 12.0, 1500.0},
		std::vector<double>{2.0,  3.0, 4.0,  1.0, 1000.0},
		std::vector<double>{3.0,  2.0,		5.0, 10.0,  300.0},
		std::vector<double>{1.0, -1.0/6.0, 0.0,  0.0,    0.0}
	};

	std::vector<std::vector<double>> v2{
		std::vector<double>{1, 2, 6},
		std::vector<double>{2, 1, 8},
		std::vector<double>{-1, 1, 1},
		std::vector<double>{0, 1, 2},
	};

	SimplexMethod sm(v, std::vector<double>{12, 5, 15, 10},std::vector<double>{ 0, 0, 0, 0});

	try
	{
		sm.execute();
		auto b_v = sm.getBasicVector();

		for (const auto& it : b_v)
			std::cout << it << " ";
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}

	SimplexMethod sm2(v2, std::vector<double>{3, 2}, std::vector<double>{ 0, 0, 0, 0});

	try
	{
		sm2.execute();
		auto b_v = sm2.getBasicVector();

		for (const auto& it : b_v)
			std::cout << it << " ";

	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}

	return 0;
}
