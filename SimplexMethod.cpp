#include "SimplexMethod.h"


double SimplexMethod::gaussianElimination(std::pair<int, int> index_target_element, std::pair<int, int> index_resolving_element)
{
	// row col
	return  (matrix.at(index_target_element.first).at(index_target_element.second) * matrix.at(index_resolving_element.first).at(index_resolving_element.second) -
		matrix.at(index_target_element.first).at(index_resolving_element.second) * matrix.at(index_resolving_element.first).at(index_target_element.second)) /
		matrix.at(index_resolving_element.first).at(index_resolving_element.second);
}


std::vector<double> SimplexMethod::calcRelativeEstimates()
{
	std::vector<double> relative_estimates; // (matrix.at(0).size());

	std::vector<double> tmp;

	for (int col = 0; col < matrix.at(0).size() - 1; col++)
	{
		for (int row = 0; row < matrix.size(); row++)
			tmp.push_back(matrix.at(row).at(col));

		relative_estimates.push_back(scalarProduct(basic_variables, tmp) - target_variables.at(col));
		tmp.clear();
	}

	return relative_estimates;
}


double SimplexMethod::calcQ()
{
	std::vector<double> tmp;
	for (int row = 0; row < matrix.size(); row++)
		tmp.push_back(matrix.at(row).at(matrix.at(0).size() - 1));

	return scalarProduct(basic_variables, tmp);
}


void SimplexMethod::complementMatrix()
{
	auto relative_estimates = calcRelativeEstimates();
	double q = calcQ();

	relative_estimates.push_back(q);
	matrix.push_back(relative_estimates);

	//memory clear
	for (auto& it : matrix)
		it.shrink_to_fit();

	matrix.shrink_to_fit();
}


void SimplexMethod::printMatrix()
{
	for (const auto& col : matrix)
	{
		for (const auto& row : col)
			std::cout << std::setw(10) << row;
		std::cout << std::endl;
	}
}


void SimplexMethod::execute()
{
	complementMatrix();

	std::cout << "Исходная матрица" << std::endl;
	printMatrix();

	std::cout << std::endl << "Базисные переменные: " << std::endl;
	for (std::size_t i = 0; i < basic_variables.size(); i++)
		std::cout << "x_" << names_basic_variables[i] + 1 << " = " << basic_variables[i] << std::endl;

	std::cout << std::endl;

	std::cout << "Небазисные переменные: " << std::endl;
	for (std::size_t i = 0; i < target_variables.size(); i++)
		std::cout << "x_" << names_target_variables[i] + 1 << " = " << target_variables[i] << std::endl;

	int i = 0;
	while (formationNewSimplexMatrix())
	{
		std::cout << std::endl;
		std::cout << "Итерация " << i << std::endl;
		printMatrix();

		std::cout << std::endl << "Базисные переменные: " << std::endl;
		for (std::size_t i = 0; i < basic_variables.size(); i++)
			std::cout << "x_" << names_basic_variables[i] + 1 << " = " << basic_variables[i] << std::endl;

		std::cout << std::endl;

		std::cout << "Небазисные переменные: " << std::endl;
		for (std::size_t i = 0; i < target_variables.size(); i++)
			std::cout << "x_" << names_target_variables[i] + 1 << " = " << target_variables[i] << std::endl;

		i++;
	}

	std::cout << std::endl;
	std::cout << "Итерация " << i << std::endl;
	printMatrix();

	std::cout << std::endl << "Базисные переменные: " << std::endl;
	for (std::size_t i = 0; i < basic_variables.size(); i++)
		std::cout << "x_" << names_basic_variables[i] + 1 << " = " << basic_variables[i] << std::endl;

	std::cout << std::endl;

	std::cout << "Небазисные переменные: " << std::endl;
	for (std::size_t i = 0; i < target_variables.size(); i++)
		std::cout << "x_" << names_target_variables[i] + 1 << " = " << target_variables[i] << std::endl;

	std::cout << std::endl << "f-max = " << matrix.at(matrix.size() - 1).at(matrix.at(matrix.size() - 1).size() - 1) << " (через симплекс-матрицу)" << std::endl;
	
	std::vector<double> a_0{};
	for (std::size_t row = 0; row < matrix.size() - 1; row++)
		a_0.push_back(matrix[row][3]);

	std::cout << "f-max = " << scalarProduct(basic_variables, a_0) << " (через скалярное произведение A0 на Cb)" << std::endl;

}


int SimplexMethod::findResolutionCol()
{
	auto min_iter = std::min_element(matrix.at(matrix.size() - 1).begin(), matrix.at(matrix.size() - 1).end());
	int index = std::distance(matrix.at(matrix.size() - 1).begin(), min_iter);

	return index;
}


int SimplexMethod::findsolutionRow(int index_resolution_col)
{
	std::pair<int, double> min_index{ -1, -1 };

	for (int i = 0; i < matrix.size(); i++)
		if (matrix.at(i).at(index_resolution_col) > 0 && matrix.at(i).at(matrix.at(i).size() - 1) / matrix.at(i).at(index_resolution_col) != 0)
		{
			min_index.first = i;
			min_index.second = matrix.at(i).at(matrix.at(i).size() - 1) / matrix.at(i).at(index_resolution_col);
			break;
		}

	if (min_index.first == -1)
		return -1;

	for (int i = 0; i < matrix.size(); i++)
		if (matrix.at(i).at(index_resolution_col) > 0 && matrix.at(i).at(matrix.at(i).size() - 1) / matrix.at(i).at(index_resolution_col) < min_index.second && matrix.at(i).at(matrix.at(i).size() - 1) / matrix.at(i).at(index_resolution_col) != 0)
		{
			min_index.first = i;
			min_index.second = matrix.at(i).at(matrix.at(i).size() - 1) / matrix.at(i).at(index_resolution_col);
		}

	return min_index.first;
}


bool SimplexMethod::formationNewSimplexMatrix()
{
	int resolution_col = findResolutionCol();				//Номер разрешающего столбца
	int resolution_row = findsolutionRow(resolution_col);	//Номер разрешающей строки

	std::cout << std::endl;
	std::cout << "Разрешающий столбец: " << resolution_col << std::endl;
	std::cout << "Разрешающая строка: " << resolution_row << std::endl;


	if (resolution_row < 0 || resolution_row == -1)
		throw std::runtime_error("The problem is not solved.");

	//swap(target_variables.at(resolution_col), basic_variables.at(resolution_row));

	//row col
	std::pair<int, int> a_rs_coords(resolution_row, resolution_col); //Индекс разрешающего элемента;

	//Расчет элементов по правилу прямоугольника (включая относительные оценки и значение целевой функции Q) - в самом начале, чтобы не хранить предыдущее состоянии элементов
	for (int row = 0; row < matrix.size(); row++)
		for (int col = 0; col < matrix.at(row).size(); col++)
			if (col != resolution_col && row != resolution_row)
				matrix.at(row).at(col) = gaussianElimination(std::make_pair(row, col), a_rs_coords);

	//Замена переменных вектора
	swapValues(resolution_row, resolution_col);
	swapNames(resolution_row, resolution_col);

	//Элементы разрешающей строки делятся на разрешающий элемент
	for (int i = 0; i < matrix.at(resolution_row).size(); i++)
		if (i != resolution_col)
			matrix.at(resolution_row).at(i) = matrix.at(resolution_row).at(i) / matrix.at(a_rs_coords.first).at(a_rs_coords.second);

	//Элементы разрешающего столбца делятся на разрешающий элемент и меняют знак
	for (int i = 0; i < matrix.size(); i++)
		if (i != resolution_row)
			matrix.at(i).at(resolution_col) = (-1) * (matrix.at(i).at(resolution_col) / matrix.at(a_rs_coords.first).at(a_rs_coords.second));

	//Разрешающий элемент заменяется на обратный
	matrix.at(a_rs_coords.first).at(a_rs_coords.second) = 1 / matrix.at(a_rs_coords.first).at(a_rs_coords.second);

	//Условие выхода, true - есть отрицательные элементы
	for (const auto& it : matrix.at(matrix.size() - 1))
		if (it < 0)
			return true;

	return false;
}


void SimplexMethod::swapValues(int first, int second)
{
	double tmp = basic_variables.at(first);
	basic_variables.at(first) = target_variables.at(second);
	target_variables.at(second) = tmp;
}

void SimplexMethod::swapNames(int first, int second)
{
	int tmp = names_basic_variables.at(first);
	names_basic_variables.at(first) = names_target_variables.at(second);
	names_target_variables.at(second) = tmp;
}


template<typename T>
T SimplexMethod::scalarProduct(const std::vector<T>& v1, const std::vector<T>& v2)
{
	if (v1.size() != v2.size())
		throw std::runtime_error("The dimension of the vectors is different.");

	T result = 0;

	for (int i = 0; i < v1.size(); i++)
		result += v1.at(i) * v2.at(i);

	return result;
}


std::vector<double> SimplexMethod::getBasicVector()
{
	return this->basic_variables;
}


std::vector<double> SimplexMethod::getTargetVector()
{
	return this->target_variables;
}


std::vector<int> SimplexMethod::getIndexBasicVector()
{
	return this->names_basic_variables;
}


std::vector<int> SimplexMethod::getIndexTargetVector()
{
	return this->names_target_variables;
}


std::vector<std::vector<double>> SimplexMethod::getBasicMatrix()
{
	return this->basic;
}


std::vector<std::vector<double>> SimplexMethod::getMatrix()
{
	return this->matrix;
}