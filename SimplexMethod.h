#pragma once

#include <string>
#include <vector>

#include <exception>
#include <locale>

#include <iostream>
#include <iomanip>

#include <algorithm>

class SimplexMethod
{
public:

	SimplexMethod(std::vector<std::vector<double>> input_matrix, std::vector<double> target_variables, std::vector<double> basic_variables)
	{
		this->matrix = input_matrix;

		this->target_variables = target_variables;
		this->basic_variables = basic_variables;

		int index_variables = 0;
		for (std::size_t i = 0; i < target_variables.size(); i++, index_variables++)
			names_target_variables.push_back(index_variables);

		for (std::size_t i = 0; i < basic_variables.size(); i++, index_variables++)
			names_basic_variables.push_back(index_variables);

		basic = std::vector<std::vector<double>>(basic_variables.size(), std::vector<double>(basic_variables.size()));

		for (std::size_t row = 0; row < basic_variables.size(); row++)
			basic[row][row] = 1;
	}


	/// <summary>
	/// Запуск метода
	/// </summary>
	void execute();


	std::vector<double> getBasicVector();

	std::vector<double> getTargetVector();

	std::vector<int> getIndexBasicVector();

	std::vector<int> getIndexTargetVector();

	std::vector<std::vector<double>> getBasicMatrix();

	std::vector<std::vector<double>> getMatrix();


private:
	/// <summary>
	/// Симплекс-матрица
	/// </summary>
	std::vector<std::vector<double>> matrix;

	/// <summary>
	/// Целевой вектор
	/// </summary>
	std::vector<double> target_variables;

	/// <summary>
	/// Базисный вектор
	/// </summary>
	std::vector<double> basic_variables;

	std::vector <int> names_basic_variables;
	std::vector <int> names_target_variables;

	std::vector<std::vector<double>> basic;


	/// <summary>
	///  Новое значение элемента по методу прямоугольника (Исключение Гаусса).
	///  Вычисляется по формуле: 
	/// 		    a_ij * a_rs - a_is * a_rj
	///		a_ij = ---------------------------  ,
	///				        a_rs
	///		где a_ij - целевой элемент;
	///			a_rs - разрешающий элемент;
	///			a_is - элемент, стоящий на пересечении разрешающего столбца и целевой строки;
	///			a_rj - элемент, стоящий на пересечении разрешающей строки и целевого столбца.
	/// </summary>
	/// <param name="index_target_element">: расположение целевого элемента</param>
	/// <param name="index_resolving_element">: расположение разрешающего элемента</param>
	/// <returns>Новое значение</returns>	
	double gaussianElimination(std::pair<int, int> index_target_element, std::pair<int, int> index_resolving_element);

	/// <summary>
	/// Поиск разрешающего столбца
	/// </summary>
	/// <retuens>Индекс разрешающего столбца</returns>
	int findResolutionCol();

	/// <summary>
	/// Поиск разрешающей строки
	/// </summary>
	/// <param name="index_resolution_col">: индекс разрешающего столбца</param>
	/// <returns>Индекс разрешающей строки</returns>
	int findsolutionRow(int index_resolution_col);

	/// <summary>
	/// Расчет новой симплекс-матрицы
	/// </summary>
	/// <returns>True - если есть отрицательные, false - только положительные</returns>
	bool formationNewSimplexMatrix();

	/// <summary>
	/// Расчет относительных оценок
	/// </summary>
	/// <returns>Вектор с относительными оценками по столбцам</returns>
	std::vector<double> calcRelativeEstimates();

	/// <summary>
	/// Расчет значения целевой функции
	/// </summary>
	/// <returns> Значение целевой функции</returns>
	double calcQ();

	/// <summary>
	/// Дособрать исходную матрицу.
	/// Рассчитать относительные оценки и значение целевой функции
	/// </summary>
	void complementMatrix();

	/// <summary>
	/// Поменять местами значения
	/// </summary>
	void swapValues(int first, int second);

	void swapNames(int first, int second);

	/// <summary>
	/// Скалярное произведение векторов
	/// </summary>
	/// <returns>Число - скалярное произведение векторов</returns>
	template<typename T>
	T scalarProduct(const std::vector<T>& v1, const std::vector<T>& v2);

	/// <summary>
	/// Вывести матрицу
	/// </summary>
	void printMatrix();

};

