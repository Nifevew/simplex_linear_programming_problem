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
	/// ������ ������
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
	/// ��������-�������
	/// </summary>
	std::vector<std::vector<double>> matrix;

	/// <summary>
	/// ������� ������
	/// </summary>
	std::vector<double> target_variables;

	/// <summary>
	/// �������� ������
	/// </summary>
	std::vector<double> basic_variables;

	std::vector <int> names_basic_variables;
	std::vector <int> names_target_variables;

	std::vector<std::vector<double>> basic;


	/// <summary>
	///  ����� �������� �������� �� ������ �������������� (���������� ������).
	///  ����������� �� �������: 
	/// 		    a_ij * a_rs - a_is * a_rj
	///		a_ij = ---------------------------  ,
	///				        a_rs
	///		��� a_ij - ������� �������;
	///			a_rs - ����������� �������;
	///			a_is - �������, ������� �� ����������� ������������ ������� � ������� ������;
	///			a_rj - �������, ������� �� ����������� ����������� ������ � �������� �������.
	/// </summary>
	/// <param name="index_target_element">: ������������ �������� ��������</param>
	/// <param name="index_resolving_element">: ������������ ������������ ��������</param>
	/// <returns>����� ��������</returns>	
	double gaussianElimination(std::pair<int, int> index_target_element, std::pair<int, int> index_resolving_element);

	/// <summary>
	/// ����� ������������ �������
	/// </summary>
	/// <retuens>������ ������������ �������</returns>
	int findResolutionCol();

	/// <summary>
	/// ����� ����������� ������
	/// </summary>
	/// <param name="index_resolution_col">: ������ ������������ �������</param>
	/// <returns>������ ����������� ������</returns>
	int findsolutionRow(int index_resolution_col);

	/// <summary>
	/// ������ ����� ��������-�������
	/// </summary>
	/// <returns>True - ���� ���� �������������, false - ������ �������������</returns>
	bool formationNewSimplexMatrix();

	/// <summary>
	/// ������ ������������� ������
	/// </summary>
	/// <returns>������ � �������������� �������� �� ��������</returns>
	std::vector<double> calcRelativeEstimates();

	/// <summary>
	/// ������ �������� ������� �������
	/// </summary>
	/// <returns> �������� ������� �������</returns>
	double calcQ();

	/// <summary>
	/// ��������� �������� �������.
	/// ���������� ������������� ������ � �������� ������� �������
	/// </summary>
	void complementMatrix();

	/// <summary>
	/// �������� ������� ��������
	/// </summary>
	void swapValues(int first, int second);

	void swapNames(int first, int second);

	/// <summary>
	/// ��������� ������������ ��������
	/// </summary>
	/// <returns>����� - ��������� ������������ ��������</returns>
	template<typename T>
	T scalarProduct(const std::vector<T>& v1, const std::vector<T>& v2);

	/// <summary>
	/// ������� �������
	/// </summary>
	void printMatrix();

};

