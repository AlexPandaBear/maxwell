#pragma once

#include <cmath>
#include <vector>
#include <utility>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <tbb/parallel_for.h>

class Matrix
{
protected:
	size_t m_nb_rows;
	size_t m_nb_columns;

	std::unique_ptr<double[]> ptr_coefs;

public:
	Matrix(size_t nb_rows, size_t nb_columns);
	Matrix(size_t size);
	virtual ~Matrix();

	double get_coef(size_t line_index, size_t column_index) const;
	void set_coef(size_t line_index, size_t column_index, double value);
	void add_to_coef(size_t line_index, size_t column_index, double value);
	void multiply_coef(size_t line_index, size_t column_index, double value);

	void reset_as_diagonal(double value);
	void reset_as_uniform(double value);

	size_t get_nb_rows() const;
	size_t get_nb_columns() const;
	std::pair<size_t, size_t> get_shape() const;

	void perform_matrix_vector_product(std::vector<double> const& V, std::vector<double>& R) const;
	double perform_gauss_seidel_iteration(std::vector<double>& X, std::vector<double> const& B) const;

	void display() const;
};