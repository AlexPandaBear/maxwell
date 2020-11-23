#include "Matrix.hxx"

Matrix::Matrix(size_t nb_rows, size_t nb_columns) :
	m_nb_rows(nb_rows),
	m_nb_columns(nb_columns),
	ptr_coefs(new double[nb_rows*nb_columns]) {}

Matrix::Matrix(size_t size) : Matrix(size, size) {}

Matrix::~Matrix() {}

double Matrix::get_coef(size_t line_index, size_t column_index) const
{
	return ptr_coefs[m_nb_columns*line_index + column_index];
}

void Matrix::set_coef(size_t line_index, size_t column_index, double value)
{
	ptr_coefs[m_nb_columns*line_index + column_index] = value;
}

void Matrix::add_to_coef(size_t line_index, size_t column_index, double value)
{
	ptr_coefs[m_nb_columns*line_index + column_index] += value;
}

void Matrix::multiply_coef(size_t line_index, size_t column_index, double value)
{
	ptr_coefs[m_nb_columns*line_index + column_index] *= value;
}

void Matrix::reset_as_diagonal(double value)
{
	if (m_nb_rows != m_nb_columns)
	{
		throw std::logic_error("Square matrix required");
	}

	reset_as_uniform(0.);

	for (size_t i = 0; i < m_nb_rows; i++)
	{
		set_coef(i, i, value);
	}
}

void Matrix::reset_as_uniform(double value)
{
	for (size_t i = 0; i < m_nb_rows; i++)
	{
		for (size_t j = 0; j < m_nb_columns; j++)
		{
			set_coef(i, j, value);
		}
	}
}

size_t Matrix::get_nb_rows() const
{
	return m_nb_rows;
}

size_t Matrix::get_nb_columns() const
{
	return m_nb_columns;
}

std::pair<size_t, size_t> Matrix::get_shape() const
{
	return {m_nb_rows, m_nb_columns};
}

void Matrix::perform_matrix_vector_product(std::vector<double> const& V, std::vector<double>& R) const
{
	if (m_nb_rows != R.size() or m_nb_columns != V.size())
	{
		throw std::invalid_argument("Vectors sizes don't agree with matrix shape");
	}


	tbb::parallel_for((size_t) 0, m_nb_rows, [&](size_t k)
	//for (size_t k = 0; k < m_nb_rows; k++)
	{
		R[k] = 0.;

		for (size_t i = 0; i < m_nb_columns; i++)
		{
			R[k] += get_coef(k,i)*V[i];
		}
	}
	);
}

double Matrix::perform_gauss_seidel_iteration(std::vector<double>& X, std::vector<double> const& B) const
{
	if (m_nb_rows != m_nb_columns)
	{
		throw std::logic_error("Square matrix required");
	}

	if (X.size() != m_nb_columns or B.size() != m_nb_rows)
	{
		throw std::invalid_argument("Vector size does not agree with matrix shape");
	}

	double tmp(0.);
	double residual_l2_norm2(0.);

	for (size_t k = 0; k < m_nb_rows; k++)
	{
		tmp = B[k];

		for (size_t i  = 0; i < k; i++)
		{
			tmp -= get_coef(k,i)*X[i];
		}

		for (size_t i  = k+1; i < m_nb_columns; i++)
		{
			tmp -= get_coef(k,i)*X[i];
		}

		tmp /= get_coef(k,k);
		residual_l2_norm2 += pow(tmp - X[k], 2);
		X[k] = tmp;
	}

	return residual_l2_norm2;
}