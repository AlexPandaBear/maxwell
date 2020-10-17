#include "SparseMatrix.hxx"

SparseMatrix::SparseMatrix(size_t size) :
	m_size(size) {}

SparseMatrix::~SparseMatrix() {}

size_t SparseMatrix::get_size() const
{
	return m_size;
}

void SparseMatrix::set_boundary_equation(size_t line_nb)
{
	for (std::map<std::pair<size_t, size_t>, double>::iterator it = m_coefs.begin(); it != m_coefs.end(); it++)
	{
		if (it->first.first == line_nb)
		{
			m_coefs.erase(it);
		}
	}

	set_coef(line_nb, line_nb, 1.);
}

void SparseMatrix::set_coef(size_t i, size_t j, double coef)
{
	try
	{
		m_coefs.at(std::pair<size_t, size_t>(i, j)) = coef;
	}
	catch (std::out_of_range)
	{
		m_coefs.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t>(i, j), coef));
	}
}

void SparseMatrix::add_to_coef(size_t i, size_t j, double delta_coef)
{
	try
	{
		m_coefs.at(std::pair<size_t, size_t>(i, j)) += delta_coef;
	}
	catch (std::out_of_range)
	{
		m_coefs.insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t>(i, j), delta_coef));
	}
}

double SparseMatrix::get_coef(size_t i, size_t j) const
{
	try
	{
		return m_coefs.at(std::pair<size_t, size_t>(i, j));
	}
	catch (std::out_of_range)
	{
		return 0.;
	}
}

std::vector<std::pair<size_t, double>> SparseMatrix::get_reduced_line_wo_diag(size_t i) const
{
	std::vector<std::pair<size_t, double>> line;
	double coef;

	for (size_t j = 0; j < m_size; j++)
	{
		coef = get_coef(i, j);
		if (coef != 0.)
		{
			line.push_back(std::pair<size_t, double>(j, coef));
		}
	}

	return line;
}