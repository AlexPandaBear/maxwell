#include "SparseMatrix.hxx"

SparseMatrix::SparseMatrix(size_t size) :
	m_size(size) {}

SparseMatrix::SparseMatrix(SparseMatrix const& m) :
	m_size(m.get_size())
{
	std::vector<std::pair<size_t, double>> line;

	for (size_t l = 0; l < m_size; l++)
	{
		line = m.get_reduced_line_w_diag(l);

		for (std::pair<size_t, double> p : line)
		{
			set_coef(l, p.first, p.second);
		}
	}
}

SparseMatrix::~SparseMatrix() {}

size_t SparseMatrix::get_size() const
{
	return m_size;
}

void SparseMatrix::set_boundary_equation(size_t line_nb)
{
	for (size_t j = 0; j < line_nb; j++)
	{
		set_coef(line_nb, j, 0.);
	}

	set_coef(line_nb, line_nb, 1.);

	for (size_t j = line_nb+1; j < m_size; j++)
	{
		set_coef(line_nb, j, 0.);
	}
}

void SparseMatrix::set_coef(size_t i, size_t j, double coef)
{
	if (coef == 0.)
	{
		m_coefs.erase({i,j});
	}
	else
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

std::vector<std::pair<size_t, double>> SparseMatrix::get_reduced_line_w_diag(size_t i) const
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

std::vector<std::pair<size_t, double>> SparseMatrix::get_reduced_line_wo_diag(size_t i) const
{
	std::vector<std::pair<size_t, double>> line;
	double coef;

	for (size_t j = 0; j < i; j++)
	{
		coef = get_coef(i, j);
		if (coef != 0.)
		{
			line.push_back(std::pair<size_t, double>(j, coef));
		}
	}

	for (size_t j = i+1; j < m_size; j++)
	{
		coef = get_coef(i, j);
		if (coef != 0.)
		{
			line.push_back(std::pair<size_t, double>(j, coef));
		}
	}

	return line;
}

void SparseMatrix::display() const
{
	for (size_t i = 0; i < m_size; i++)
	{
		for (size_t j = 0; j < m_size; j++)
		{
			std::cout << get_coef(i, j) << "  ";
		}
		std::cout << std::endl;
	}
}

void SparseMatrix::operator+=(SparseMatrix const& m)
{
	if (m.get_size() != m_size)
	{
		throw std::invalid_argument("Matrix sizes must be equal to perform summation");
	}

	std::vector<std::pair<size_t, double>> line;

	for (size_t l = 0; l < m_size; l++)
	{
		line = m.get_reduced_line_w_diag(l);

		for (std::pair<size_t, double> p : line)
		{
			add_to_coef(l, p.first, p.second);
		}
	}
}

void SparseMatrix::operator-=(SparseMatrix const& m)
{
	if (m.get_size() != m_size)
	{
		throw std::invalid_argument("Matrix sizes must be equal to perform difference");
	}

	std::vector<std::pair<size_t, double>> line;

	for (size_t l = 0; l < m_size; l++)
	{
		line = m.get_reduced_line_w_diag(l);

		for (std::pair<size_t, double> p : line)
		{
			add_to_coef(l, p.first, -p.second);
		}
	}
}

void SparseMatrix::operator*=(double a)
{
	for (std::map<std::pair<size_t, size_t>, double>::iterator it = m_coefs.begin(); it != m_coefs.end(); it++)
	{
		it->second *= a;
	}
}

void SparseMatrix::operator/=(double a)
{
	if (a == 0)
	{
		throw std::invalid_argument("Matrix can't be divided by zero");
	}

	for (std::map<std::pair<size_t, size_t>, double>::iterator it = m_coefs.begin(); it != m_coefs.end(); it++)
	{
		it->second /= a;
	}
}

SparseMatrix SparseMatrix::operator+(SparseMatrix const& m) const
{
	SparseMatrix copy(*this);
	copy += m;
	return copy;
}

SparseMatrix SparseMatrix::operator-(SparseMatrix const& m) const
{
	SparseMatrix copy(*this);
	copy -= m;
	return copy;
}

SparseMatrix SparseMatrix::operator*(double a) const
{
	SparseMatrix copy(*this);
	copy *= a;
	return copy;
}

SparseMatrix SparseMatrix::operator/(double a) const
{
	SparseMatrix copy(*this);
	copy /= a;
	return copy;
}
