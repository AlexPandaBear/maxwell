#pragma once

#include <map>

class SparseMatrix
{
private:
	size_t m_size;
	std::map<std::pair<size_t, size_t>, double> m_coefs;

public:
	SparseMatrix(size_t size);
	~SparseMatrix();

	size_t get_size() const;
	void set_coef(size_t i, size_t j, double coef);
	void add_to_coef(size_t i, size_t j, double delta_coef);
	double get_coef(size_t i, size_t j) const;
};