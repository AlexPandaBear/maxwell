#pragma once

#include <string>
#include "DataKeeper.hxx"

class DataProcessor
{
private:
	DataKeeper& m_data;

	bool m_energy_ready;
	bool m_poynting_ready;

	std::vector<ScalarField> m_energy;
	std::vector<VectorField> m_poynting;
	std::vector<ScalarField> m_poynting_norm;

	void compute_energy();
	void compute_poynting();
	void compute_all();

public:
	DataProcessor(DataKeeper& data);
	~DataProcessor();

	ScalarField const& get_energy_density(size_t step);
	VectorField const& get_poynting_vector(size_t step);
	ScalarField const& get_poynting_vector_norm(size_t step);
};