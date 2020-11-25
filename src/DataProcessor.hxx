#pragma once

#include <memory>
#include <string>
#include "DataKeeper.hxx"

class DataProcessor
{
private:
	DataKeeper const& m_data;

	bool m_energy_ready;
	bool m_poynting_ready;

	std::unique_ptr<std::vector<ScalarField>> ptr_energy;
	std::unique_ptr<std::vector<VectorField>> ptr_poynting;
	std::unique_ptr<std::vector<ScalarField>> ptr_poynting_norm;

	void compute_energy();
	void compute_poynting();
	void compute_all();

public:
	DataProcessor(DataKeeper const& data);
	~DataProcessor();

	ScalarField const& get_energy_density(size_t step);
	VectorField const& get_poynting_vector(size_t step);
	ScalarField const& get_poynting_vector_norm(size_t step);
};