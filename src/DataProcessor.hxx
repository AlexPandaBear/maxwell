#pragma once

#include <string>
#include "DataKeeper.hxx"
#include "UnsteadyField.hxx"

class DataProcessor
{
private:
	DataKeeper& m_data;

	bool m_energy_ready;
	std::vector<ScalarField> m_energy;

	void compute_energy();
	void compute_all();

public:
	DataProcessor(DataKeeper& data);
	~DataProcessor();

	ScalarField const& get_energy_density(size_t step);
};