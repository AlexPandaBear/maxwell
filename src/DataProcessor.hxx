#pragma once

#include <string>
#include "DataKeeper.hxx"
#include "UnsteadyField.hxx"

class DataProcessor
{
private:
	DataKeeper& m_data;

	bool m_energy_ready;
	UnsteadyField<double> m_energy;

	void compute_energy();
	void compute_all();

public:
	DataProcessor(DataKeeper& data);
	~DataProcessor();

	UnsteadyField<double> const& get_energy_density();
};