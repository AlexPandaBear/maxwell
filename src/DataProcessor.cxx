#include "DataProcessor.hxx"

DataProcessor::DataProcessor(DataKeeper const& data) :
	m_data(data),
	m_energy_ready(false),
	m_poynting_ready(false),
	ptr_energy(new std::vector<ScalarField>(0, ScalarField(0))),
	ptr_poynting(new std::vector<VectorField>(0, VectorField(0))),
	ptr_poynting_norm(new std::vector<ScalarField>(0, ScalarField(0))) {}

DataProcessor::~DataProcessor() {}

void DataProcessor::compute_energy()
{
	size_t nb_steps(m_data.get_nb_steps());
	size_t nb_nodes(m_data.get_nb_nodes());

	double epsilon(m_data.get_epsilon());
	double inv_mu(1. / m_data.get_mu());

	ptr_energy.reset(new std::vector<ScalarField>(nb_steps+1, ScalarField(nb_nodes)));

	for (size_t t = 0; t < nb_steps+1; t++)
	{
		for (size_t n = 0; n < nb_nodes; n++)
		{
			(*ptr_energy)[t].set_value(n, 0.5*(epsilon * m_data.get_E(t, n).compute_squared_norm() + inv_mu * m_data.get_B(t, n).compute_squared_norm()));
		}
	}

	m_energy_ready = true;
}

void DataProcessor::compute_poynting()
{
	size_t nb_steps(m_data.get_nb_steps());
	size_t nb_nodes(m_data.get_nb_nodes());

	double inv_mu(1. / m_data.get_mu());

	ptr_poynting.reset(new std::vector<VectorField>(nb_steps+1, VectorField(nb_nodes)));
	ptr_poynting_norm.reset(new std::vector<ScalarField>(nb_steps+1, ScalarField(nb_nodes)));

	Vec3D Pi;

	for (size_t t = 0; t < nb_steps+1; t++)
	{
		for (size_t n = 0; n < nb_nodes; n++)
		{
			Pi = Vec3D::cross_product(m_data.get_E(t, n), m_data.get_B(t, n)) * inv_mu;
			(*ptr_poynting)[t].set_value(n, Pi);
			(*ptr_poynting_norm)[t].set_value(n, Pi.compute_norm());
		}
	}

	m_poynting_ready = true;
}

void DataProcessor::compute_all()
{
	compute_energy();
}

ScalarField const& DataProcessor::get_energy_density(size_t step)
{
	if (!m_energy_ready)
	{
		compute_energy();
	}

	return (*ptr_energy)[step];
}

VectorField const& DataProcessor::get_poynting_vector(size_t step)
{
	if (!m_poynting_ready)
	{
		compute_poynting();
	}

	return (*ptr_poynting)[step];
}

ScalarField const& DataProcessor::get_poynting_vector_norm(size_t step)
{
	if (!m_poynting_ready)
	{
		compute_poynting();
	}

	return (*ptr_poynting_norm)[step];
}