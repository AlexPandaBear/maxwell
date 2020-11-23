#include "DataProcessor.hxx"

DataProcessor::DataProcessor(DataKeeper& data) :
	m_data(data),
	m_energy_ready(false),
	m_poynting_ready(false),
	m_energy(std::vector<ScalarField>(0, ScalarField(0))),
	m_poynting(std::vector<VectorField>(0, VectorField(0))),
	m_poynting_norm(std::vector<ScalarField>(0, ScalarField(0))) {}

DataProcessor::~DataProcessor() {}

void DataProcessor::compute_energy()
{
	size_t nb_steps(m_data.get_nb_steps());
	size_t nb_nodes(m_data.get_nb_nodes());

	double epsilon(m_data.get_epsilon());
	double inv_mu(1. / m_data.get_mu());

	m_energy = std::vector<ScalarField>(nb_steps+1, ScalarField(nb_nodes));

	for (size_t t = 0; t < nb_steps+1; t++)
	{
		for (size_t n = 0; n < nb_nodes; n++)
		{
			m_energy[t].set_value(n, 0.5*(epsilon * m_data.get_E(t, n).compute_squared_norm() + inv_mu * m_data.get_B(t, n).compute_squared_norm()));
		}
	}

	m_energy_ready = true;
}

void DataProcessor::compute_poynting()
{
	size_t nb_steps(m_data.get_nb_steps());
	size_t nb_nodes(m_data.get_nb_nodes());

	double inv_mu(1. / m_data.get_mu());

	m_poynting = std::vector<VectorField>(nb_steps+1, VectorField(nb_nodes));
	m_poynting_norm = std::vector<ScalarField>(nb_steps+1, ScalarField(nb_nodes));

	Vec3D Pi;

	for (size_t t = 0; t < nb_steps+1; t++)
	{
		for (size_t n = 0; n < nb_nodes; n++)
		{
			Pi = Vec3D::cross_product(m_data.get_E(t, n), m_data.get_B(t, n)) * inv_mu;
			m_poynting[t].set_value(n, Pi);
			m_poynting_norm[t].set_value(n, Pi.compute_norm());
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

	return m_energy[step];
}

VectorField const& DataProcessor::get_poynting_vector(size_t step)
{
	if (!m_poynting_ready)
	{
		compute_poynting();
	}

	return m_poynting[step];
}

ScalarField const& DataProcessor::get_poynting_vector_norm(size_t step)
{
	if (!m_poynting_ready)
	{
		compute_poynting();
	}

	return m_poynting_norm[step];
}