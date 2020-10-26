#include "DataKeeper.hxx"

DataKeeper::DataKeeper() {}

DataKeeper::DataKeeper(std::string file) {} //TODO

DataKeeper::~DataKeeper() {}

double DataKeeper::get_epsilon() const
{
	return m_eps;
}

double DataKeeper::get_mu() const
{
	return m_mu;
}

double DataKeeper::get_t_max() const
{
	return m_t_max;
}

double DataKeeper::get_theta() const
{
	return m_theta;
}

double DataKeeper::get_accuracy() const
{
	return m_accuracy;
}

size_t DataKeeper::get_max_nb_iterations() const
{
	return m_max_nb_iterations;
}

size_t DataKeeper::get_nb_steps() const
{
	return m_nb_steps;
}

size_t DataKeeper::get_nb_nodes() const
{
	return m_nb_nodes;
}

size_t DataKeeper::get_nb_cells() const
{
	return m_nb_cells;
}

void DataKeeper::reset_dimensions(size_t nb_steps, size_t nb_nodes, size_t nb_cells)
{
	m_nb_steps = nb_steps;
	m_nb_nodes = nb_nodes;
	m_nb_cells = nb_cells;

	m_rho = std::vector<Field<double>>(nb_steps+1, Field<double>(nb_cells));
	m_j = std::vector<Field<Vec3D>>(nb_steps+1, Field<Vec3D>(nb_cells));
	m_E = std::vector<Field<Vec3D>>(nb_steps+1, Field<Vec3D>(nb_nodes));
	m_B = std::vector<Field<Vec3D>>(nb_steps+1, Field<Vec3D>(nb_nodes));

	erase_BCs();
}

void DataKeeper::set_epsilon(double epsilon)
{
	m_eps = epsilon;
}

void DataKeeper::set_mu(double mu)
{
	m_mu = mu;
}

void DataKeeper::set_t_max(double t_max)
{
	m_t_max = t_max;
}

void DataKeeper::set_theta(double theta)
{
	m_theta = theta;
}

void DataKeeper::set_accuracy(double accuracy)
{
	m_accuracy = accuracy;
}

void DataKeeper::set_max_nb_iterations(size_t max_nb_iterations)
{
	m_max_nb_iterations = max_nb_iterations;
}

double DataKeeper::get_rho(size_t t, size_t node_nb) const
{
	return m_rho[t].get_value(node_nb);
}

void DataKeeper::set_rho(size_t t, size_t node_nb, double rho)
{
	m_rho[t].set_value(node_nb, rho);
}

Vec3D DataKeeper::get_j(size_t t, size_t node_nb) const
{
	return m_j[t].get_value(node_nb);
}

void DataKeeper::set_j(size_t t, size_t node_nb, Vec3D j)
{
	m_j[t].set_value(node_nb, j);
}

Vec3D DataKeeper::get_E(size_t t, size_t node_nb) const
{
	return m_E[t].get_value(node_nb);
}

void DataKeeper::set_E(size_t t, size_t node_nb, Vec3D E)
{
	m_E[t].set_value(node_nb, E);
}

Vec3D DataKeeper::get_B(size_t t, size_t node_nb) const
{
	return m_B[t].get_value(node_nb);
}

void DataKeeper::set_B(size_t t, size_t node_nb, Vec3D B)
{
	m_B[t].set_value(node_nb, B);
}

void DataKeeper::erase_BCs()
{
	m_BC = std::vector<BoundaryCondition>();
}

void DataKeeper::add_BC(size_t node_nb, Vec3D E, Vec3D B)
{
	m_BC.push_back(BoundaryCondition(node_nb, E, B));
}

Vec3D DataKeeper::get_boundary_condition_E(size_t node_nb) const
{
	for (size_t i = 0; i < m_BC.size(); i++)
	{
		if (m_BC[i].get_node_nb() == node_nb)
		{
			return m_BC[i].get_E();
		}
	}

	throw std::invalid_argument("No boundary condition defined for node " + node_nb);
}

Vec3D DataKeeper::get_boundary_condition_B(size_t node_nb) const
{
	for (size_t i = 0; i < m_BC.size(); i++)
	{
		if (m_BC[i].get_node_nb() == node_nb)
		{
			return m_BC[i].get_B();
		}
	}

	throw std::invalid_argument("No boundary condition defined for node " + node_nb);
}

bool DataKeeper::boundary_node(size_t node_nb) const
{
	for (size_t i = 0; i < m_BC.size(); i++)
	{
		if (m_BC[i].get_node_nb() == node_nb)
		{
			return true;
		}
	}

	return false;
}

bool DataKeeper::check_ID_BC_compatibility() const
{
	for (size_t n = 0; n < m_BC.size();n++)
	{
		if (m_E[0].get_value(n) != m_BC[n].get_E() || m_B[0].get_value(n) != m_BC[n].get_B())
		{
			return false;
		}
	}

	return true;
}

void DataKeeper::save(std::string file) const {} //TODO