#include "ElectrostaticDataKeeper.hxx"

ElectrostaticDataKeeper::ElectrostaticDataKeeper():
	m_eps0(8.85418782*pow(10., -12)),
	m_nb_nodes(0),
	m_accuracy(0.01),
	m_max_nb_iterations(100),
	m_rho(ScalarField(0)),
	m_eps(ScalarField(0)),
	m_E(VectorField(0)) {}

ElectrostaticDataKeeper::~ElectrostaticDataKeeper() {}

size_t ElectrostaticDataKeeper::get_nb_nodes() const
{
	return m_nb_nodes;
}

void ElectrostaticDataKeeper::set_nb_nodes(size_t nb_nodes)
{
	m_nb_nodes = nb_nodes;
	m_rho = ScalarField(nb_nodes);
	m_eps = ScalarField(nb_nodes, m_eps0);
	m_E = VectorField(nb_nodes);
}

double ElectrostaticDataKeeper::get_epsilon0() const
{
	return m_eps0;
}

void ElectrostaticDataKeeper::set_epsilon0(double epsilon0)
{
	m_eps0 = epsilon0;
}

size_t ElectrostaticDataKeeper::get_max_nb_iterations() const
{
	return m_max_nb_iterations;
}

void ElectrostaticDataKeeper::set_max_nb_iterations(size_t max_nb_iterations)
{
	m_max_nb_iterations = max_nb_iterations;
}

double ElectrostaticDataKeeper::get_accuracy() const
{
	return m_accuracy;
}

void ElectrostaticDataKeeper::set_accuracy(double accuracy)
{
	m_accuracy = accuracy;
}

double ElectrostaticDataKeeper::get_rho(size_t node_nb) const
{
	return m_rho.get_value(node_nb);
}

double ElectrostaticDataKeeper::get_epsilon(size_t node_nb) const
{
	return m_eps.get_value(node_nb);
}

Vec3D ElectrostaticDataKeeper::get_E(size_t node_nb) const
{
	return m_E.get_value(node_nb);
}

ScalarField const& ElectrostaticDataKeeper::get_rho() const
{
	return m_rho;
}

ScalarField const& ElectrostaticDataKeeper::get_epsilon() const
{
	return m_eps;
}

VectorField const& ElectrostaticDataKeeper::get_E() const
{
	return m_E;
}

void ElectrostaticDataKeeper::set_rho(size_t node_nb, double rho)
{
	m_rho.set_value(node_nb, rho);
}

void ElectrostaticDataKeeper::set_epsilon_r(size_t node_nb, double epsilon_r)
{
	m_eps.set_value(node_nb, epsilon_r*m_eps0);
}

void ElectrostaticDataKeeper::set_E(size_t node_nb, Vec3D const& E)
{
	m_E.set_value(node_nb, E);
}