#include "SimManager.hxx"

SimManager::SimManager() {}

SimManager::SimManager(std::string simulation_name) :
	m_mesh(Mesh3D(simulation_name + ".m")),
	m_data(DataKeeper(simulation_name + ".d")) {}

SimManager::~SimManager() {}

void SimManager::set_constants(double epsilon, double mu)
{
	m_data.set_epsilon(epsilon);
	m_data.set_mu(mu);
}

void SimManager::set_simulation_parameters(double t_max, size_t nb_steps, double theta, double accuracy, size_t max_nb_iterations, size_t nb_pts, size_t nb_cells)
{
	m_data.set_t_max(t_max);
	m_data.set_theta(theta);
	m_data.set_accuracy(accuracy);
	m_data.set_max_nb_iterations(max_nb_iterations);
	m_data.reset_dimensions(nb_steps, nb_pts, nb_cells);
	m_data.erase_BCs();
}

void SimManager::define_initial_state(std::vector<double> rho0, std::vector<Vec3D> j0, std::vector<Vec3D> E0, std::vector<Vec3D> B0)
{
	for (size_t c = 0; c < m_data.get_nb_cells(); c++)
	{
		m_data.set_rho(0, c, rho0[c]);
		m_data.set_j(0, c, j0[c]);
	}

	for (size_t n = 0; n < m_data.get_nb_nodes(); n++)
	{
		m_data.set_E(0, n, E0[n]);
		m_data.set_B(0, n, B0[n]);
	}
}

void SimManager::add_boundary_condition(size_t node_nb, Vec3D E, Vec3D B)
{
	m_data.add_BC(node_nb, E, B);
}

void SimManager::simulate()
{
	UnsteadyMaxwellKernel kernel(m_mesh, m_data);
	kernel.simulate();
}

void SimManager::save(std::string simulation_name) const
{
	m_mesh.save(simulation_name + ".m");
	m_data.save(simulation_name + ".d");
}