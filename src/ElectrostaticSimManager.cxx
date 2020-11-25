#include "ElectrostaticSimManager.hxx"

ElectrostaticSimManager::ElectrostaticSimManager(bool verbose) :
	m_verbose(verbose)
{
	if (verbose)
	{
		std::cout << "Creating an electrostatic simulation environment [ OK ]" << std::endl;
	}
}

ElectrostaticSimManager::~ElectrostaticSimManager() {}

void ElectrostaticSimManager::set_epsilon0(double epsilon0)
{
	m_data.set_epsilon0(epsilon0);
}

void ElectrostaticSimManager::set_accuracy(double accuracy)
{
	m_data.set_accuracy(accuracy);
}

void ElectrostaticSimManager::set_max_nb_iterations(size_t max_nb_iterarions)
{
	m_data.set_max_nb_iterations(max_nb_iterarions);
}

void ElectrostaticSimManager::generate_cube_mesh(double x_min, double x_max, size_t nx,
												double y_min, double y_max, size_t ny,
												double z_min, double z_max, size_t nz)
{
	if (m_verbose)
	{
		std::cout << "Building the mesh ";
	}

	m_data.set_nb_nodes(nx*ny*nz);
	m_mesh.generate_grid_mesh(x_min, x_max, y_min, y_max, z_min, z_max, nx, ny, nz);
	
	if (m_verbose)
	{
		std::cout << "[ OK ]" << std::endl;
	}
}

void ElectrostaticSimManager::set_rho(size_t node_id, double rho)
{
	if (node_id >= m_data.get_nb_nodes())
	{
		throw std::invalid_argument("Wrong node number");
	}

	m_data.set_rho(node_id, rho);
}

void ElectrostaticSimManager::set_rho_field(ScalarField const& rho)
{
	if (rho.get_nb_nodes() != m_data.get_nb_nodes())
	{
		throw std::invalid_argument("Wrong number of nodes in this field");
	}

	for (size_t n = 0; n < m_data.get_nb_nodes(); n++)
	{
		m_data.set_rho(n, rho.get_value(n));
	}
}

void ElectrostaticSimManager::set_epsilon_r_field(ScalarField const& epsilon_r)
{
	if (epsilon_r.get_nb_nodes() != m_data.get_nb_nodes())
	{
		throw std::invalid_argument("Wrong number of nodes in this field");
	}

	for (size_t n = 0; n < m_data.get_nb_nodes(); n++)
	{
		m_data.set_epsilon_r(n, epsilon_r.get_value(n));
	}
}

void ElectrostaticSimManager::simulate()
{
	ElectrostaticKernel kernel(m_mesh, m_data, m_verbose);
	kernel.simulate_phi();
	kernel.derive_E();
}

size_t ElectrostaticSimManager::get_nb_nodes() const
{
	return m_data.get_nb_nodes();
}

Vec3D ElectrostaticSimManager::get_node_xyz(size_t node_id) const
{
	return m_mesh.get_node_xyz(node_id);
}

VectorField const& ElectrostaticSimManager::get_mesh() const
{
	return m_mesh.get_all_nodes_xyz();
}

ScalarField const& ElectrostaticSimManager::get_rho() const
{
	return m_data.get_rho();
}

ScalarField const& ElectrostaticSimManager::get_epsilon() const
{
	return m_data.get_epsilon();
}

VectorField const& ElectrostaticSimManager::get_E() const
{
	return m_data.get_E();
}
