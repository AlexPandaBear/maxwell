#include "SimManager.hxx"

SimManager::SimManager() :
	m_mesh(),
	m_data(),
	m_processor(m_data) 
{
	std::cout << "Creation of the simulation environment [ OK ]" << std::endl;
}
/*
SimManager::SimManager(std::string simulation_name) :
	m_mesh(Mesh3D(simulation_name + ".m")),
	m_data(DataKeeper(simulation_name + ".d")),
	m_processor(DataProcessor(simulation_name + ".p")) {}
*/
SimManager::~SimManager() {}

void SimManager::set_constants(double epsilon, double mu, double sigma)
{
	std::cout << "Definition of the physical constants ";
	
	m_data.set_epsilon(epsilon);
	m_data.set_mu(mu);
	m_data.set_sigma(sigma);
	
	std::cout << "[ OK ]" << std::endl;
}

void SimManager::set_simulation_parameters(double t_max, size_t nb_steps, double theta, double accuracy, size_t max_nb_iterations, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, size_t nx, size_t ny, size_t nz)
{
	std::cout << "Creation of the mesh and the time samples ";
	
	m_data.set_t_max(t_max);
	m_data.set_theta(theta);
	m_data.set_accuracy(accuracy);
	m_data.set_max_nb_iterations(max_nb_iterations);
	m_mesh.generate_grid_mesh(x_min, x_max, y_min, y_max, z_min, z_max, nx, ny, nz);
	m_data.reset_dimensions(nb_steps, m_mesh.get_nb_nodes(), m_mesh.get_nb_cells());
	m_data.erase_BCs();
	
	std::cout << "[ OK ]" << std::endl;
}

void SimManager::define_initial_state(VectorField const& E0, VectorField const& B0)
{
	m_data.set_E(0, E0);
	m_data.set_B(0, B0);
}

void SimManager::define_initial_state_background_values(double rho0, Vec3D j0, Vec3D E0, Vec3D B0)
{
	std::cout << "Definition of the initial state ";
	
	for (size_t n = 0; n < m_data.get_nb_nodes(); n++)
	{
		m_data.set_rho(0, n, rho0);
		m_data.set_j(0, n, j0);
		m_data.set_E(0, n, E0);
		m_data.set_B(0, n, B0);
	}

	std::cout << "[ OK ]" << std::endl;
}

double SimManager::generate_random_double(double min, double max) const
{
	return (max - min) * ( (double) rand() / (double) RAND_MAX ) + min;
}

double SimManager::compute_distance_to_axis(Vec3D pt, Vec3D axis_pt_A, Vec3D axis_pt_B) const
{
	Vec3D axis(axis_pt_B - axis_pt_A);
	axis.normalize();
	Vec3D proj(axis * Vec3D::dot_product(pt - axis_pt_A, axis));
	return (pt - axis_pt_A - proj).compute_norm();
}

double SimManager::estimate_wire_intersection(Vec3D pt_A, Vec3D pt_B, double radius, size_t cell_nb, size_t sample_size) const
{
	double coef(0.);
	double alpha, beta, gamma, delta;
	Vec3D random_pt;
	const Cell* ptr_cell(m_mesh.get_cell_ptr(cell_nb)); 

	for (size_t i = 0; i < sample_size; i++)
	{
		alpha = generate_random_double(0., 1.);
		beta = generate_random_double(0., 1. - alpha);
		gamma = generate_random_double(0., 1. - (alpha + beta));
		delta = 1. - (alpha + beta + gamma);

		random_pt = ptr_cell->get_node_xyz(0)*alpha + ptr_cell->get_node_xyz(1)*beta + ptr_cell->get_node_xyz(2)*gamma + ptr_cell->get_node_xyz(3)*delta;

		if (compute_distance_to_axis(random_pt, pt_A, pt_B) > radius)
		{
			coef += 1.;
		}
	}

	return coef/sample_size;
}

void SimManager::add_wire(std::vector<Vec3D> wire_skeleton, double wire_radius, double wire_current)
{
	Vec3D vec_j;
	double j(wire_current / (M_PI*wire_radius*wire_radius));

	double coef;

	for (size_t s = 0; s < wire_skeleton.size()-1; s++)
	{
		vec_j = wire_skeleton[s+1] - wire_skeleton[s];
		vec_j.normalize();
		vec_j *= j;

		for (size_t c = 0; c < m_data.get_nb_cells(); c++)
		{
			coef = estimate_wire_intersection(wire_skeleton[s+1], wire_skeleton[s], wire_radius, c, 100);

			if (coef > 0.)
			{
				m_data.set_j(0, c, vec_j*coef);
			}
		}
	}
}

void SimManager::add_boundary_condition(size_t node_nb)
{
	m_data.add_BC(node_nb);
}

void SimManager::lock_all_boundary_nodes()
{
	std::cout << "Locking as constant Dirichlet boundary conditions all the frontier ";

	for (size_t n : m_mesh.get_boundary_nodes())
	{
		m_data.add_BC(n);
	}

	std::cout << "[ OK ]" << std::endl;
}

void SimManager::simulate()
{
	std::cout << "Creation of the kernel ";
	UnsteadyMaxwellKernel kernel(m_mesh, m_data);
	std::cout << "[ OK ]" << std::endl;

	kernel.simulate();
}

/*
void SimManager::save(std::string simulation_name) const
{
	m_mesh.save(simulation_name + ".m");
	m_data.save(simulation_name + ".d");
}
*/

size_t SimManager::get_nb_nodes() const
{
	return m_mesh.get_nb_nodes();
}

size_t SimManager::get_nb_cells() const
{
	return m_mesh.get_nb_cells();
}

std::vector<size_t> SimManager::get_node_ids(size_t cell_id) const
{
	Cell cell = m_mesh.get_cell(cell_id);
	return {cell.get_global_node_id(0), cell.get_global_node_id(1), cell.get_global_node_id(2), cell.get_global_node_id(3)};
}

Vec3D SimManager::get_node_xyz(size_t node_id) const
{
	return m_mesh.get_node_xyz(node_id);
}

VectorField const& SimManager::get_mesh()
{
	return m_mesh.get_all_nodes_xyz();
}

std::vector<double> const& SimManager::get_time()
{
	return m_data.get_time();
}

ScalarField const& SimManager::get_rho(size_t step)
{
	return m_data.get_rho(step);
}

VectorField const& SimManager::get_j(size_t step)
{
	return m_data.get_j(step);
}

VectorField const& SimManager::get_E(size_t step)
{
	return m_data.get_E(step);
}

VectorField const& SimManager::get_B(size_t step)
{
	return m_data.get_B(step);
}

ScalarField const& SimManager::get_energy_density(size_t step)
{
	return m_processor.get_energy_density(step);
}

VectorField const& SimManager::get_poynting_vector(size_t step)
{
	return m_processor.get_poynting_vector(step);
}

ScalarField const& SimManager::get_poynting_vector_norm(size_t step)
{
	return m_processor.get_poynting_vector_norm(step);
}