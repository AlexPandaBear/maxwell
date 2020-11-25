#include "ElectrostaticKernel.hxx"

ElectrostaticKernel::ElectrostaticKernel(Mesh3D& mesh, ElectrostaticDataKeeper& data) :
	m_mesh(mesh),
	m_data(data),
	mat_A(data.get_nb_nodes()),
	vec_X(std::vector<double>(data.get_nb_nodes(), 0.)),
	vec_B(std::vector<double>(data.get_nb_nodes(), 0.)) {}

ElectrostaticKernel::~ElectrostaticKernel() {}

void ElectrostaticKernel::build_matrix_A()
{
	size_t node_ids[4];
	const Cell* ptr_cell;

	for (size_t k = 0; k < m_mesh.get_nb_cells(); k++)
	{
		ptr_cell = m_mesh.get_cell_ptr(k);

		for (size_t n = 0; n < 4; n++)
		{
			node_ids[n] = ptr_cell->get_global_node_id(n);
		}

		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 4; j++)
			{
				mat_A.set_coef(node_ids[i], node_ids[j], ptr_cell->compute_int3d_grad_phi_grad_psi(i, j));
			}
		}
	}
}

void ElectrostaticKernel::build_vector_B()
{
	size_t node_ids[4];
	const Cell* ptr_cell;

	for (size_t k = 0; k < m_mesh.get_nb_cells(); k++)
	{
		ptr_cell = m_mesh.get_cell_ptr(k);

		for (size_t n = 0; n < 4; n++)
		{
			node_ids[n] = ptr_cell->get_global_node_id(n);
		}

		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 4; j++)
			{
				vec_B[node_ids[j]] += ptr_cell->compute_int3d_phi_psi(i, j) * m_data.get_rho(node_ids[i])/m_data.get_epsilon(node_ids[i]);
				vec_B[node_ids[j]] -= Vec3D::dot_product(ptr_cell->compute_int2d_phi_psi_n(i, j), m_data.get_E(node_ids[i]));
			}
		}
	}
}

void ElectrostaticKernel::simulate_phi()
{
	const size_t max_nb_iterations(m_data.get_max_nb_iterations());
	const double accuracy(m_data.get_accuracy());

	build_matrix_A();
	build_vector_B();

	double residual;
	size_t it(0);

	do
	{
		residual = mat_A.perform_gauss_seidel_iteration(vec_X, vec_B);
		it++;
	}
	while (residual > accuracy and it < max_nb_iterations);

	if (it == max_nb_iterations and residual > accuracy)
	{
		throw std::runtime_error("Maximum number of iterations reached before convergence");
	}
}

void ElectrostaticKernel::derive_E()
{
	size_t node_ids[4];
	const Cell* ptr_cell;
	Vec3D E_cell;

	for (size_t k = 0; k < m_mesh.get_nb_cells(); k++)
	{
		ptr_cell = m_mesh.get_cell_ptr(k);
		E_cell.set(0., 0., 0.);

		for (size_t n = 0; n < 4; n++)
		{
			node_ids[n] = ptr_cell->get_global_node_id(n);
			E_cell += ptr_cell->get_surface(n) * vec_X[node_ids[n]];
		}

		E_cell /= 3.;

		for (size_t n = 0; n < 4; n++)
		{
			m_data.set_E(node_ids[n], m_data.get_E(node_ids[n]) + E_cell);
		}
	}

	std::vector<size_t> neighbor_cells_id;
	double volume;

	for (size_t n = 0; n < m_data.get_nb_nodes(); n++)
	{
		neighbor_cells_id = m_mesh.get_neighbor_cells_id(n);
		volume = 0.;

		for (size_t k : neighbor_cells_id)
		{
			volume += m_mesh.get_cell_ptr(k)->get_volume();
		}

		m_data.set_E(n, m_data.get_E(n)/volume);
	}
}