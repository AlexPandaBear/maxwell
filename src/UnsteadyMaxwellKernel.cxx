#include "UnsteadyMaxwellKernel.hxx"

UnsteadyMaxwellKernel::UnsteadyMaxwellKernel(Mesh3D& mesh, DataKeeper& data) :
	m_eps(data.get_epsilon()),
	m_mu(data.get_mu()),
	m_sigma(data.get_sigma()),
	m_nb_steps(data.get_nb_steps()),
	m_t_max(data.get_t_max()),
	m_dt(m_t_max/m_nb_steps),
	m_theta(data.get_theta()),
	m_accuracy(data.get_accuracy()),
	m_max_nb_iterations(data.get_max_nb_iterations()),
	m_nb_nodes(data.get_nb_nodes()),
	m_nb_cells(data.get_nb_cells()),
	m_mesh(mesh),
	m_data(data),
	mat_L(m_nb_nodes, 6),
	mat_R(m_nb_nodes, 6),
	vec_X(m_nb_nodes * 6),
	vec_B(m_nb_nodes * 6) {}

UnsteadyMaxwellKernel::~UnsteadyMaxwellKernel() {}

bool UnsteadyMaxwellKernel::check_divergence_in_ID() const
{
	return true;

	double div_E, rho_eps, div_B;
	Vec3D E_sum, B_sum;
	const Cell* ptr_cell;
	size_t id_nodes[4];

	std::cout << "Check of initial data divergence conditions ";

	for (size_t k = 0;  k < m_nb_cells; k++)
	{
		ptr_cell = m_mesh.get_cell_ptr(k);

		for (size_t n = 0; n < 4; n++)
		{
			id_nodes[n] = ptr_cell->get_global_node_id(n);
		}

		div_E = 0.;

		E_sum = m_data.get_E(0, id_nodes[0]) + m_data.get_E(0, id_nodes[1]) + m_data.get_E(0, id_nodes[2]);
		div_E += Vec3D::dot_product(E_sum, ptr_cell->get_surface(3));

		E_sum = m_data.get_E(0, id_nodes[0]) + m_data.get_E(0, id_nodes[1]) + m_data.get_E(0, id_nodes[3]);
		div_E += Vec3D::dot_product(E_sum, ptr_cell->get_surface(2));

		E_sum = m_data.get_E(0, id_nodes[0]) + m_data.get_E(0, id_nodes[2]) + m_data.get_E(0, id_nodes[3]);
		div_E += Vec3D::dot_product(E_sum, ptr_cell->get_surface(1));

		E_sum = m_data.get_E(0, id_nodes[1]) + m_data.get_E(0, id_nodes[2]) + m_data.get_E(0, id_nodes[3]);
		div_E += Vec3D::dot_product(E_sum, ptr_cell->get_surface(0));

		
		rho_eps = 3.*ptr_cell->get_volume() * ( m_data.get_rho(0, id_nodes[0]) + m_data.get_rho(0, id_nodes[1]) + m_data.get_rho(0, id_nodes[2]) + m_data.get_rho(0, id_nodes[3]) ) / (4.*m_eps);
		

		div_B = 0.;

		B_sum = m_data.get_B(0, id_nodes[0]) + m_data.get_B(0, id_nodes[1]) + m_data.get_B(0, id_nodes[2]);
		div_B += Vec3D::dot_product(B_sum, ptr_cell->get_surface(3));

		B_sum = m_data.get_B(0, id_nodes[0]) + m_data.get_B(0, id_nodes[1]) + m_data.get_B(0, id_nodes[3]);
		div_B += Vec3D::dot_product(B_sum, ptr_cell->get_surface(2));

		B_sum = m_data.get_B(0, id_nodes[0]) + m_data.get_B(0, id_nodes[2]) + m_data.get_B(0, id_nodes[3]);
		div_B += Vec3D::dot_product(B_sum, ptr_cell->get_surface(1));

		B_sum = m_data.get_B(0, id_nodes[1]) + m_data.get_B(0, id_nodes[2]) + m_data.get_B(0, id_nodes[3]);
		div_B += Vec3D::dot_product(B_sum, ptr_cell->get_surface(0));

		if ((div_E != rho_eps) || (div_B != 0))
		{
			std::cout << "[ FAILED ]" << std::endl;
			return false;
		}
	}

	std::cout << "[ OK ]" << std::endl;
	return true;
}

void UnsteadyMaxwellKernel::build_matrices()
{
	std::cout << "Building of the time stepping matrices ";

	const Cell* ptr_cell;

	size_t id_nodes[4];
	Matrix block_L(6), block_R(6);
	double A;
	Vec3D G, S, C;

	double theta_dt = m_theta * m_dt;
	double c2_theta_dt = theta_dt / (m_eps * m_mu);
	double one_minus_theta_dt = (1. - m_theta) * m_dt;
	double c2_one_minus_theta_dt = one_minus_theta_dt / (m_eps * m_mu);

	for (size_t k = 0; k < m_nb_cells; k++)
	{
		ptr_cell = m_mesh.get_cell_ptr(k);

		for (size_t i = 0; i < 4; i++)
		{
			id_nodes[i] = ptr_cell->get_global_node_id(i);
		}

		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 4; j++)
			{
				A = ptr_cell->compute_int3d_phi_psi(i, j);
				G = ptr_cell->compute_int3d_phi_grad_psi(i);
				S = ptr_cell->compute_int2d_phi_psi_n(i, j);
				C = G + S;

				block_L.reset_as_diagonal(A);
				block_R.reset_as_diagonal(A);

				for (size_t c = 0; c < 3; c++)
				{
					block_L.multiply_coef(c, c, 1. + m_theta * m_dt * m_mu * m_sigma);
					block_R.multiply_coef(c, c ,1. - (1. - m_theta) * m_dt * m_mu * m_sigma);
				}

				block_L.set_coef(0, 4, - c2_theta_dt * C.get_z());
				block_L.set_coef(0, 5, c2_theta_dt * C.get_y());
				block_L.set_coef(1, 3, c2_theta_dt * C.get_z());
				block_L.set_coef(1, 5, - c2_theta_dt * C.get_x());
				block_L.set_coef(2, 3, - c2_theta_dt * C.get_y());
				block_L.set_coef(2, 4, c2_theta_dt * C.get_x());

				block_L.set_coef(3, 1, - theta_dt * C.get_z());
				block_L.set_coef(3, 2, theta_dt * C.get_y());
				block_L.set_coef(4, 0, theta_dt * C.get_z());
				block_L.set_coef(4, 2, - theta_dt * C.get_x());
				block_L.set_coef(5, 0, - theta_dt * C.get_y());
				block_L.set_coef(5, 1, theta_dt * C.get_x());

				block_R.set_coef(0, 4, - c2_one_minus_theta_dt * C.get_z());
				block_R.set_coef(0, 5, c2_one_minus_theta_dt * C.get_y());
				block_R.set_coef(1, 3, c2_one_minus_theta_dt * C.get_z());
				block_R.set_coef(1, 5, - c2_one_minus_theta_dt * C.get_x());
				block_R.set_coef(2, 3, - c2_one_minus_theta_dt * C.get_y());
				block_R.set_coef(2, 4, c2_one_minus_theta_dt * C.get_x());

				block_R.set_coef(3, 1, - one_minus_theta_dt * C.get_z());
				block_R.set_coef(3, 2, one_minus_theta_dt * C.get_y());
				block_R.set_coef(4, 0, one_minus_theta_dt * C.get_z());
				block_R.set_coef(4, 2, - one_minus_theta_dt * C.get_x());
				block_R.set_coef(5, 0, - one_minus_theta_dt * C.get_y());
				block_R.set_coef(5, 1, one_minus_theta_dt * C.get_x());

				mat_L.add_to_interaction_bloc(id_nodes[i], id_nodes[j], block_L);
				mat_R.add_to_interaction_bloc(id_nodes[i], id_nodes[j], block_R);
			}
		}
	}

	for (size_t n : m_data.get_dirichlet_nodes())
	{
		for (int i = 0; i < 10; i++)
		{
			//TODO
		}
	}
	
	std::cout << "[ OK ]" << std::endl;
}

void UnsteadyMaxwellKernel::initialize_vector_X()
{
	Vec3D tmp;

	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		tmp = m_data.get_E(0, n);
		vec_X[6*n] = tmp.get_x();
		vec_X[6*n + 1] = tmp.get_y();
		vec_X[6*n + 2] = tmp.get_z();

		tmp = m_data.get_B(0, n);
		vec_X[6*n + 3] = tmp.get_x();
		vec_X[6*n + 4] = tmp.get_y();
		vec_X[6*n + 5] = tmp.get_z();
	}
}

void UnsteadyMaxwellKernel::simulate()
{
	if (!check_divergence_in_ID())
	{
		throw std::invalid_argument("Initial data is not solution of Maxwell-Gauss and Maxwell-Thompson equations");
	}

	build_matrices();
	initialize_vector_X();

	double mean_energy_diff;
	size_t k;

	for (size_t t = 0; t < m_nb_steps; t++)
	{
		std::cout << "\rSimulation (step " << t+1 << " / " << m_nb_steps << " -- " << 100*(t+1)/m_nb_steps << " % completed) " << std::flush;
		
		m_data.set_time(t+1, (t+1)*m_dt);
		mat_R.perform_matrix_vector_product(vec_X, vec_B);
		k = 0;

		do
		{
			mean_energy_diff = mat_L.perform_gauss_seidel_iteration(vec_X, vec_B);
			k += 1;
		}
		while (mean_energy_diff > m_accuracy and k < m_max_nb_iterations);

		if (k == m_max_nb_iterations and mean_energy_diff > m_accuracy)
		{
			std::cout << "[ FAILED ]" << std::endl;
			throw std::runtime_error("Maximun number of iterations reached while computing step " + t+1);
		}

		for (size_t n = 0; n < m_nb_nodes; n++)
		{
			m_data.set_E(t+1, n, Vec3D(vec_X[6*n], vec_X[6*n + 1], vec_X[6*n + 2]));
			m_data.set_B(t+1, n, Vec3D(vec_X[6*n + 3], vec_X[6*n + 4], vec_X[6*n + 5]));
		}
	}

	std::cout << "[ OK ]" << std::endl;
}