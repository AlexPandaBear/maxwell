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
	mat_M(10*m_nb_nodes),
	mat_A(10*m_nb_nodes),
	mat_L(10*m_nb_nodes),
	mat_R(10*m_nb_nodes),
	vec_X(10*m_nb_nodes, 0.),
	vec_B(10*m_nb_nodes, 0.) {}

UnsteadyMaxwellKernel::~UnsteadyMaxwellKernel() {}

bool UnsteadyMaxwellKernel::boundary_node(size_t node_nb) const
{
	return m_data.boundary_node(node_nb);
}

bool UnsteadyMaxwellKernel::check_divergence_in_ID() const
{
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
	std::cout << "Building of the mass and stiffness matrices ";

	size_t id_nodes[4];
	const Cell* ptr_cell;
	size_t i_block, j_block;

	Vec3D S, N;
	double alpha, Sx, Sy, Sz, Nx, Ny, Nz, c2, epsilon, sigma;

	for (size_t k = 0; k < m_nb_cells; k++)
	{
		ptr_cell = m_mesh.get_cell_ptr(k);

		for (int i = 0; i < 4; i++)
		{
			id_nodes[i] = ptr_cell->get_global_node_id(i);
		}



		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 4; j++)
			{
				i_block = id_nodes[i];
				j_block = id_nodes[j];



				alpha = 0.2 * ptr_cell->get_volume();
				
				S = Vec3D(0., 0., 0.);

				for (int n = 0; n < 4; n++)
				{
					S += ptr_cell->get_surface(n);
				}

				S -= ptr_cell->get_surface(i);

				if (j != i)
				{
					S -= ptr_cell->get_surface(j);
					S /= 2.;

					alpha /= 2.;
				}

				S /= 6.;

				Sx = S.get_x();
				Sy = S.get_y();
				Sz = S.get_z();
				
				N = ptr_cell->get_surface(i) / (-12.);

				Nx = N.get_x();
				Ny = N.get_y();
				Nz = N.get_z();

				c2 = 1./(m_data.get_epsilon()*m_data.get_mu());
				epsilon = m_data.get_epsilon();
				sigma = m_data.get_sigma();



				mat_M.add_to_coef(10*i_block, 10*j_block, alpha);
				mat_A.add_to_coef(10*i_block, 10*j_block + 1, Sx-Nx);
				mat_A.add_to_coef(10*i_block, 10*j_block + 2, Sy-Ny);
				mat_A.add_to_coef(10*i_block, 10*j_block + 3, Sz-Nz);

				mat_A.add_to_coef(10*i_block + 1, 10*j_block + 1, alpha);
				mat_A.add_to_coef(10*i_block + 1, 10*j_block + 4, -sigma*alpha);

				mat_A.add_to_coef(10*i_block + 2, 10*j_block + 2, alpha);
				mat_A.add_to_coef(10*i_block + 2, 10*j_block + 5, -sigma*alpha);

				mat_A.add_to_coef(10*i_block + 3, 10*j_block + 3, alpha);
				mat_A.add_to_coef(10*i_block + 3, 10*j_block + 6, -sigma*alpha);

				mat_M.add_to_coef(10*i_block + 4, 10*j_block + 4, alpha);
				mat_A.add_to_coef(10*i_block + 4, 10*j_block + 1, alpha/epsilon);
				mat_A.add_to_coef(10*i_block + 4, 10*j_block + 8, -c2*(Sz+Nz));
				mat_A.add_to_coef(10*i_block + 4, 10*j_block + 9, c2*(Sy+Ny));

				mat_M.add_to_coef(10*i_block + 5, 10*j_block + 5, alpha);
				mat_A.add_to_coef(10*i_block + 5, 10*j_block + 2, alpha/epsilon);
				mat_A.add_to_coef(10*i_block + 5, 10*j_block + 7, c2*(Sz+Nz));
				mat_A.add_to_coef(10*i_block + 5, 10*j_block + 9, -c2*(Sx+Nx));

				mat_M.add_to_coef(10*i_block + 6, 10*j_block + 6, alpha);
				mat_A.add_to_coef(10*i_block + 6, 10*j_block + 3, alpha/epsilon);
				mat_A.add_to_coef(10*i_block + 6, 10*j_block + 7, -c2*(Sy+Ny));
				mat_A.add_to_coef(10*i_block + 6, 10*j_block + 8, c2*(Sx+Nx));

				mat_M.add_to_coef(10*i_block + 7, 10*j_block + 7, alpha);
				mat_A.add_to_coef(10*i_block + 7, 10*j_block + 5, Sz+Nz);
				mat_A.add_to_coef(10*i_block + 7, 10*j_block + 6, -(Sy+Ny));

				mat_M.add_to_coef(10*i_block + 8, 10*j_block + 8, alpha);
				mat_A.add_to_coef(10*i_block + 8, 10*j_block + 4, -(Sz+Nz));
				mat_A.add_to_coef(10*i_block + 8, 10*j_block + 6, Sx+Nx);

				mat_M.add_to_coef(10*i_block + 9, 10*j_block + 9, alpha);
				mat_A.add_to_coef(10*i_block + 9, 10*j_block + 4, Sy+Ny);
				mat_A.add_to_coef(10*i_block + 9, 10*j_block + 5, -(Sx+Nx));
			}
		}
	}

	std::cout << "[ OK ]" << std::endl;

	std::cout << "Building of the time stepping matrices ";
	
	mat_L = mat_M + mat_A*m_theta*m_dt;
	mat_R = mat_M - mat_A*(1-m_theta)*m_dt;

	for (size_t n : m_data.get_dirichlet_nodes())
	{
		for (int i = 0; i < 10; i++)
		{
			mat_L.set_boundary_equation(10*n + i);
			mat_R.set_boundary_equation(10*n + i);
		}
	}
	
	std::cout << "[ OK ]" << std::endl;
}

void UnsteadyMaxwellKernel::initialize_vector_X()
{
	Vec3D tmp;

	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		vec_X[10*n] = m_data.get_rho(0, n);
		
		tmp = m_data.get_j(0, n);
		vec_X[10*n + 1] = tmp.get_x();
		vec_X[10*n + 2] = tmp.get_y();
		vec_X[10*n + 3] = tmp.get_z();

		tmp = m_data.get_E(0, n);
		vec_X[10*n + 4] = tmp.get_x();
		vec_X[10*n + 5] = tmp.get_y();
		vec_X[10*n + 6] = tmp.get_z();

		tmp = m_data.get_B(0, n);
		vec_X[10*n + 7] = tmp.get_x();
		vec_X[10*n + 8] = tmp.get_y();
		vec_X[10*n + 9] = tmp.get_z();
	}
}

void UnsteadyMaxwellKernel::build_vector_B(size_t t)
{
	double sum;
	std::vector<std::pair<size_t, double>> R_line;

	for (size_t i = 0; i < 10*m_nb_nodes; i++)
	{
		sum = 0.;
		R_line = mat_R.get_reduced_line_w_diag(i);

		for (size_t j = 0; j < R_line.size(); j++)
		{
			sum += R_line[j].second * vec_X[R_line[j].first];
		}

		vec_B[i] = sum;
	}
}

double UnsteadyMaxwellKernel::GS_step(size_t i) const
{
	double new_Xi(vec_B[i]);
	std::vector<std::pair<size_t, double>> L_line(mat_L.get_reduced_line_wo_diag(i));

	for (size_t j = 0; j < L_line.size(); j++)
	{
		new_Xi -= L_line[j].second * vec_X[L_line[j].first];
	}

	return new_Xi/mat_L.get_coef(i, i);
}

double UnsteadyMaxwellKernel::GS_iteration()
{
	double new_Xi;
	double energy_diff(0.);

	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		for (int i = 0; i < 4; i++) //updating rho and j
		{
			vec_X[10*n + i] = GS_step(10*n + i);
		}

		for (size_t i = 4; i < 7; i++) //updating E
		{
			new_Xi = GS_step(10*n + i);
			energy_diff += m_eps*(new_Xi*new_Xi - vec_X[10*n + i]*vec_X[10*n + i]);
			vec_X[10*n + i] = new_Xi;
		}

		for (size_t i = 4; i < 9; i++) //updating B
		{
			new_Xi = GS_step(10*n + i);
			energy_diff += (new_Xi*new_Xi - vec_X[10*n + i]*vec_X[10*n + i])/m_mu;
			vec_X[10*n + i] = new_Xi;
		}
	}

	return 0.5*energy_diff/m_nb_nodes;
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
		build_vector_B(t);
		k = 0;

		do
		{
			mean_energy_diff = GS_iteration();
			k += 1;
		}
		while (mean_energy_diff > m_accuracy && k < m_max_nb_iterations);

		if (k == m_max_nb_iterations)
		{
			std::cout << "[ FAILED ]" << std::endl;
			throw std::runtime_error("Maximun number of iterations reached while computing step " + t+1);
		}

		for (size_t n = 0; n < m_nb_nodes; n++)
		{
			m_data.set_rho(t+1, n, vec_X[10*n]);
			m_data.set_j(t+1, n, Vec3D(vec_X[10*n + 1], vec_X[10*n + 2], vec_X[10*n + 3]));
			m_data.set_E(t+1, n, Vec3D(vec_X[10*n + 4], vec_X[10*n + 5], vec_X[10*n + 6]));
			m_data.set_B(t+1, n, Vec3D(vec_X[10*n + 7], vec_X[10*n + 8], vec_X[10*n + 9]));
		}
	}

	std::cout << "[ OK ]" << std::endl;
}