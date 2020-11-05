#include "UnsteadyMaxwellKernel.hxx"

UnsteadyMaxwellKernel::UnsteadyMaxwellKernel(Mesh3D& mesh, DataKeeper& data) :
	m_eps(data.get_epsilon()),
	m_mu(data.get_mu()),
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
	mat_A(6*m_nb_nodes),
	vec_X(6*m_nb_nodes, 0.),
	vec_B(6*m_nb_nodes, 0.) {}

UnsteadyMaxwellKernel::~UnsteadyMaxwellKernel() {}

bool UnsteadyMaxwellKernel::boundary_node(size_t node_nb) const
{
	return m_data.boundary_node(node_nb);
}

bool UnsteadyMaxwellKernel::check_divergence_in_ID() const
{
	double divE, divB;
	Cell cell(m_mesh.get_cell(0));

	for (size_t c = 0; c < m_nb_cells; c++)
	{
		divE = 0.;
		divB = 0.;
		cell = m_mesh.get_cell(c);

		divE += m_data.get_E(0, cell.get_global_node_id(0)).get_x() * cell.get_ddx_lamb0();
		divE += m_data.get_E(0, cell.get_global_node_id(0)).get_y() * cell.get_ddy_lamb0();
		divE += m_data.get_E(0, cell.get_global_node_id(0)).get_z() * cell.get_ddz_lamb0();

		divE += m_data.get_E(0, cell.get_global_node_id(1)).get_x() * cell.get_ddx_lamb1();
		divE += m_data.get_E(0, cell.get_global_node_id(1)).get_y() * cell.get_ddy_lamb1();
		divE += m_data.get_E(0, cell.get_global_node_id(1)).get_z() * cell.get_ddz_lamb1();

		divE += m_data.get_E(0, cell.get_global_node_id(2)).get_x() * cell.get_ddx_lamb2();
		divE += m_data.get_E(0, cell.get_global_node_id(2)).get_y() * cell.get_ddy_lamb2();
		divE += m_data.get_E(0, cell.get_global_node_id(2)).get_z() * cell.get_ddz_lamb2();

		divE += m_data.get_E(0, cell.get_global_node_id(3)).get_x() * cell.get_ddx_lamb3();
		divE += m_data.get_E(0, cell.get_global_node_id(3)).get_y() * cell.get_ddy_lamb3();
		divE += m_data.get_E(0, cell.get_global_node_id(3)).get_z() * cell.get_ddz_lamb3();


		divB += m_data.get_B(0, cell.get_global_node_id(0)).get_x() * cell.get_ddx_lamb0();
		divB += m_data.get_B(0, cell.get_global_node_id(0)).get_y() * cell.get_ddy_lamb0();
		divB += m_data.get_B(0, cell.get_global_node_id(0)).get_z() * cell.get_ddz_lamb0();

		divB += m_data.get_B(0, cell.get_global_node_id(1)).get_x() * cell.get_ddx_lamb1();
		divB += m_data.get_B(0, cell.get_global_node_id(1)).get_y() * cell.get_ddy_lamb1();
		divB += m_data.get_B(0, cell.get_global_node_id(1)).get_z() * cell.get_ddz_lamb1();

		divB += m_data.get_B(0, cell.get_global_node_id(2)).get_x() * cell.get_ddx_lamb2();
		divB += m_data.get_B(0, cell.get_global_node_id(2)).get_y() * cell.get_ddy_lamb2();
		divB += m_data.get_B(0, cell.get_global_node_id(2)).get_z() * cell.get_ddz_lamb2();

		divB += m_data.get_B(0, cell.get_global_node_id(3)).get_x() * cell.get_ddx_lamb3();
		divB += m_data.get_B(0, cell.get_global_node_id(3)).get_y() * cell.get_ddy_lamb3();
		divB += m_data.get_B(0, cell.get_global_node_id(3)).get_z() * cell.get_ddz_lamb3();

		if (divE != m_data.get_rho(0, c)/m_eps || divB != 0.)
		{
			return false;
		}
	}

	return true;
}

void UnsteadyMaxwellKernel::build_matrix_A()
{
	double coef_E(m_theta*m_dt/(m_mu*m_eps));
	double coef_B(m_theta*m_dt);

	std::vector<size_t> cells_id;
	size_t nb_cells;

	Cell cell(m_mesh.get_cell(0));
	size_t S0, S1, S2, S3;
	
	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		if (boundary_node(n))
		{
			mat_A.set_boundary_equation(6*n);
			mat_A.set_boundary_equation(6*n+1);
			mat_A.set_boundary_equation(6*n+2);
			mat_A.set_boundary_equation(6*n+3);
			mat_A.set_boundary_equation(6*n+4);
			mat_A.set_boundary_equation(6*n+5);
		}

		else
		{
			cells_id = m_mesh.get_neighbor_cells_id(n);
			nb_cells = cells_id.size();

			for (size_t c = 0; c < nb_cells; c++)
			{
				cell = m_mesh.get_cell(cells_id[c]);
				
				S0 = cell.get_global_node_id(0);
				S1 = cell.get_global_node_id(1);
				S2 = cell.get_global_node_id(2);
				S3 = cell.get_global_node_id(3);
				
				//finite difference term
				for (size_t i = 0; i < 6; i++) //common part to all 6 equations
				{
					mat_A.add_to_coef(6*n + i, 6*S0 + i, 0.25);
					mat_A.add_to_coef(6*n + i, 6*S1 + i, 0.25);
					mat_A.add_to_coef(6*n + i, 6*S2 + i, 0.25);
					mat_A.add_to_coef(6*n + i, 6*S3 + i, 0.25);
				}

				//rotB term
				//Ex equation
				mat_A.add_to_coef(6*n, 6*S0+5, -coef_E*cell.get_ddy_lamb0());
				mat_A.add_to_coef(6*n, 6*S1+5, -coef_E*cell.get_ddy_lamb1());
				mat_A.add_to_coef(6*n, 6*S2+5, -coef_E*cell.get_ddy_lamb2());
				mat_A.add_to_coef(6*n, 6*S3+5, -coef_E*cell.get_ddy_lamb3());

				mat_A.add_to_coef(6*n, 6*S0+4, coef_E*cell.get_ddz_lamb0());
				mat_A.add_to_coef(6*n, 6*S1+4, coef_E*cell.get_ddz_lamb1());
				mat_A.add_to_coef(6*n, 6*S2+4, coef_E*cell.get_ddz_lamb2());
				mat_A.add_to_coef(6*n, 6*S3+4, coef_E*cell.get_ddz_lamb3());

				//Ey equation
				mat_A.add_to_coef(6*n+1, 6*S0+3, -coef_E*cell.get_ddz_lamb0());
				mat_A.add_to_coef(6*n+1, 6*S1+3, -coef_E*cell.get_ddz_lamb1());
				mat_A.add_to_coef(6*n+1, 6*S2+3, -coef_E*cell.get_ddz_lamb2());
				mat_A.add_to_coef(6*n+1, 6*S3+3, -coef_E*cell.get_ddz_lamb3());

				mat_A.add_to_coef(6*n+1, 6*S0+5, coef_E*cell.get_ddx_lamb0());
				mat_A.add_to_coef(6*n+1, 6*S1+5, coef_E*cell.get_ddx_lamb1());
				mat_A.add_to_coef(6*n+1, 6*S2+5, coef_E*cell.get_ddx_lamb2());
				mat_A.add_to_coef(6*n+1, 6*S3+5, coef_E*cell.get_ddx_lamb3());

				//Ez equation
				mat_A.add_to_coef(6*n+2, 6*S0+4, -coef_E*cell.get_ddx_lamb0());
				mat_A.add_to_coef(6*n+2, 6*S1+4, -coef_E*cell.get_ddx_lamb1());
				mat_A.add_to_coef(6*n+2, 6*S2+4, -coef_E*cell.get_ddx_lamb2());
				mat_A.add_to_coef(6*n+2, 6*S3+4, -coef_E*cell.get_ddx_lamb3());

				mat_A.add_to_coef(6*n+2, 6*S0+3, coef_E*cell.get_ddy_lamb0());
				mat_A.add_to_coef(6*n+2, 6*S1+3, coef_E*cell.get_ddy_lamb1());
				mat_A.add_to_coef(6*n+2, 6*S2+3, coef_E*cell.get_ddy_lamb2());
				mat_A.add_to_coef(6*n+2, 6*S3+3, coef_E*cell.get_ddy_lamb3());


				//Bx equation
				mat_A.add_to_coef(6*n+3, 6*S0+2, coef_B*cell.get_ddy_lamb0());
				mat_A.add_to_coef(6*n+3, 6*S1+2, coef_B*cell.get_ddy_lamb1());
				mat_A.add_to_coef(6*n+3, 6*S2+2, coef_B*cell.get_ddy_lamb2());
				mat_A.add_to_coef(6*n+3, 6*S3+2, coef_B*cell.get_ddy_lamb3());

				mat_A.add_to_coef(6*n+3, 6*S0+1, -coef_B*cell.get_ddz_lamb0());
				mat_A.add_to_coef(6*n+3, 6*S1+1, -coef_B*cell.get_ddz_lamb1());
				mat_A.add_to_coef(6*n+3, 6*S2+1, -coef_B*cell.get_ddz_lamb2());
				mat_A.add_to_coef(6*n+3, 6*S3+1, -coef_B*cell.get_ddz_lamb3());

				//By equation
				mat_A.add_to_coef(6*n+4, 6*S0, coef_B*cell.get_ddz_lamb0());
				mat_A.add_to_coef(6*n+4, 6*S1, coef_B*cell.get_ddz_lamb1());
				mat_A.add_to_coef(6*n+4, 6*S2, coef_B*cell.get_ddz_lamb2());
				mat_A.add_to_coef(6*n+4, 6*S3, coef_B*cell.get_ddz_lamb3());

				mat_A.add_to_coef(6*n+4, 6*S0+2, -coef_B*cell.get_ddx_lamb0());
				mat_A.add_to_coef(6*n+4, 6*S1+2, -coef_B*cell.get_ddx_lamb1());
				mat_A.add_to_coef(6*n+4, 6*S2+2, -coef_B*cell.get_ddx_lamb2());
				mat_A.add_to_coef(6*n+4, 6*S3+2, -coef_B*cell.get_ddx_lamb3());

				//Bz equation
				mat_A.add_to_coef(6*n+5, 6*S0+1, coef_B*cell.get_ddx_lamb0());
				mat_A.add_to_coef(6*n+5, 6*S1+1, coef_B*cell.get_ddx_lamb1());
				mat_A.add_to_coef(6*n+5, 6*S2+1, coef_B*cell.get_ddx_lamb2());
				mat_A.add_to_coef(6*n+5, 6*S3+1, coef_B*cell.get_ddx_lamb3());

				mat_A.add_to_coef(6*n+5, 6*S0, -coef_B*cell.get_ddy_lamb0());
				mat_A.add_to_coef(6*n+5, 6*S1, -coef_B*cell.get_ddy_lamb1());
				mat_A.add_to_coef(6*n+5, 6*S2, -coef_B*cell.get_ddy_lamb2());
				mat_A.add_to_coef(6*n+5, 6*S3, -coef_B*cell.get_ddy_lamb3());
			}
		}
	}
}

void UnsteadyMaxwellKernel::build_vector_B(size_t t)
{
	double coef_E((1.-m_theta)*m_dt/(m_mu*m_eps));
	double coef_B((1.-m_theta)*m_dt);

	std::vector<size_t> cells_id;
	size_t nb_cells;

	const Cell* ptr_cell;
	size_t S0, S1, S2, S3;

	Vec3D E0, E1, E2, E3;
	Vec3D B0, B1, B2, B3;
	
	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		if (boundary_node(n))
		{
			E0 = m_data.get_boundary_condition_E(n);
			B0 = m_data.get_boundary_condition_B(n);

			vec_B[6*n] = E0.get_x();
			vec_B[6*n+1] = E0.get_y();
			vec_B[6*n+2] = E0.get_z();
			vec_B[6*n+3] = B0.get_x();
			vec_B[6*n+4] = B0.get_y();
			vec_B[6*n+5] = B0.get_z();
		}

		else
		{		
			cells_id = m_mesh.get_neighbor_cells_id(n);
			nb_cells = cells_id.size();

			for (size_t c = 0; c < nb_cells; c++)
			{
				ptr_cell = m_mesh.get_cell_ptr(cells_id[c]);

				S0 = ptr_cell->get_global_node_id(0);
				S1 = ptr_cell->get_global_node_id(1);
				S2 = ptr_cell->get_global_node_id(2);
				S3 = ptr_cell->get_global_node_id(3);

				E0 = m_data.get_E(t, S0);
				E1 = m_data.get_E(t, S1);
				E2 = m_data.get_E(t, S2);
				E3 = m_data.get_E(t, S3);

				B0 = m_data.get_B(t, S0);
				B1 = m_data.get_B(t, S1);
				B2 = m_data.get_B(t, S2);
				B3 = m_data.get_B(t, S3);


				//fintite difference term
				vec_B[6*n] += 0.25*E0.get_x();
				vec_B[6*n] += 0.25*E1.get_x();
				vec_B[6*n] += 0.25*E2.get_x();
				vec_B[6*n] += 0.25*E3.get_x();

				vec_B[6*n+1] += 0.25*E0.get_y();
				vec_B[6*n+1] += 0.25*E1.get_y();
				vec_B[6*n+1] += 0.25*E2.get_y();
				vec_B[6*n+1] += 0.25*E3.get_y();

				vec_B[6*n+2] += 0.25*E0.get_z();
				vec_B[6*n+2] += 0.25*E1.get_z();
				vec_B[6*n+2] += 0.25*E2.get_z();
				vec_B[6*n+2] += 0.25*E3.get_z();

				vec_B[6*n+3] += 0.25*B0.get_x();
				vec_B[6*n+3] += 0.25*B1.get_x();
				vec_B[6*n+3] += 0.25*B2.get_x();
				vec_B[6*n+3] += 0.25*B3.get_x();

				vec_B[6*n+4] += 0.25*B0.get_y();
				vec_B[6*n+4] += 0.25*B1.get_y();
				vec_B[6*n+4] += 0.25*B2.get_y();
				vec_B[6*n+4] += 0.25*B3.get_y();

				vec_B[6*n+5] += 0.25*B0.get_z();
				vec_B[6*n+5] += 0.25*B1.get_z();
				vec_B[6*n+5] += 0.25*B2.get_z();
				vec_B[6*n+5] += 0.25*B3.get_z();


				//rotB term
				//Ex equation
				vec_B[6*n] += coef_E*ptr_cell->get_ddy_lamb0()*B0.get_z();
				vec_B[6*n] += coef_E*ptr_cell->get_ddy_lamb1()*B1.get_z();
				vec_B[6*n] += coef_E*ptr_cell->get_ddy_lamb2()*B2.get_z();
				vec_B[6*n] += coef_E*ptr_cell->get_ddy_lamb3()*B3.get_z();

				vec_B[6*n] -= coef_E*ptr_cell->get_ddz_lamb0()*B0.get_y();
				vec_B[6*n] -= coef_E*ptr_cell->get_ddz_lamb1()*B1.get_y();
				vec_B[6*n] -= coef_E*ptr_cell->get_ddz_lamb2()*B2.get_y();
				vec_B[6*n] -= coef_E*ptr_cell->get_ddz_lamb3()*B3.get_y();

				//Ey equation
				vec_B[6*n+1] += coef_E*ptr_cell->get_ddz_lamb0()*B0.get_x();
				vec_B[6*n+1] += coef_E*ptr_cell->get_ddz_lamb1()*B1.get_x();
				vec_B[6*n+1] += coef_E*ptr_cell->get_ddz_lamb2()*B2.get_x();
				vec_B[6*n+1] += coef_E*ptr_cell->get_ddz_lamb3()*B3.get_x();

				vec_B[6*n+1] -= coef_E*ptr_cell->get_ddx_lamb0()*B0.get_z();
				vec_B[6*n+1] -= coef_E*ptr_cell->get_ddx_lamb1()*B1.get_z();
				vec_B[6*n+1] -= coef_E*ptr_cell->get_ddx_lamb2()*B2.get_z();
				vec_B[6*n+1] -= coef_E*ptr_cell->get_ddx_lamb3()*B3.get_z();

				//Ez equation
				vec_B[6*n+2] += coef_E*ptr_cell->get_ddx_lamb0()*B0.get_y();
				vec_B[6*n+2] += coef_E*ptr_cell->get_ddx_lamb1()*B1.get_y();
				vec_B[6*n+2] += coef_E*ptr_cell->get_ddx_lamb2()*B2.get_y();
				vec_B[6*n+2] += coef_E*ptr_cell->get_ddx_lamb3()*B3.get_y();

				vec_B[6*n+2] -= coef_E*ptr_cell->get_ddy_lamb0()*B0.get_x();
				vec_B[6*n+2] -= coef_E*ptr_cell->get_ddy_lamb1()*B1.get_x();
				vec_B[6*n+2] -= coef_E*ptr_cell->get_ddy_lamb2()*B2.get_x();
				vec_B[6*n+2] -= coef_E*ptr_cell->get_ddy_lamb3()*B3.get_x();

				//Bx equation
				vec_B[6*n+3] -= coef_B*ptr_cell->get_ddy_lamb0()*E0.get_z();
				vec_B[6*n+3] -= coef_B*ptr_cell->get_ddy_lamb1()*E1.get_z();
				vec_B[6*n+3] -= coef_B*ptr_cell->get_ddy_lamb2()*E2.get_z();
				vec_B[6*n+3] -= coef_B*ptr_cell->get_ddy_lamb3()*E3.get_z();

				vec_B[6*n+3] += coef_B*ptr_cell->get_ddz_lamb0()*E0.get_y();
				vec_B[6*n+3] += coef_B*ptr_cell->get_ddz_lamb1()*E1.get_y();
				vec_B[6*n+3] += coef_B*ptr_cell->get_ddz_lamb2()*E2.get_y();
				vec_B[6*n+3] += coef_B*ptr_cell->get_ddz_lamb3()*E3.get_y();

				//By equation
				vec_B[6*n+4] -= coef_B*ptr_cell->get_ddz_lamb0()*E0.get_x();
				vec_B[6*n+4] -= coef_B*ptr_cell->get_ddz_lamb1()*E1.get_x();
				vec_B[6*n+4] -= coef_B*ptr_cell->get_ddz_lamb2()*E2.get_x();
				vec_B[6*n+4] -= coef_B*ptr_cell->get_ddz_lamb3()*E3.get_x();

				vec_B[6*n+4] += coef_B*ptr_cell->get_ddx_lamb0()*E0.get_z();
				vec_B[6*n+4] += coef_B*ptr_cell->get_ddx_lamb1()*E1.get_z();
				vec_B[6*n+4] += coef_B*ptr_cell->get_ddx_lamb2()*E2.get_z();
				vec_B[6*n+4] += coef_B*ptr_cell->get_ddx_lamb3()*E3.get_z();

				//Bz equation
				vec_B[6*n+5] -= coef_B*ptr_cell->get_ddx_lamb0()*E0.get_y();
				vec_B[6*n+5] -= coef_B*ptr_cell->get_ddx_lamb1()*E1.get_y();
				vec_B[6*n+5] -= coef_B*ptr_cell->get_ddx_lamb2()*E2.get_y();
				vec_B[6*n+5] -= coef_B*ptr_cell->get_ddx_lamb3()*E3.get_y();

				vec_B[6*n+5] += coef_B*ptr_cell->get_ddy_lamb0()*E0.get_x();
				vec_B[6*n+5] += coef_B*ptr_cell->get_ddy_lamb1()*E1.get_x();
				vec_B[6*n+5] += coef_B*ptr_cell->get_ddy_lamb2()*E2.get_x();
				vec_B[6*n+5] += coef_B*ptr_cell->get_ddy_lamb3()*E3.get_x();
			}
		}
	}
}

double UnsteadyMaxwellKernel::GS_step(size_t i) const
{
	double new_Xi(vec_B[i]);
	std::vector<std::pair<size_t, double>> line(mat_A.get_reduced_line_wo_diag(i));

	for (size_t j = 0; j < line.size(); j++)
	{
		new_Xi -= line[j].second * vec_X[line[j].first];
	}

	return new_Xi/mat_A.get_coef(i, i);
}

double UnsteadyMaxwellKernel::GS_iteration()
{
	double new_Xi;
	double energy_diff(0.);

	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		for (size_t i = 0; i < 3; i++) //updating E
		{
			new_Xi = GS_step(6*n + i);
			energy_diff += m_eps*(new_Xi*new_Xi - vec_X[6*n + i]*vec_X[6*n + i]);
			vec_X[6*n + i] = new_Xi;
		}

		for (size_t i = 3; i < 6; i++) //updating B
		{
			new_Xi = GS_step(6*n + i);
			energy_diff += (new_Xi*new_Xi - vec_X[6*n + i]*vec_X[6*n + i])/m_mu;
			vec_X[6*n + i] = new_Xi;
		}
	}

	return 0.5*energy_diff/m_nb_nodes;
}

void UnsteadyMaxwellKernel::display_progression(size_t step) const
{
	std::cout << "\r - Simulating (step " << step+1 << " / " << m_nb_steps << " -- " << 100*(step+1)/m_nb_steps << " % completed)     " << std::flush;
}

void UnsteadyMaxwellKernel::simulate()
{
	if (!check_divergence_in_ID())
	{
		throw std::invalid_argument("Initial data is not solution of Maxwell-Gauss and Maxwell-Thompson equations");
	}

	build_matrix_A();

	double mean_energy_diff;

	for (size_t t = 0; t < m_nb_steps; t++)
	{
		display_progression(t);

		build_vector_B(t);
		size_t k(0);

		do
		{
			mean_energy_diff = GS_iteration();
			k += 1;
		}
		while (mean_energy_diff > m_accuracy && k < m_max_nb_iterations);

		if (k == m_max_nb_iterations)
		{
			throw std::runtime_error("Maximun number of iterations reached while computing step " + t+1);
		}

		for (size_t n = 0; n < m_nb_nodes; n++)
		{
			m_data.set_E(t+1, n, Vec3D(vec_X[6*n], vec_X[6*n + 1], vec_X[6*n + 2]));
			m_data.set_B(t+1, n, Vec3D(vec_X[6*n + 3], vec_X[6*n + 4], vec_X[6*n + 5]));
		}
	}
}