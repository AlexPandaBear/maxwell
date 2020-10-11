#include "UnsteadyMaxwellKernel.hxx"

UnsteadyMaxwellKernel::UnsteadyMaxwellKernel(Mesh& mesh, DataKeeper& data) :
	m_eps(data.get_epsilon()),
	m_mu(data.get_mu()),
	m_nb_steps(data.get_nb_steps()),
	m_t_max(data.get_t_max()),
	m_theta(data.get_theta()),
	m_nb_nodes(data.get_nb_nodes()),
	m_nb_cells(data.get_nb_triangles()),
	m_mesh(mesh),
	m_data(data),
	mat_A(6*m_nb_nodes),
	vec_X(6*m_nb_nodes, 0.),
	vec_B(6*m_nb_nodes, 0.) {}

UnsteadyMaxwellKernel::~UnsteadyMaxwellKernel() {}

void UnsteadyMaxwellKernel::check_ID_BC_compatibility() const;

void UnsteadyMaxwellKernel::check_divergence_in_ID() const;

void UnsteadyMaxwellKernel::build_matrix_A()
{
	double coef_E(m_theta*m_dt/(m_mu*m_eps));
	double coef_B(m_theta*m_dt);

	std::vector<size_t> cells_id;
	size_t nb_cells;

	Cell cell;
	size_t S0, S1, S2, S3;
	
	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		cells_id = m_mesh.get_neighbor_cells_id(n);
		nb_cells = cells_id.size();

		for (size_t i = 0; i < 6; i++)
		{
			mat_A.set_coef(6*n + i, 6*n + i, nb_cells); // diagonal term for all 6 equations
		}

		for (size_t c = 0; c < nb_cells; c++)
		{
			cell = m_mesh.get_cell(cells_id[c]);
			S0 = cell.get_node_id(0);
			S1 = cell.get_node_id(1);
			S2 = cell.get_node_id(2);
			S3 = cell.get_node_id(3);

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

void UnsteadyMaxwellKernel::simulate()
{
	if (!check_divergence_in_ID() || !check_ID_BC_compatibility())
	{
		throw ...;
	}

	build_matrix_A();

	for (size_t t = 0; t < m_nb_steps; t++)
	{
		//compute values of E, B, rotE, rotB for each cell at n
		//build AX=B for E, B at each cell at n+1
		//solve for E,B at each node w/ BC
	}
}