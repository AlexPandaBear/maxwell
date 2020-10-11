#include "Mesh.hxx"

Mesh::Mesh() {}

Mesh::~Mesh() {}

void Mesh::generate_grid_mesh(double x_min, double x_max, double y_min, double y_max, size_t nx, size_t ny)
{
	m_nb_pts = nx*ny;
	m_nb_trgls = 2*(nx-1)*(ny-1);

	m_pts_x = std::vector<double>(m_nb_pts, 0.);
	m_pts_y = std::vector<double>(m_nb_pts, 0.);

	double dx((x_max-x_min)/(nx-1));
	double dy((y_max-y_min)/(ny-1));

	for (size_t i = 0; i < nx; i++)
	{
		for (size_t j = 0; j < ny; j++)
		{
			m_pts_x[i*ny + j] = x_min + i*dx;
			m_pts_y[i*ny + j] = y_min + j*dy;
		}
	}

	m_trgls = std::vector<Triangle>(m_nb_trgls, Triangle());

	for (size_t i = 0; i < nx-1; i++)
	{
		for (size_t j = 0; i < ny-1; j++)
		{
			m_trgls[i*2*(ny-1) + 2*j].set_Sid(i*ny + j, (i+1)*ny + j, (i+1)*ny + j+1);
			m_trgls[i*2*(ny-1) + 2*j].set_edges((j==0), (i+1==nx-1), false);

			m_trgls[i*2*(ny-1) + 2*j + 1].set_Sid(i*ny + j, i*ny + j+1, (i+1)*ny + j+1);
			m_trgls[i*2*(ny-1) + 2*j + 1].set_edges((i==0), (j+1==ny-1), false);
		}
	}

	size_t S1, S2, S3;

	for (size_t k = 0; k < m_nb_trgls; k++)
	{
		S1 = m_trgls[k].get_S1_id();
		S2 = m_trgls[k].get_S2_id();
		S3 = m_trgls[k].get_S3_id();

		m_trgls[k].set_coordinates(m_pts_x[S1], m_pts_x[S2], m_pts_x[S3], m_pts_y[S1], m_pts_y[S2], m_pts_y[S3]);
		m_trgls[k].compute();

	}
}