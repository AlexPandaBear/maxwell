#include "Mesh3D.hxx"

Mesh3D::Mesh3D() :
	m_nb_nodes(0),
	m_nb_cells(0),
	m_nodes_xyz(VectorField(0)),
	m_cells(0, Cell(m_nodes_xyz)) {}

//Mesh3D::Mesh3D(std::string file) {} //TODO

Mesh3D::~Mesh3D() {}

void Mesh3D::generate_grid_mesh(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, size_t nx, size_t ny, size_t nz)
{
	m_nb_nodes = nx*ny*nz;
	m_nb_cells = 5*(nx-1)*(ny-1)*(nz-1);

	double dx((x_max-x_min)/(nx-1));
	double dy((y_max-y_min)/(ny-1));
	double dz((z_max-z_min)/(nz-1));


	m_nodes_xyz = VectorField(m_nb_nodes);

	size_t id(0);

	for (size_t i = 0; i < nx; i++)
	{
		for (size_t j = 0; j < ny; j++)
		{
			for (size_t k = 0; k < nz; k++)
			{
				m_nodes_xyz.set_value(id, x_min+i*dx, y_min+j*dy, z_min+k*dz);
				id++;
			}
		}
	}


	std::cout << "hello" << std::endl;
	m_cells = std::vector<Cell>(m_nb_cells, Cell(m_nodes_xyz, 0, 0, 0, 0));
	id = 0;

	for (size_t i = 0; i < nx-1; i++)
	{
		for (size_t j = 0; j < ny-1; j++)
		{
			for (size_t k = 0; k < nz-1; k++)
			{
				m_cells[id].set_nodes_id(i*ny*nz + j*nz + k, i*ny*nz + j*nz + k+1, i*ny*nz + (j+1)*nz + k, (i+1)*ny*nz + j*nz + k);
				id++;

				m_cells[id].set_nodes_id(i*ny*nz + j*nz + k+1, i*ny*nz + (j+1)*nz + k+1, i*ny*nz + (j+1)*nz + k, (i+1)*ny*nz + (j+1)*nz + k+1);
				id++;

				m_cells[id].set_nodes_id(i*ny*nz + (j+1)*nz + k, (i+1)*ny*nz + j*nz + k, (i+1)*ny*nz + (j+1)*nz + k, (i+1)*ny*nz + (j+1)*nz + k+1);
				id++;

				m_cells[id].set_nodes_id((i+1)*ny*nz + j*nz + k, i*ny*nz + j*nz + k+1, (i+1)*ny*nz + j*nz + (k+1), (i+1)*ny*nz + (j+1)*nz + k+1);
				id++;

				m_cells[id].set_nodes_id((i+1)*ny*nz + j*nz + k, i*ny*nz + (j+1)*nz + k, i*ny*nz + j*nz + (k+1), (i+1)*ny*nz + (j+1)*nz + k+1);
				id++;
			}
		}
	}
}

size_t Mesh3D::get_nb_nodes() const
{
	return m_nb_nodes;
}

size_t Mesh3D::get_nb_cells() const
{
	return m_nb_cells;
}

Vec3D Mesh3D::get_node_xyz(size_t node_id)
{
	return m_nodes_xyz.get_value(node_id);
}

Cell const& Mesh3D::get_cell(size_t cell_id) const
{
	return m_cells[cell_id];
}

const Cell* Mesh3D::get_cell_ptr(size_t cell_id) const
{
	return &m_cells[cell_id];
}

std::vector<size_t> Mesh3D::get_neighbor_cells_id(size_t node_id) const
{
	std::vector<size_t> neighbors;

	for (size_t i = 0; i < m_nb_cells; i++)
	{
		if (m_cells[i].contains(node_id))
		{
			neighbors.push_back(i);
		}
	}

	return neighbors;
}

//void Mesh3D::save(std::string file) const {} //TODO

VectorField const& Mesh3D::get_all_nodes_xyz() const
{
	return m_nodes_xyz;
}