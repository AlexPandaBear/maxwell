#include "Mesh3D.hxx"

Mesh3D::Mesh3D() {}

Mesh3D::Mesh3D(std::string file) {} //TODO

Mesh3D::~Mesh3D() {}

void Mesh3D::generate_grid_mesh(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, size_t nx, size_t ny, size_t nz)
{
	m_nb_nodes = nx*ny*nz;
	m_nb_cells = 5*(nx-1)*(ny-1)*(nz-1);

	double dx((x_max-x_min)/(nx-1));
	double dy((y_max-y_min)/(ny-1));
	double dz((z_max-z_min)/(nz-1));

	m_nodes = std::vector<Node>();
	size_t id(0);

	Node n;
	for (size_t i = 0; i < nx; i++)
	{
		for (size_t j = 0; j < ny; j++)
		{
			for (size_t k = 0; k < nz; k++)
			{
				n = Node(id, x_min+i*dx, y_min+j*dy, z_min+k*dz);
				n.display();
				m_nodes.push_back(n);
				id++;
			}
		}
	}


	m_cells = std::vector<Cell>();

	for (size_t i = 0; i < nx-1; i++)
	{
		for (size_t j = 0; j < ny-1; j++)
		{
			for (size_t k = 0; k < nz-1; k++)
			{
				m_cells.push_back(Cell(m_nodes[i*ny*nz + j*nz + k], m_nodes[i*ny*nz + j*nz + k+1], m_nodes[i*ny*nz + (j+1)*nz + k], m_nodes[(i+1)*ny*nz + j*nz + k]));
				m_cells.push_back(Cell(m_nodes[i*ny*nz + j*nz + k+1], m_nodes[i*ny*nz + (j+1)*nz + k+1], m_nodes[i*ny*nz + (j+1)*nz + k], m_nodes[(i+1)*ny*nz + (j+1)*nz + k+1]));
				m_cells.push_back(Cell(m_nodes[i*ny*nz + (j+1)*nz + k], m_nodes[(i+1)*ny*nz + j*nz + k], m_nodes[(i+1)*ny*nz + (j+1)*nz + k], m_nodes[(i+1)*ny*nz + (j+1)*nz + k+1]));
				m_cells.push_back(Cell(m_nodes[(i+1)*ny*nz + j*nz + k], m_nodes[i*ny*nz + j*nz + k+1], m_nodes[(i+1)*ny*nz + j*nz + (k+1)], m_nodes[(i+1)*ny*nz + (j+1)*nz + k+1]));
				m_cells.push_back(Cell(m_nodes[(i+1)*ny*nz + j*nz + k], m_nodes[i*ny*nz + (j+1)*nz + k], m_nodes[i*ny*nz + j*nz + (k+1)], m_nodes[(i+1)*ny*nz + (j+1)*nz + k+1]));
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

Node Mesh3D::get_node(size_t node_nb) const
{
	return m_nodes[node_nb];
}

Cell const& Mesh3D::get_cell(size_t cell_nb) const
{
	return m_cells[cell_nb];
}

const Cell* Mesh3D::get_cell_ptr(size_t cell_nb) const
{
	return &m_cells[cell_nb];
}

std::vector<size_t> Mesh3D::get_neighbor_cells_id(size_t node_nb) const
{
	std::vector<size_t> neighbors;

	for (size_t i = 0; i < m_cells.size(); i++)
	{
		if (m_cells[i].contains(node_nb))
		{
			neighbors.push_back(i);
		}
	}

	return neighbors;
}

void Mesh3D::save(std::string file) const {} //TODO

Vec3D Mesh3D::get_node_xyz(size_t node_nb) const
{
	return get_node(node_nb).get_xyz();
}

Field<Vec3D> Mesh3D::get_all_nodes_xyz() const
{
	Field<Vec3D> f(m_nb_nodes);

	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		f.set_value(n, get_node(n).get_xyz());
	}

	return f;
}