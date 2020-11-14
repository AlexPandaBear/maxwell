#include "Cell.hxx"

Cell::Cell(VectorField const& nodes) :
	m_nodes_ref(nodes) {}

Cell::Cell(VectorField const& nodes, size_t N0_id, size_t N1_id, size_t N2_id, size_t N3_id) :
	m_nodes_ref(nodes),
	m_N0_id(N0_id),
	m_N1_id(N1_id),
	m_N2_id(N2_id),
	m_N3_id(N3_id)
{
	compute_nodes();
}

Cell::Cell(Cell const& cell) :
	m_nodes_ref(cell.get_nodes_ref()),
	m_N0_id(cell.get_global_node_id(0)),
	m_N1_id(cell.get_global_node_id(1)),
	m_N2_id(cell.get_global_node_id(2)),
	m_N3_id(cell.get_global_node_id(3))
{
	compute_nodes();
}

Cell::~Cell() {}

void Cell::compute_nodes()
{
	Vec3D N0(m_nodes_ref.get_value(m_N0_id));
	Vec3D N1(m_nodes_ref.get_value(m_N1_id));
	Vec3D N2(m_nodes_ref.get_value(m_N2_id));
	Vec3D N3(m_nodes_ref.get_value(m_N3_id));


	m_S_012 = Vec3D::cross_product(N1-N0, N2-N0) * 0.5;
	m_S_013 = Vec3D::cross_product(N1-N0, N3-N0) * 0.5;
	m_S_023 = Vec3D::cross_product(N2-N0, N3-N0) * 0.5;
	m_S_123 = Vec3D::cross_product(N2-N1, N3-N1) * 0.5;


	if (Vec3D::dot_product(m_S_012, N3-N0) > 0.)
	{
		m_S_012 *= -1.;
	}

	if (Vec3D::dot_product(m_S_013, N2-N0) > 0.)
	{
		m_S_013 *= -1.;
	}

	if (Vec3D::dot_product(m_S_023, N1-N0) > 0.)
	{
		m_S_023 *= -1.;
	}

	if (Vec3D::dot_product(m_S_123, N0-N1) > 0.)
	{
		m_S_123 *= -1.;
	}


	m_volume = Vec3D::dot_product(m_S_012, N0-N3)/3.;
}

void Cell::set_global_node_id(size_t local_node_id, size_t global_node_id)
{
	switch (local_node_id)
	{
		case 0: m_N0_id = global_node_id;
		case 1: m_N1_id = global_node_id;
		case 2: m_N2_id = global_node_id;
		case 3: m_N3_id = global_node_id;
		default: throw std::invalid_argument("Invalid local node number");
	}

	compute_nodes();
}

void Cell::set_nodes_id(size_t N0_id, size_t N1_id, size_t N2_id, size_t N3_id)
{
	m_N0_id = N0_id;
	m_N1_id = N1_id;
	m_N2_id = N2_id;
	m_N3_id = N3_id;

	compute_nodes();
}

size_t Cell::get_global_node_id(size_t local_node_id) const
{
	switch (local_node_id)
	{
		case 0: return m_N0_id;
		case 1: return m_N1_id;
		case 2: return m_N2_id;
		case 3: return m_N3_id;
		default: throw std::invalid_argument("Invalid local node number");
	}
}

Vec3D Cell::get_node_xyz(size_t local_node_id) const
{
	switch (local_node_id)
	{
		case 0: return m_nodes_ref.get_value(m_N0_id);
		case 1: return m_nodes_ref.get_value(m_N1_id);
		case 2: return m_nodes_ref.get_value(m_N2_id);
		case 3: return m_nodes_ref.get_value(m_N3_id);
		default: throw std::invalid_argument("Invalid local node number");
	}
}
/*
void Cell::set_nodes_ref(VectorField const& nodes)
{
	m_nodes_ref = nodes;
}
*/
VectorField const& Cell::get_nodes_ref() const
{
	return m_nodes_ref;
}

bool Cell::contains(size_t global_node_id) const
{
	if (global_node_id == m_N0_id || global_node_id == m_N1_id || global_node_id == m_N2_id || global_node_id == m_N3_id)
	{
		return true;
	}

	return false;
}

double Cell::get_volume() const
{
	return m_volume;
}

Vec3D Cell::get_surface(size_t opposite_node_local_id) const
{
	switch (opposite_node_local_id)
	{
		case 0: return m_S_123;
		case 1: return m_S_023;
		case 2: return m_S_013;
		case 3: return m_S_012;
		default: throw std::invalid_argument("Invalid local node number");
	}
}
/*
double Cell::get_ddx_lamb(size_t local_node_id) const
{
	switch (local_node_id)
	{
		case 0: return m_ddx_lamb0;
		case 1: return m_ddx_lamb1;
		case 2: return m_ddx_lamb2;
		case 3: return m_ddx_lamb3;
		default: throw std::invalid_argument("Invalid local node number");
	}
}

double Cell::get_ddy_lamb(size_t local_node_id) const
{
	switch (local_node_id)
	{
		case 0: return m_ddy_lamb0;
		case 1: return m_ddy_lamb1;
		case 2: return m_ddy_lamb2;
		case 3: return m_ddy_lamb3;
		default: throw std::invalid_argument("Invalid local node number");
	}
}

double Cell::get_ddz_lamb(size_t local_node_id) const
{
	switch (local_node_id)
	{
		case 0: return m_ddz_lamb0;
		case 1: return m_ddz_lamb1;
		case 2: return m_ddz_lamb2;
		case 3: return m_ddz_lamb3;
		default: throw std::invalid_argument("Invalid local node number");
	}
}

void Cell::display() const
{
	m_N0.display();
	m_N1.display();
	m_N2.display();
	m_N3.display();
}
*/
Cell Cell::operator=(Cell const& cell)
{
	Cell copy(cell);
	return copy;
}