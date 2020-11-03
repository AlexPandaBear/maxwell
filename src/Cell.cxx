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
	Vec3D N0(m_nodes_ref.get_value(m_N0_id)), N1(m_nodes_ref.get_value(m_N1_id)), N2(m_nodes_ref.get_value(m_N2_id)), N3(m_nodes_ref.get_value(m_N3_id));

	double x0(N0.get_x()), y0(N0.get_y()), z0(N0.get_z());
	double x1(N1.get_x()), y1(N1.get_y()), z1(N1.get_z());
	double x2(N2.get_x()), y2(N2.get_y()), z2(N2.get_z());
	double x3(N3.get_x()), y3(N3.get_y()), z3(N3.get_z());


	m_volume = ( (x1-x0)*(y2-y0)*(z3-z0) + (x3-x0)*(y1-y0)*(z2-z0) + (x2-x0)*(y3-y0)*(z1-z0) - (x1-x0)*(y3-y0)*(z2-z0) - (x2-x0)*(y1-y0)*(z3-z0) - (x3-x0)*(y2-y0)*(z1-z0) )/6.;


	double coef(1./(6.*m_volume));

	m_ddx_lamb0 = coef * ( (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1) );
	m_ddy_lamb0 = coef * ( (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1) );
	m_ddz_lamb0 = coef * ( (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1) );

	if (m_ddx_lamb0*(x0-x1) + m_ddy_lamb0*(y0-y1) + m_ddz_lamb0*(z0-z1) < 0.)
	{
		m_ddx_lamb0 = - m_ddx_lamb0;
		m_ddy_lamb0 = - m_ddy_lamb0;
		m_ddz_lamb0 = - m_ddz_lamb0;
	}


	m_ddx_lamb1 = coef * ( (y3-y2)*(z0-z2) - (z3-z2)*(y0-y2) );
	m_ddy_lamb1 = coef * ( (z3-z2)*(x0-x2) - (x3-x2)*(z0-z2) );
	m_ddz_lamb1 = coef * ( (x3-x2)*(y0-y2) - (y3-y2)*(x0-x2) );

	if (m_ddx_lamb1*(x1-x2) + m_ddy_lamb1*(y1-y2) + m_ddz_lamb1*(z1-z2) < 0.)
	{
		m_ddx_lamb1 = - m_ddx_lamb1;
		m_ddy_lamb1 = - m_ddy_lamb1;
		m_ddz_lamb1 = - m_ddz_lamb1;
	}


	m_ddx_lamb2 = coef * ( (y0-y3)*(z1-z3) - (z0-z3)*(y1-y3) );
	m_ddy_lamb2 = coef * ( (z0-z3)*(x1-x3) - (x0-x3)*(z1-z3) );
	m_ddz_lamb2 = coef * ( (x0-x3)*(y1-y3) - (y0-y3)*(x1-x3) );

	if (m_ddx_lamb2*(x2-x3) + m_ddy_lamb2*(y2-y3) + m_ddz_lamb2*(z2-z3) < 0.)
	{
		m_ddx_lamb2 = - m_ddx_lamb2;
		m_ddy_lamb2 = - m_ddy_lamb2;
		m_ddz_lamb2 = - m_ddz_lamb2;
	}


	m_ddx_lamb3 = coef * ( (y1-y0)*(z2-z0) - (z1-z0)*(y2-y0) );
	m_ddy_lamb3 = coef * ( (z1-z0)*(x2-x0) - (x1-x0)*(z2-z0) );
	m_ddz_lamb3 = coef * ( (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0) );

	if (m_ddx_lamb3*(x3-x0) + m_ddy_lamb3*(y3-y0) + m_ddz_lamb3*(z3-z0) < 0.)
	{
		m_ddx_lamb3 = - m_ddx_lamb3;
		m_ddy_lamb3 = - m_ddy_lamb3;
		m_ddz_lamb3 = - m_ddz_lamb3;
	}
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

double Cell::get_ddx_lamb0() const
{
	return m_ddx_lamb0;
}

double Cell::get_ddx_lamb1() const
{
	return m_ddx_lamb1;
}

double Cell::get_ddx_lamb2() const
{
	return m_ddx_lamb2;
}

double Cell::get_ddx_lamb3() const
{
	return m_ddx_lamb3;
}

double Cell::get_ddy_lamb0() const
{
	return m_ddy_lamb0;
}

double Cell::get_ddy_lamb1() const
{
	return m_ddy_lamb1;
}

double Cell::get_ddy_lamb2() const
{
	return m_ddy_lamb2;
}

double Cell::get_ddy_lamb3() const
{
	return m_ddy_lamb3;
}

double Cell::get_ddz_lamb0() const
{
	return m_ddz_lamb0;
}

double Cell::get_ddz_lamb1() const
{
	return m_ddz_lamb1;
}

double Cell::get_ddz_lamb2() const
{
	return m_ddz_lamb2;
}

double Cell::get_ddz_lamb3() const
{
	return m_ddz_lamb3;
}
/*
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