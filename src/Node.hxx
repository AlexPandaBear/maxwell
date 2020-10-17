#pragma once

#include <cstddef>

class Node
{
private:
	size_t m_id;

	double m_x;
	double m_y;
	double m_z;

	bool m_boundary;

public:
	Node();
	Node(size_t id, double x, double y, double z, bool boundary = false);
	~Node();

	size_t get_id() const;
	double get_x() const;
	double get_y() const;
	double get_z() const;
	bool get_boudary() const;

	void set_x(double x);
	void set_y(double y);
	void set_z(double z);
	void set_boundary(bool boundary);
};