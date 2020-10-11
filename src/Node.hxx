#pragma once

class Node
{
private:
	double m_x;
	double m_y;
	double m_z;

	bool boundary;

public:
	Node();
	Node(double x, double y, double z, bool boundary);
	~Node();

	double get_x() const;
	double get_y() const;
	double get_z() const;
	bool get_boudary() const;

	void set_x(double x);
	void set_y(double y);
	void set_z(double z);
	void set_boundary(bool boundary);
};