#pragma once

#include "ElectrostaticKernel.hxx"

class ElectrostaticSimManager
{
private:
	const bool m_verbose;

	Mesh3D m_mesh;
	ElectrostaticDataKeeper m_data;

public:
	ElectrostaticSimManager(bool verbose = true);
	~ElectrostaticSimManager();

	void set_epsilon0(double epsilon0);
	void set_accuracy(double accuracy);
	void set_max_nb_iterations(size_t max_nb_iterarions);

	void generate_cube_mesh(double x_min, double x_max, size_t nx,
							double y_min, double y_max, size_t ny,
							double z_min, double z_max, size_t nz);

	void set_rho(size_t node_id, double rho);

	void set_rho_field(ScalarField const& rho);
	void set_epsilon_r_field(ScalarField const& epsilon_r);

	void simulate();

	size_t get_nb_nodes() const;
	Vec3D get_node_xyz(size_t node_id) const;

	VectorField const& get_mesh() const;
	ScalarField const& get_rho() const;
	ScalarField const& get_epsilon() const;
	ScalarField const& get_phi() const;
	VectorField const& get_E() const;	
};