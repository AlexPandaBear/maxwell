#pragma once

#include "UnsteadyMaxwellKernel.hxx"

class SimManager
{
private:
	Mesh3D m_mesh;
	DataKeeper m_data;

	double generate_random_double(double min, double max) const;
	double compute_distance_to_axis(Vec3D pt, Vec3D axis_pt_A, Vec3D axis_pt_B) const;
	double estimate_wire_intersection(Vec3D pt_A, Vec3D pt_B, double radius, size_t cell_nb, size_t sample_size) const;

public:
	SimManager();
	SimManager(std::string simulation_name);
	~SimManager();

	void set_constants(double epsilon, double mu);
	void set_simulation_parameters(double t_max, size_t nb_steps, double theta, double accuracy, size_t max_nb_iterations, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, size_t nx, size_t ny, size_t nz);
	void define_initial_state_background_values(double rho0, Vec3D j0, Vec3D E0, Vec3D B0);
	void add_wire(std::vector<Vec3D> wire_skeleton, double wire_radius, double wire_current);
	void add_boundary_condition(size_t node_nb, Vec3D E, Vec3D B);

	void simulate();

	void save(std::string simulation_name) const;
};