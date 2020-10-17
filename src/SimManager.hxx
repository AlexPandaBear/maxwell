#pragma once

#include "UnsteadyMaxwellKernel.hxx"

class SimManager
{
private:
	Mesh3D m_mesh;
	DataKeeper m_data;

public:
	SimManager();
	SimManager(std::string simulation_name);
	~SimManager();

	void set_constants(double epsilon, double mu);
	void set_simulation_parameters(double t_max, size_t nb_steps, double theta, double accuracy, size_t max_nb_iterations, size_t nb_nodes, size_t nb_cells);
	void define_initial_state(std::vector<double> rho0, std::vector<Vec3D> j0, std::vector<Vec3D> E0, std::vector<Vec3D> B0);
	void add_boundary_condition(size_t node_nb, Vec3D E, Vec3D B);

	void simulate();

	void save(std::string simulation_name) const;	
};