#pragma once

#include "Mesh.hxx"
#include "DataKeeper.hxx"

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
	void set_simulation_parameters(double t_max, size_t nb_steps, double theta, size_t nb_nodes, size_t nb_cells);
	void define_initial_state(std::vector<double> Ex_0, std::vector<double> Ey_0, std::vector<double> Ez_0, std::vector<double> Bx_0, std::vector<double> By_0, std::vector<double> Bz_0);
	void add_boundary_condition(size_t node_nb, Vec3D E, Vec3D B);

	void simulate();

	void save(std::string simulation_name);	
};