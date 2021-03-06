#pragma once

#include "UnsteadyMaxwellKernel.hxx"
#include "DataProcessor.hxx"

class SimManager
{
private:
	Mesh3D m_mesh;
	DataKeeper m_data;
	DataProcessor m_processor;

	double generate_random_double(double min, double max) const;
	double compute_distance_to_axis(Vec3D pt, Vec3D axis_pt_A, Vec3D axis_pt_B) const;
	double estimate_wire_intersection(Vec3D pt_A, Vec3D pt_B, double radius, size_t cell_nb, size_t sample_size) const;

public:
	SimManager();
	//SimManager(std::string simulation_name);
	~SimManager();

	void set_constants(double epsilon, double mu, double sigma);
	void set_simulation_parameters(double t_max, size_t nb_steps, double theta, double accuracy, size_t max_nb_iterations, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, size_t nx, size_t ny, size_t nz);
	
	void define_initial_state(VectorField const& E0, VectorField const& B0);
	void define_initial_state_background_values(double rho0, Vec3D j0, Vec3D E0, Vec3D B0);
	void add_wire(std::vector<Vec3D> wire_skeleton, double wire_radius, double wire_current);
	
	void add_boundary_condition(size_t node_nb);
	void lock_all_boundary_nodes();

	void simulate();

	//void save(std::string simulation_name) const;

	size_t get_nb_nodes() const;
	size_t get_nb_cells() const;
	std::vector<size_t> get_node_ids(size_t cell_id) const;

	Vec3D get_node_xyz(size_t node_id) const;

	VectorField const& get_mesh();
	std::vector<double> const& get_time();

	ScalarField const& get_rho(size_t step);
	VectorField const& get_j(size_t step);
	VectorField const& get_E(size_t step);
	VectorField const& get_B(size_t step);
	
	ScalarField const& get_energy_density(size_t step);
	VectorField const& get_poynting_vector(size_t step);
	ScalarField const& get_poynting_vector_norm(size_t step);
};