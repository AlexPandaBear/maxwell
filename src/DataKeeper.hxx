#pragma once

#include <stdexcept>
#include <vector>
#include <list>
#include <iostream>
#include <string>
#include "Field.hxx"
#include "Vec3D.hxx"
#include "BoundaryCondition.hxx"

class DataKeeper
{
private:
	double m_eps, m_mu;

	size_t m_nb_steps;
	double m_t_max;
	double m_theta;
	double m_accuracy;
	size_t m_max_nb_iterations;

	size_t m_nb_nodes;
	size_t m_nb_cells;

	std::vector<Field<double>> m_rho;
	std::vector<Field<Vec3D>> m_j;
	std::vector<Field<Vec3D>> m_E, m_B;

	std::vector<BoundaryCondition> m_BC;

public:
	DataKeeper();
	DataKeeper(std::string file);
	~DataKeeper();

	double get_epsilon() const;
	double get_mu() const;
	double get_t_max() const;
	double get_theta() const;
	double get_accuracy() const;
	size_t get_max_nb_iterations() const;

	size_t get_nb_steps() const;
	size_t get_nb_nodes() const;
	size_t get_nb_cells() const;

	void reset_dimensions(size_t nb_steps, size_t nb_nodes, size_t nb_cells);

	void set_epsilon(double epsilon);
	void set_mu(double mu);
	void set_t_max(double t_max);
	void set_theta(double theta);
	void set_accuracy(double accuracy);
	void set_max_nb_iterations(size_t max_nb_iterations);

	double get_rho(size_t t, size_t node_nb) const;
	void set_rho(size_t t, size_t node_nb, double rho);

	Vec3D get_j(size_t t, size_t node_nb) const;
	void set_j(size_t t, size_t node_nb, Vec3D j);

	Vec3D get_E(size_t t, size_t node_nb) const;
	void set_E(size_t t, size_t node_nb, Vec3D E);
	
	Vec3D get_B(size_t t, size_t node_nb) const;
	void set_B(size_t t, size_t node_nb, Vec3D B);
	
	void erase_BCs();
	void add_BC(size_t node_nb, Vec3D E, Vec3D B);
	Vec3D get_boundary_condition_E(size_t node_nb) const;
	Vec3D get_boundary_condition_B(size_t node_nb) const;

	bool boundary_node(size_t node_nb) const;

	bool check_ID_BC_compatibility() const;

	void save(std::string file) const;
};