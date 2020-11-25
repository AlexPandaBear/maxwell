#pragma once

#include <stdexcept>
#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <memory>
#include "ScalarField.hxx"
#include "VectorField.hxx"
#include "Vec3D.hxx"

class DataKeeper
{
private:
	double m_eps, m_mu, m_sigma;

	size_t m_nb_steps;
	double m_t_max;
	double m_theta;
	double m_accuracy;
	size_t m_max_nb_iterations;

	size_t m_nb_nodes;
	size_t m_nb_cells;

	std::vector<double> m_time;

	std::unique_ptr<std::vector<ScalarField>> ptr_rho;
	std::unique_ptr<std::vector<VectorField>> ptr_j, ptr_E, ptr_B;

	std::vector<size_t> m_BC;

public:
	DataKeeper();
	//DataKeeper(std::string file);
	~DataKeeper();

	double get_epsilon() const;
	double get_mu() const;
	double get_sigma() const;
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
	void set_sigma(double sigma);
	void set_t_max(double t_max);
	void set_theta(double theta);
	void set_accuracy(double accuracy);
	void set_max_nb_iterations(size_t max_nb_iterations);

	double get_time(size_t step) const;
	std::vector<double> const& get_time() const;
	void set_time(size_t step, double time);

	double get_rho(size_t t, size_t node_nb) const;
	ScalarField const& get_rho(size_t step) const;
	void set_rho(size_t t, size_t node_nb, double rho);

	Vec3D get_j(size_t t, size_t node_nb) const;
	VectorField const& get_j(size_t step) const;
	void set_j(size_t t, size_t node_nb, Vec3D j);

	Vec3D get_E(size_t t, size_t node_nb) const;
	VectorField const& get_E(size_t step) const;
	void set_E(size_t t, size_t node_nb, Vec3D E);
	void set_E(size_t t, VectorField const& E);
	
	Vec3D get_B(size_t t, size_t node_nb) const;
	VectorField const& get_B(size_t step) const;
	void set_B(size_t t, size_t node_nb, Vec3D B);
	void set_B(size_t t, VectorField const& B);

	void erase_BCs();
	void add_BC(size_t node_nb);
	Vec3D get_boundary_condition_E(size_t node_nb);
	Vec3D get_boundary_condition_B(size_t node_nb);
	bool boundary_node(size_t node_nb) const;
	std::vector<size_t> const& get_dirichlet_nodes() const;

	//void save(std::string file) const;
};