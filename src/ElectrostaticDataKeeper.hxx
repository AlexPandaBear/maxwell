#pragma once

#include "ScalarField.hxx"
#include "VectorField.hxx"

class ElectrostaticDataKeeper
{
private:
	double m_eps0;
	size_t m_nb_nodes;
	double m_accuracy;
	size_t m_max_nb_iterations;

	ScalarField m_rho, m_eps;
	VectorField m_E;

public:
	ElectrostaticDataKeeper();
	~ElectrostaticDataKeeper();

	size_t get_nb_nodes() const;
	void set_nb_nodes(size_t nb_nodes);

	double get_epsilon0() const;
	void set_epsilon0(double epsilon0);

	size_t get_max_nb_iterations() const;
	void set_max_nb_iterations(size_t max_nb_iterations);

	double get_accuracy() const;
	void set_accuracy(double accuracy);	

	double get_rho(size_t node_nb) const;
	double get_epsilon(size_t node_nb) const;
	Vec3D get_E(size_t node_nb) const;

	ScalarField const& get_rho() const;
	ScalarField const& get_epsilon() const;
	VectorField const& get_E() const;

	void set_rho(size_t node_nb, double rho);
	void set_epsilon_r(size_t node_nb, double epsilon_r);
	void set_E(size_t node_nb, Vec3D const& E);
};