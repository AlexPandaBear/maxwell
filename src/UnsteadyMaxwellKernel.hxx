#pragma once

#include "Mesh3D.hxx"
#include "DataKeeper.hxx"
#include "MatrixFEM.hxx"

class UnsteadyMaxwellKernel
{
private:
	const double m_eps, m_mu, m_sigma;

	const size_t m_nb_steps;
	const double m_t_max;
	const double m_dt;
	const double m_theta;
	const double m_accuracy;
	const size_t m_max_nb_iterations;

	const size_t m_nb_nodes, m_nb_cells;

	const Mesh3D& m_mesh;
	DataKeeper& m_data;

	MatrixFEM mat_L, mat_R;
	std::vector<double> vec_X, vec_B;

	bool check_divergence_in_ID() const;

	void build_matrices();

	void initialize_vector_X();

	void display_progression(size_t step) const;

public:
	UnsteadyMaxwellKernel(Mesh3D& mesh, DataKeeper& data);
	~UnsteadyMaxwellKernel();

	void simulate();	
};