#pragma once

#include <DataKeeper.hxx>
#include <SparseMatrix.hxx>

class UnsteadyMaxwellKernel
{
private:
	const double m_eps, m_mu;

	const size_t m_nb_steps;
	const double m_t_max;
	const double m_theta;
	const double m_accuracy;
	const size_t m_max_nb_iterations;

	const size_t m_nb_nodes, m_nb_trgls;

	const Mesh3D& m_mesh;
	DataKeeper& m_data;

	SparseMatrix mat_A;
	std::vector<double> vec_X, vec_B;

	void check_ID_BC_compatibility() const;
	void check_divergence_in_ID() const;

	bool boundary_node(size_t node_nb) const;

	void build_matrix_A();
	void build_vector_B(size_t t);

	double GS_step(size_t i) const;
	double GS_iteration();

public:
	UnsteadyMaxwellKernel(Mesh3D& mesh, DataKeeper& data);
	~UnsteadyMaxwellKernel();

	void simulate();	
};