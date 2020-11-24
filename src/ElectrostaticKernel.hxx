#pragma once

#include <stdexcept>
#include "Mesh3D.hxx"
#include "ElectrostaticDataKeeper.hxx"
#include "MatrixFEM.hxx"

class ElectrostaticKernel
{
private:
	const Mesh3D& m_mesh;
	ElectrostaticDataKeeper& m_data;

	Matrix mat_A;
	std::vector<double> vec_X, vec_B;

	void build_matrix_A();
	void build_vector_B();

public:
	ElectrostaticKernel(Mesh3D& mesh, ElectrostaticDataKeeper data);
	~ElectrostaticKernel();

	void simulate_phi();
	void derive_E();	
};