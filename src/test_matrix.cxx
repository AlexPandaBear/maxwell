#include "Matrix.hxx"

int main(int argc, char const *argv[])
{
	std::cout << "Testing matricial system resolution..." << std::endl;

	Matrix A(10);
	std::vector<double> X(10, 0.), B(10, 10.);

	for (size_t i = 0; i < 10; i++)
	{
		for (size_t j = 0; j < 10; j++)
		{
			A.set_coef(i, j, 1.);
		}

		A.set_coef(i, i, 11.);
	}

	std::cout << "\nInitial state of the system:" << std::endl;
	A.display_matricial_system(X, B);

	double residual;

	do
	{
		residual = A.perform_gauss_seidel_iteration(X, B);
	}
	while (residual > 0.00001);

	std::cout << "\nSystem after Gauss-Seidel resolution:" << std::endl;
	A.display_matricial_system(X, B);




	std::cout << "\nTesting matrix vector product..." << std::endl;

	std::vector<double> R(10, 0.);
	A.perform_matrix_vector_product(X, R);
	A.display_matricial_system(X, R);

	return 0;
}