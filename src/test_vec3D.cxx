#include "Vec3D.hxx"

int main(int argc, char const *argv[])
{
	Vec3D U(1., 2., 3.), V(-1., 0., 1.);
	
	std::cout << "U = ";
	U.display();

	std::cout << "V = ";
	V.display();

	std::cout << "Testing dot product: U.V = " << Vec3D::dot_product(U, V) << std::endl;

	std::cout << "Testing cross product: U x V = ";
	Vec3D::cross_product(U, V).display();


	return 0;
}