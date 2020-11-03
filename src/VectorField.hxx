#pragma once

#include <vector>
#include "Vec3D.hxx"

/**
 * A class representing the values of a vector field
 */
class VectorField
{
private:
	size_t m_nb_nodes;
	std::vector<double> m_values;

public:
	/**
	 * The constructor of the class
	 *
	 * @param nb_nodes The number of nodes of the field
	 */
	VectorField(size_t nb_nodes);

	/**
	 * The destructor of the class
	 */
	~VectorField();

	/**
	 * A method to reset the number of nodes of the field
	 *
	 * @param nb_nodes The new number of nodes for the field
	 *
	 * @warning All data previously stored in the instance will be lost
	 */
	void reset_nb_nodes(size_t nb_nodes);

	/**
	 * A method to get the number of nodes of the field
	 *
	 * @returns the number of values stored in this instance
	 */
	size_t get_nb_nodes() const;

	/**
	 * A method to set the value at a specific node
	 *
	 * @param node_nb The index of the node which value will be set
	 *
	 * @param x_value The value along the x-axis which will be given to the node
	 *
	 * @param y_value The value along the y-axis which will be given to the node
	 *
	 * @param z_value The value along the z-axis which will be given to the node
	 *
	 * @warning This method will overwrite the previously stored value
	 */
	void set_value(size_t node_nb, double x_value, double y_value, double z_value);

	/**
	 * A method to set the value at a specific node
	 *
	 * @param node_nb The index of the node which value will be set
	 *
	 * @param value The value which will be given to the node
	 *
	 * @warning This method will overwrite the previously stored value
	 */
	void set_value(size_t node_nb, Vec3D value);

	/**
	 * A method to get the value along the x-axis of the field at a specific node
	 *
	 * @param node_nb The index of the node which value will be read
	 *
	 * @returns The value along the x-axis of the field at this node
	 */
	double get_x_value(size_t node_nb) const;

	/**
	 * A method to get the value along the y-axis of the field at a specific node
	 *
	 * @param node_nb The index of the node which value will be read
	 *
	 * @returns The value along the y-axis of the field at this node
	 */
	double get_y_value(size_t node_nb) const;

	/**
	 * A method to get the value along the z-axis of the field at a specific node
	 *
	 * @param node_nb The index of the node which value will be read
	 *
	 * @returns The value along the z-axis of the field at this node
	 */
	double get_z_value(size_t node_nb) const;

	/**
	 * A method to get access to the value of the field at a specific node as a Vec3D object
	 *
	 * @param node_nb The index of the node which value will be read
	 *
	 * @returns A Vec3D object with access to the value of the field at this node
	 */
	Vec3D get_value(size_t node_nb) const;

	/**
	 * A method to get a pointer to the raw data vector
	 *
	 * @returns A pointer to the first element of the std::vector where the values are stored
	 *
	 * @warning Should only be used for Python-C++ interfacing purposes (buffer-protocol)
	 */
	double* get_ptr();

	/**
	 * A method to get the memory size of each value of the field
	 *
	 * @returns The memory size of each element of the std::vector where the values are stored
	 *
	 * @warning Should only be used for Python-C++ interfacing purposes (buffer-protocol)
	 */
	size_t get_element_size() const;
};