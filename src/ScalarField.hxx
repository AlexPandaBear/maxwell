#pragma once

#include <vector>

/**
 * A class representing the values of a scalar field
 */
class ScalarField
{
private:
	size_t m_nb_nodes;
	std::vector<double> m_values;

public:
	/**
	 * A constructor of the class creating a nil field
	 *
	 * @param nb_nodes The number of nodes of the field
	 */
	ScalarField(size_t nb_nodes);

	/**
	 * A constructor of the class creating an instance already filled with an uniform value
	 *
	 * @param nb_nodes The number of nodes of the field
	 *
	 * @param value The uniform value of the field
	 */
	ScalarField(size_t nb_nodes, double value);

	/**
	 * The destructor of the class
	 */
	~ScalarField();

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
	 * @param value The value which will be given to the node
	 *
	 * @warning This method will overwrite the previously stored value
	 */
	void set_value(size_t node_nb, double value);

	/**
	 * A method to get the value of the field at a specific node
	 *
	 * @param node_nb The index of the node which value will be read
	 *
	 * @returns The value of the field at this node
	 */
	double get_value(size_t node_nb) const;

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