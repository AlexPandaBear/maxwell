#pragma once

#include <vector>

/**
 * A class template representing an unsteady field
 */
template<typename T>
class UnsteadyField
{
private:
	size_t m_nb_steps;
	size_t m_nb_nodes;
	std::vector<T> m_values;

public:
	/**
	 * The constructor of the class
	 *
	 * @param nb_steps The number of time steps the instance can store
	 *
	 * @param nb_nodes The number of nodes in the filed
	 *
	 * @warning N steps means N+1 snapshots of the field stored in the instance
	 */
	UnsteadyField(size_t nb_steps, size_t nb_nodes) :
		m_nb_steps(nb_steps),
		m_nb_nodes(nb_nodes),
		m_values(std::vector<T>((nb_steps+1)*nb_nodes)) {}

	/**
	 * The destructor of the class
	 */
	~UnsteadyField() {}

	/**
	 * A method to reset the sizes of both dimensions
	 *
	 * @param nb_steps The new number of time steps the instance will store
	 *
	 * @param nb_nodes The new number of nodes in the field
	 *
	 * @warning All data stored in the instance will be lost when calling this method
	 *
	 * @warning N steps means N+1 snapshots of the field stored in the instance
	 */
	void reset_size(size_t nb_steps, size_t nb_nodes)
	{
		m_nb_steps = nb_steps;
		m_nb_nodes = nb_nodes;
		m_values = std::vector<T>((nb_steps+1)*nb_nodes);
	}

	/**
	 * A method to get the number of time steps the instance stores
	 *
	 * @returns The number of time steps
	 *
	 * @warning N steps means N+1 snapshots of the field stored in the instance
	 */
	size_t get_nb_steps() const
	{
		return m_nb_steps;
	}

	/**
	 * A method to get the number of nodes in the field represented by the instance
	 *
	 * @returns The number of nodes
	 */
	size_t get_nb_nodes() const
	{
		return m_nb_nodes;
	}

	/**
	 * A method to set the value of the field at a specific time step and a specific node
	 *
	 * @param step The number of the time step at which the value will be stored
	 *
	 * @param node_nb The index of the node which value will be set
	 *
	 * @param value The value to store
	 *
	 * @warning This method will overwrite the previously stored value
	 */
	void set_value(size_t step, size_t node_nb, T value)
	{
		m_values[m_nb_nodes*step + node_nb] = value;
	}

	/**
	 * A method to read the value of the field at a specific time step and a specific node
	 *
	 * @param step The number of the time step at which the value will be read
	 *
	 * @param node_nb The index of the node which value will be read
	 */
	T get_value(size_t step, size_t node_nb) const
	{
		return m_values[m_nb_nodes*step + node_nb];
	}

	/**
	 * A method to get a pointer to the raw data vector
	 *
	 * @returns A pointer to the first element of the std::vector where the values are stored
	 *
	 * @warning Should only be used for Python-C++ interfacing purposes (buffer-protocol)
	 */
	std::vector<T>* get_ptr()
	{
		return &m_values;
	}

	/**
	 * A method to get the memory size of each value of the field
	 *
	 * @returns The memory size of each element of the std::vector where the values are stored
	 *
	 * @warning Should only be used for Python-C++ interfacing purposes (buffer-protocol)
	 */
	size_t get_element_size() const
	{
		return sizeof(T);
	}
};