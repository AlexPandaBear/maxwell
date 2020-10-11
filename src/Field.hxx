#pragma once

#include <vector>

template<typename T>
class Field
{
private:
	size_t m_nb_nodes;
	std::vector<T> m_values;

public:
	Field() {}

	~Field() {}

	void reset_nb_nodes(size_t nb_nodes)
	{
		m_nb_nodes = nb_nodes;
		m_values = std::vector<T>(nb_nodes);
	}

	void set_value(size_t node_nb, T value)
	{
		m_values[node_nb] = value;
	}

	T get_value(size_t node_nb) const
	{
		return m_values[node_nb];
	}	
};