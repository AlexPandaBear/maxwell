/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-std=c++11', '-I/usr/local/lib/python3.7/dist-packages/pybind11-2.4.3-py3.7.egg/']
cfg['dependencies'] = ['SimManager.hxx', 'UnsteadyMaxwellKernel.hxx', 'SparseMatrix.hxx', 'Mesh3D.hxx', 'Cell.hxx', 'DataKeeper.hxx', 'DataProcessor.hxx', 'Vec3D.hxx', 'ScalarField.hxx', 'VectorField.hxx']
cfg['sources'] = ['SimManager.cxx', 'UnsteadyMaxwellKernel.cxx', 'SparseMatrix.cxx', 'Mesh3D.cxx', 'Cell.cxx', 'DataKeeper.cxx', 'DataProcessor.cxx', 'Vec3D.cxx', 'ScalarField.cxx', 'VectorField.cxx']
%>
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "SimManager.hxx"
#include "Vec3D.hxx"

namespace py = pybind11;

PYBIND11_MODULE(_maxwell, m)
{
	py::class_<SimManager>(m, "SM")
		.def(py::init<>())
		.def("setConstants", &SimManager::set_constants,
			py::arg("epsilon"),
			py::arg("mu"),
			py::arg("sigma"))
		.def("setSimulationParameters", &SimManager::set_simulation_parameters,
			py::arg("t_max"),
			py::arg("nb_steps"),
			py::arg("theta"),
			py::arg("accuracy"),
			py::arg("max_nb_iterations"),
			py::arg("x_min"),
			py::arg("x_max"),
			py::arg("y_min"),
			py::arg("y_max"),
			py::arg("z_min"),
			py::arg("z_max"),
			py::arg("nx"),
			py::arg("ny"),
			py::arg("nz"))
		.def("defineInitialStateBackgroundValues", &SimManager::define_initial_state_background_values,
			py::arg("rho0"),
			py::arg("j0"),
			py::arg("E0"),
			py::arg("B0"))
		.def("addWire", &SimManager::add_wire,
			py::arg("wire_skeleton"),
			py::arg("wire_radius"),
			py::arg("wire_current"))
		.def("addBoundaryCondition", &SimManager::add_boundary_condition,
			py::arg("node_nb"))
		.def("lockAllBoundaryNodes", &SimManager::lock_all_boundary_nodes)
		.def("simulate", &SimManager::simulate)
		//.def("save", &SimManager::save,
		//	py::arg("simulation_name"))
		.def("getNbNodes", &SimManager::get_nb_nodes)
		.def("getNbCells", &SimManager::get_nb_cells)
		.def("getNodeIds", &SimManager::get_node_ids)
		.def("getMesh", &SimManager::get_mesh)
		.def("getTime", &SimManager::get_time)
		.def("getRho", &SimManager::get_rho,
			py::arg("step"))
		.def("getJ", &SimManager::get_j,
			py::arg("step"))
		.def("getE", &SimManager::get_E,
			py::arg("step"))
		.def("getB", &SimManager::get_B,
			py::arg("step"))
		.def("getEnergyDensity", &SimManager::get_energy_density,
			py::arg("step"))
		.def("getPoyntingVector", &SimManager::get_poynting_vector,
			py::arg("step"));

	py::class_<Vec3D>(m, "Vec3D")
		.def(py::init<>())
		.def(py::init<double, double, double>(),
			py::arg("x"),
			py::arg("y"),
			py::arg("z"))
		.def("getX", &Vec3D::get_x)
		.def("getY", &Vec3D::get_y)
		.def("getZ", &Vec3D::get_z)
		.def("setX", &Vec3D::set_x,
			py::arg("x"))
		.def("setY", &Vec3D::set_y,
			py::arg("y"))
		.def("setZ", &Vec3D::set_z,
			py::arg("z"))
		.def("set", &Vec3D::set,
			py::arg("x"),
			py::arg("y"),
			py::arg("z"));

	py::class_<ScalarField>(m, "ScalarField", py::buffer_protocol())
		.def_buffer([](ScalarField& f) -> py::buffer_info
			{
				return py::buffer_info(
						f.get_ptr(),
						f.get_element_size(),
						py::format_descriptor<double>::format(),
						1,
						py::detail::any_container<long int>({f.get_nb_nodes()}),
						py::detail::any_container<long int>({f.get_element_size()})
					);
			});

	py::class_<VectorField>(m, "VectorField", py::buffer_protocol())
		.def("getVector", &VectorField::get_value,
			py::arg("node_nb"))
		.def_buffer([](VectorField& f) -> py::buffer_info
			{
				return py::buffer_info(
						f.get_ptr(),
						f.get_element_size(),
						py::format_descriptor<double>::format(),
						2,
						py::detail::any_container<long int>({f.get_nb_nodes(), 3}),
						py::detail::any_container<long int>({f.get_element_size()*3, f.get_element_size()})
					);
			});
}