/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-std=c++11', '-I/usr/local/lib/python3.7/dist-packages/pybind11-2.4.3-py3.7.egg/']
cfg['dependencies'] = ['SimManager.hxx', 'UnsteadyMaxwellKernel.hxx', 'SparseMatrix.hxx', 'Mesh3D.hxx', 'BoundaryCondition.hxx', 'Cell.hxx', 'Node.hxx', 'DataKeeper.hxx', 'Field.hxx']
cfg['sources'] = ['SimManager.cxx', 'UnsteadyMaxwellKernel.cxx', 'SparseMatrix.cxx', 'Mesh3D.cxx', 'BoundaryCondition.cxx', 'Cell.cxx', 'Node.cxx', 'DataKeeper.cxx']
%>
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "SimManager.hxx"

namespace py = pybind11;

PYBIND11_MODULE(_maxwell, m)
{
	py::class_<SimManager>(m, "SM")
		.def(py::init<>())
		.def("setConstants", &SimManager::set_constants,
			py::arg("epsilon"),
			py::arg("mu"))
		.def("setSimulationParameters", &SimManager::set_simulation_parameters,
			py::arg("t_max"),
			py::arg("nb_steps"),
			py::arg("theta"),
			py::arg("accuracy"),
			py::arg("max_nb_iterations"),
			py::arg("nb_nodes"),
			py::arg("nb_cells"))
		.def("defineInitialState", &SimManager::define_initial_state,
			py::arg("rho0"),
			py::arg("j0"),
			py::arg("E0"),
			py::arg("B0"))
		.def("addBoundaryCondition", &SimManager::add_boundary_condition,
			py::arg("node_nb"),
			py::arg("E"),
			py::arg("B"))
		.def("simulate", &SimManager::simulate)
		.def("save", &SimManager::save,
			py::arg("simulation_name"));

	py::class_<Vec3D>(m, "Vec3D")
		.def(py::init<>())
		.def("getX", &Vec3D::get_x)
		.def("getY", &Vec3D::get_y)
		.def("getZ", &Vec3D::get_z)
		.def("setX", &Vec3D::set_x,
			py::arg("x"))
		.def("setY", &Vec3D::set_y,
			py::arg("y"))
		.def("setZ", &Vec3D::set_z,
			py::arg("z"));
}