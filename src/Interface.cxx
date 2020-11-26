/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-std=c++11', '-I/usr/local/lib/python3.7/dist-packages/pybind11-2.4.3-py3.7.egg/']
cfg['dependencies'] = ['SimManager.hxx', 'UnsteadyMaxwellKernel.hxx', 'Mesh3D.hxx', 'Cell.hxx', 'DataKeeper.hxx', 'DataProcessor.hxx', 'Vec3D.hxx', 'ScalarField.hxx', 'VectorField.hxx', 'Matrix.hxx', 'MatrixFEM.hxx', 'ElectrostaticSimManager.hxx', 'ElectrostaticKernel.hxx', 'ElectrostaticDataKeeper.hxx']
cfg['sources'] = ['SimManager.cxx', 'UnsteadyMaxwellKernel.cxx', 'Mesh3D.cxx', 'Cell.cxx', 'DataKeeper.cxx', 'DataProcessor.cxx', 'Vec3D.cxx', 'ScalarField.cxx', 'VectorField.cxx', 'Matrix.cxx', 'MatrixFEM.cxx', 'ElectrostaticSimManager.cxx', 'ElectrostaticKernel.cxx', 'ElectrostaticDataKeeper.cxx']
%>
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "SimManager.hxx"
#include "ElectrostaticSimManager.hxx"
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
		.def("defineInitialState", &SimManager::define_initial_state,
			py::arg("E0"),
			py::arg("B0"))
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
		.def("getNodeXYZ", &SimManager::get_node_xyz,
			py::arg("node_nb"))
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
			py::arg("step"))
		.def("getPoyntingVectorNorm", &SimManager::get_poynting_vector_norm,
			py::arg("step"));

	py::class_<ElectrostaticSimManager>(m, "ESM")
		.def(py::init<bool>(),
			py::arg("verbose") = true)
		.def("setEpsilon0", &ElectrostaticSimManager::set_epsilon0,
			py::arg("epsilon0"))
		.def("setAccuracy", &ElectrostaticSimManager::set_accuracy,
			py::arg("accuracy"))
		.def("setMaxNbIterations", &ElectrostaticSimManager::set_max_nb_iterations,
			py::arg("max_nb_iterations"))
		.def("generateCubeMesh", &ElectrostaticSimManager::generate_cube_mesh,
			py::arg("x_min"),
			py::arg("x_max"),
			py::arg("nx"),
			py::arg("y_min"),
			py::arg("y_max"),
			py::arg("ny"),
			py::arg("z_min"),
			py::arg("z_max"),
			py::arg("nz"))
		.def("setRho", &ElectrostaticSimManager::set_rho,
			py::arg("node_id"),
			py::arg("rho_field"))
		.def("setRhoField", &ElectrostaticSimManager::set_rho_field,
			py::arg("rho_field"))
		.def("setEpsilonRField", &ElectrostaticSimManager::set_epsilon_r_field,
			py::arg("epsilon_r_field"))
		.def("simulate", &ElectrostaticSimManager::simulate)
		.def("getNbNodes", &ElectrostaticSimManager::get_nb_nodes)
		.def("getNodeXYZ", &ElectrostaticSimManager::get_node_xyz,
			py::arg("node_id"))
		.def("getMesh", &ElectrostaticSimManager::get_mesh)
		.def("getRho", &ElectrostaticSimManager::get_rho)
		.def("getEpsilon", &ElectrostaticSimManager::get_epsilon)
		.def("getPhi", &ElectrostaticSimManager::get_phi)
		.def("getE", &ElectrostaticSimManager::get_E);

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
		.def(py::init<size_t>(),
			py::arg("nb_nodes"))
		.def(py::init<size_t, double>(),
			py::arg("node_nb"),
			py::arg("value"))
		.def("getValue", &ScalarField::get_value,
			py::arg("node_nb"))
		.def("setValue", &ScalarField::set_value,
			py::arg("node_nb"),
			py::arg("value"))
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
		.def(py::init<size_t>(),
			py::arg("nb_nodes"))
		.def("getVector", &VectorField::get_value,
			py::arg("node_nb"))
		.def("setVector", py::overload_cast<size_t, Vec3D>(&VectorField::set_value),
			py::arg("node_nb"),
			py::arg("vector"))
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