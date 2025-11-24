## Calculating dt from CFL

1. [CFDSolver](../../library/natrium/solver/CFDSolver.cpp):245
```
	// set time step size
	double delta_t = CFDSolverUtilities::calculateTimestep<dim>(
			*(m_problemDescription->getMesh()),
			m_configuration->getSedgOrderOfFiniteElement(), *m_stencil,
			configuration->getCFL());
```
2. [CFDSolverUtilities](../../library/natrium/utilities/CFDSolverUtilities.cpp):169
```
template<size_t dim>
double CFDSolverUtilities::calculateTimestep(
		const Mesh<dim>& tria,
		const size_t orderOfFiniteElement,
		const Stencil& stencil,
		double cFL) {
	assert(orderOfFiniteElement >= 1);
	double dx = CFDSolverUtilities::getMinimumVertexDistance<dim>(tria);
	double u = stencil.getMaxParticleVelocityMagnitude();
	// according to Hesthaven, dt ~ p^{-2}
	double dt = cFL * dx / (u * orderOfFiniteElement * orderOfFiniteElement);
	return dt;
}
```
3. [CFDSolverUtilities](../../library/natrium/utilities/CFDSolverUtilities.cpp):97
```
template<size_t dim>
double CFDSolverUtilities::getMinimumVertexDistance(const Mesh<dim>& tria) {
	// calculate minimal distance between vertices of the triangulation
	double min_vertex_distance = 100000000000.0;
	double distance;
	for (typename Mesh<dim>::active_cell_iterator cell = tria.begin_active(); cell != tria.end(); ++cell) {
		if (cell->is_locally_owned()) {
			distance = cell->minimum_vertex_distance();
            min_vertex_distance = std::min(min_vertex_distance, distance);
		}
	} // sync over all MPI processes
	return dealii::Utilities::MPI::min_max_avg(min_vertex_distance, MPI_COMM_WORLD).min;
}
```
= 0.00278996
4. [D2Q19H.h](../../library/natrium/stencils/D2Q19H.h):97
```
	virtual double getMaxParticleVelocityMagnitude() const {
        return sqrt(2)*m_scaling;
	}
```
5. [step-grid-in.cpp](../../examples/step-grid-in/step-grid-in.cpp):63
```
    double scaling = sqrt(3) * U / (Ma * sqrt(gamma*reference_temperature));
```
= 0.975900

--> $$dt = \frac{(dx_{Cell}/p^2)*CFL}{(\sqrt{2}*(\sqrt{3}*U*Ma*\sqrt{\gamma T_0}))} = 0.00278996*4/(\sqrt{2}*0.975900*4*4) = 0.0005053796$$
