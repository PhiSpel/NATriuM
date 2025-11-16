/**
 * @file step-sphere.cpp
 * @short Flow around a sphere
 * @date 15.02.2022
 * @author Dominik Wilde, UC San Diego
 */

#include <stdlib.h>
#include <sstream>
#include <mpi.h>

#include "deal.II/numerics/data_out.h"

#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/CompressibleCFDSolver.h"

#include "natrium/solver/SolverConfiguration.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/CommandLineParser.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/stencils/D3Q45.h"

#include "natrium/utilities/BasicNames.h"

#include "Sphere.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

    MPIGuard::getInstance();
    CommandLineParser parser(argc, argv);
    parser.setArgument<double>("Ma", "Mach number", 1.05); // 1.2
    parser.setArgument<int>("Re", "Reynolds number", 1000);  // 300
    parser.setArgument<int>("ref-level", "Refinement level", 4);
    parser.setArgument<int>("server-end", "Maximum server time [s]", 30000);
    //parser.setArgument<double>("tmax", "Maximum physical time [s]", 25);

    try {
        parser.importOptions();
    } catch (HelpMessageStop&){
        return 0;
    }


	pout << "Starting NATriuM step-sphere..." << endl;

    const int refLevel = parser.getArgument<int>("ref-level");

	// set Reynolds and Mach number
	const double Ma = parser.getArgument<double>("Ma");
	const double gamma = 1.4;
	const double T = 1.0;
	// increase velocity to gain correct speed
    const double Re = parser.getArgument<int>("Re");;

    // set Problem so that the right Re and Ma are achieved
    double U = 1;
    double scaling = 1;
    scaling /= sqrt(gamma);
    scaling*=sqrt(3)/Ma;
    const double viscosity = U / Re; // (because L = 1)


	// load grid
	boost::shared_ptr<Sphere> sphere = boost::make_shared<Sphere>(
				viscosity, U, refLevel);

	// configure solver
	boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<
			SolverConfiguration>();

	configuration->setUserInteraction(false);
	configuration->setOutputCheckpointInterval(10000);
	configuration->setOutputSolutionInterval(100);
	configuration->setOutputTableInterval(100);
	configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
	configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);
	configuration->setStencil(Stencil_D3Q45);
	configuration->setMachNumber(Ma);
	configuration->setSupportPoints(GAUSS_LOBATTO_CHEBYSHEV_POINTS);
	configuration->setCollisionScheme(BGK_STANDARD);
	configuration->setServerEndTime(parser.getArgument<int>("server-end"));

	configuration->setStencilScaling(scaling);
	configuration->setHeatCapacityRatioGamma(gamma);
	configuration->setPrandtlNumber(0.71);

    std::stringstream dirname;
    dirname << getenv("NATRIUM_HOME") << "/Sphere/refLevel" << refLevel << "-p" << configuration->getSedgOrderOfFiniteElement() << "-Ma" << Ma  << "-Re" << Re << "-CFL" << configuration->getCFL();
    configuration->setOutputDirectory(dirname.str());

	parser.applyToSolverConfiguration(*configuration);

	boost::shared_ptr<ProblemDescription<3> > sphereProblem= sphere;


        CompressibleCFDSolver<3> solver(configuration, sphereProblem);
        solver.run();

	if (is_MPI_rank_0()) LOG(WELCOME) << "NATriuM step-sphere terminated." << endl;

	return 0;
}
