/**
 * @file shockTube.cpp
 * @short Sod shock tube simulation
 * @date 18.11.2025
 * @author Dominik Wilde, Philipp Spelten, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include <stdlib.h>
#include <sstream>

#include "deal.II/numerics/data_out.h"
#include "deal.II/base/utilities.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/CFDSolverUtilities.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/benchmarks/ShearLayer2D.h"
#include "natrium/utilities/Info.h"
#include "natrium/utilities/CommandLineParser.h"

#include "SodShockTube.h"

using namespace natrium;

// Main function
int main(int argc, char** argv) {

  MPIGuard::getInstance(argc, argv);

  // ========================================================================
  // READ COMMAND LINE PARAMETERS
  // ========================================================================

  CommandLineParser parser(argc, argv);
  parser.addDocumentationString("shockTube", "Shocktube as described by Sod (1978)");
  parser.setArgument<int>("ref-level", "refinement of the computational grid", 0);
  parser.setArgument<double>("length", "length in x direction", 1);
  parser.setArgument<int>("nx", "number of cells in x-direction", 25);  // p=4 -> 100 grid points
  parser.setArgument<double>("tx", "transformation of the grid in x-direction (<1)", 0);
  parser.setArgument<double>("ty", "transformation or the grid in y-direction (<1)", 0);
  parser.setArgument<int>("filter", "apply filtering", 0);
  parser.setArgument<int>("filter-s", "parameter as filter", 32);
  parser.setArgument<int>("vmult", "apply vMultLimiter", 0);
  parser.setArgument<double>("visc","viscosity of the fluid",0.001);

  try {
    parser.importOptions();
  } catch (HelpMessageStop&) {
    return 0;
  }

  // ========================================================================
  // MAKE FLOW PROBLEM
  // ========================================================================

  double perturbation = 0.05;
  double kappa = 80;

  // double Ma = 0.04 / (1.0 / sqrt(3));
  // double Re;
  double u0;
  // double scaling = sqrt(3) * u0 / Ma
  double scaling = 1.0;
  double scaled_viscosity = parser.getArgument<double>("length") * parser.getArgument<double>("visc");

  boost::shared_ptr<ProblemDescription<2> > shockTube = boost::make_shared<SodShockTube>(
    parser.getArgument<double>("length"),
    scaled_viscosity,
    parser.getArgument<int>("ref-level"),
    u0, kappa,
    parser.getArgument<int>("nx"),
    perturbation,
    parser.getArgument<double>("tx"),
    parser.getArgument<double>("ty")
  );

  double t_max = parser.getArgument<double>("length")*sqrt(3.0)*0.15; //for t_phys = 0.15

  // ========================================================================
  // CONFIGURE SOLVER
  // ========================================================================

  boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<SolverConfiguration>();
  configuration->setSwitchOutputOff(false);
  configuration->setUserInteraction(false);
  configuration->setCommandLineVerbosity(ALL);
  configuration->setOutputTableInterval(10);
  configuration->setOutputSolutionInterval(10);
  configuration->setOutputCheckpointInterval(1e9);
  configuration->setOutputGlobalTurbulenceStatistics(true);

  configuration->setConvergenceThreshold(1e-10);
  configuration->setSimulationEndTime(t_max);
  configuration->setCFL(1);
  configuration->setPrandtlNumber(1.0);
	configuration->setHeatCapacityRatioGamma(1.4);

  configuration->setStencilScaling(scaling);
  configuration->setStencil(Stencil_D2Q25H);
  configuration->setSedgOrderOfFiniteElement(4);
  configuration->setCollisionScheme(BGK_STANDARD);
  configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);
  configuration->setAdvectionScheme(SEMI_LAGRANGIAN);

  configuration->setExponentialFilterAlpha(36);
  configuration->setExponentialFilterNc(3);

  //parser.applyToSolverConfiguration(*configuration);
  //configuration->setFiltering(true);
  //configuration->setFilteringScheme(EXmathPONENTIAL_FILTER);
  configuration->setVmultLimiter(bool(parser.getArgument<int>("vmult")));
  pout << "VMultLimiter is " << configuration->isVmultLimiter() << endl;
  std::stringstream dirname;
  dirname << getenv("NATRIUM_HOME") << "/shockTube";
  if (parser.hasArgument("minion-brown")) {
    dirname << "-MinionBrown";
  }
  dirname << "/N" << parser.getArgument<int>("ref-level")*2.0*parser.getArgument<int>("nx") << "-p"
    << configuration->getSedgOrderOfFiniteElement() << "-sl"
    << static_cast<int>(configuration->getAdvectionScheme()) << "-coll"
    << static_cast<int>(configuration->getCollisionScheme()) << "-int"
    << static_cast<int>(configuration->getTimeIntegrator()) << "_"
    << static_cast<int>(configuration->getDealIntegrator()) << "-CFL"
    << configuration->getCFL() << "-reg" << static_cast<int>(configuration->getRegularizationScheme())<< "-scaling"
    << configuration->getStencilScaling() << "-suppP"
    << configuration->getSupportPoints() << "-vMult"
    << configuration->isVmultLimiter() << "-visc"
    << parser.getArgument<double>("visc");
  if (parser.getArgument<int>("filter") != 0) {
    dirname << "-filter" << parser.getArgument<int>("filter") << "-filt_s"
      << parser.getArgument<int>("filter-s");
  }
  if ((parser.getArgument<double>("tx") != 0)
      or (parser.getArgument<double>("ty") != 0)) {
    dirname << "-tx" << parser.getArgument<double>("tx") << "-ty"
      << parser.getArgument<double>("ty");
  }
  if (configuration->getRegularizationScheme() != NO_REGULARIZATION){
    dirname << "-reg" << static_cast<int>(configuration->getRegularizationScheme());
  }
  configuration->setOutputDirectory(dirname.str());

  parser.applyToSolverConfiguration(*configuration);
  pout << "Simulation end time will be t_max = " << t_max << endl;

  // ========================================================================
  // RUN SOLVER
  // ========================================================================

  natrium::CompressibleCFDSolver<2> solver(configuration, shockTube);
  solver.run();

  // ========================================================================
  // FINAL OUTPUT
  // ========================================================================
  
  pout << "Simulation successful." << endl;
  return 0;
}
