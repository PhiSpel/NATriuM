//
// Created by dwilde3m on 02.12.21.
//
/*
 * step-mixingLayer.cpp
 *
 *  Created on: Dec 02, 2021
 *      Author: dominik
 */
#include <fstream>
//#include <time.h>
#include <stdlib.h>
//#include "deal.II/numerics/data_out.h"
#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/solver/SolverConfiguration.h"
#include "natrium/stencils/Stencil.h"
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/CommandLineParser.h"
#include "MixingLayer3D.h"
#include "ShearLayerStats.h"
#include <unistd.h>
#include <mpi.h>


using namespace natrium;

double shearLayerThickness = 0.093;
//int dft_points = 20;

// Main function
int main(int argc, char** argv) {
    MPIGuard::getInstance(argc, argv);
    pout << "Starting NATriuM step-mixingLayer ..." << endl;
    /////////////////////////////////////////////////
    // read from command line
    //////////////////////////////////////////////////
    CommandLineParser parser(argc, argv);
    parser.setArgument<int>("Re", "Reynolds number 1/nu", 800);
    /* TODO: Convective Mach number Mc=(U1-Uc)/c1, Uc=(U1c2+U2c1)/(c1+c2)
    Uc is the convective velocity of the large structures, U\
    and U2 are the freestream velocities, and c{ and c2 are the
    freestream sound speeds.
    Set to 0.3, 0.7, 0.9, 1.0, 1.2
    */
    parser.setArgument<double>("Ma", "Mach number", 0.3);
    parser.setArgument<double>("time", "simulation time (s)", 15);
    parser.setArgument<double>("randuscaling", "factor to scale random velocity field", 5);
    parser.setArgument<double>("uscaling", "factor to scale U1, i.e. deltaUx", 1);
    parser.setArgument<double>("CFL", "CFL number. Should be between 0.4 and 2", 1);
    parser.setArgument<double>("gamma", "Heat capacity ratio. Should be 1.4", 1.4);
    parser.setArgument<double>("ref-temp", "Reference temperature. Should be between 0.85 and 1 (lower may be more stable).", 1);
    parser.setArgument<double>("lx", "Half length in x-direction (multiples of deltaTheta0)", 150);
    parser.setArgument<double>("ly", "Half length in y-direction (multiples of deltaTheta0)", 75);
    parser.setArgument<double>("lz", "Half length in z-direction (multiples of deltaTheta0)", 40);
    parser.setArgument<int>("nout", "output vtk every nout steps", 2000);
    parser.setArgument<int>("ncheckpoint", "output checkpoint every ncheckpoint steps", 20000);
    parser.setArgument<int>("n-no-out", "do not output vtk before iteration n-no-out", -1);
    parser.setArgument<int>("nstats", "output stats every nstats steps", 20);
    parser.setArgument<string>("meshname", "name of the mesh file (shearlayer_*.txt)", "cube");
    parser.setArgument<string>("randuname", "name of the initial velocity file (random_u_*.txt)", "cube_k048_half");
    parser.setArgument<string>("bc", "Boundary condition. Choose between 'EQ_BC' (equilibrium), 'DN_BC' (do nothing),"
                                     "'FOBB_BC' (First Order Bounce Back),'ThBB_BC' (Thermal Bounce Back), 'VNeq_BC' (Velocity Non-Equilibrium Bounce Back),"
                                     "'PP_BC' (Periodic - meh)", "EQ_BC");
    parser.setArgument<int>("order", "order of finite elements", 3);
    parser.setArgument<int>("ref-level", "Refinement level of the computation grid.", 4);
    parser.setArgument<int>("grid-repetitions",
                            "Number of grid cells along each axis before global refinement; "
                            "to produce grids with refinements that are not powers of two.", 1);
    parser.setArgument<int>("restart", "Restart at iteration ...", 0);
    parser.setArgument<int>("server-end", "Maximum server time [s]", 82800);
    parser.setArgument<int>("rep-x", "Number of repetitions in x-direction (to refine the grid in steps that are not 2^N).", 8);
    parser.setArgument<int>("rep-y", "cf. rep-x", 4);
    parser.setArgument<int>("rep-z", "cf. rep-x", 2);
    parser.setArgument<double>("center", "Central part with high-res grid, choose between 0.1 and 1", 0.7);
    parser.setArgument<double>("dy-scaling", "scale dy to dy-scaling-times the element size (<1 to refine boundaries, >1 to loosen, 1 for equidistant mesh)", 3);

    try { parser.importOptions();
    } catch (HelpMessageStop&) { return 0;
    }
    auto meshname = parser.getArgument<string>("meshname");
    auto randuname = parser.getArgument<string>("randuname");
    auto bc = parser.getArgument<string>("bc");
    if ((bc != "DN_BC") and (bc != "EQ_BC") and (bc != "FOBB_BC") and (bc != "ThBB_BC") and (bc != "VNeq_BC") and (bc != "PP_BC")) {
        if (is_MPI_rank_0()) LOG(BASIC) << "Invalid boundary condition option! Fallback to default (EQ_BC)." << endl << endl;
        bc = "EQ_BC";
    }
    double randuscaling = parser.getArgument<double>("randuscaling");
    double uscaling = parser.getArgument<double>("uscaling");
    double Re = parser.getArgument<int>("Re");
    // Grid resolution
    const int ref_level = parser.getArgument<int>("ref-level");
    std::vector<unsigned int> repetitions(3);
    repetitions.at(0) = parser.getArgument<int>("rep-x");
    repetitions.at(1) = parser.getArgument<int>("rep-y");
    repetitions.at(2) = parser.getArgument<int>("rep-z");

    long nout = parser.getArgument<int>("nout");
    auto time = parser.getArgument<double>("time");
    const int restart = parser.getArgument<int>("restart");
    if ((restart > 0) and is_MPI_rank_0()) {
        LOG(WELCOME) << "==================================================="
                     << endl << "=== Starting NATriuM step-mixingLayer... ===="
                     << endl << "=== Restart iteration: " << restart << endl
                     << "===================================================" << endl;
    }

    /////////////////////////////////////////////////
    // set parameters, set up configuration object
    //////////////////////////////////////////////////
    // im Paper von Gassner und Beck ist U = 1/2pi definiert !!!!!
    // Aber ihre zeitangaben beziehen sich auf U = 1, wie bei Brachet (1991)
    // hier simulieren wir jetzt  U = 1 (ist im TGV3D modul sowieso nur so definiert)
    // und mit der Reynoldsähnlichkeit ist das kein Problem, da wir gegenüber Gassner
    // und Beck ja auch die Viskosität verändern
    // Nach einem Blick in van Rees et.al. (2011) und einer erfolgreichen Simulation in Palabos:
    // (Hier war tau = 3*nu_LB + 0.5, nu_LB = U_lattice * (N/2pi) / Re, dt = 2pi/N* U_lattice
    // Re = 1/nu,L=2pi, U = 1 und D = [0,2pi*L]^3
    const double U = 1 * uscaling;
    //const double L = 2 * M_PI;
    const double viscosity = 1.0 / Re;
    const double Ma = parser.getArgument<double>("Ma")*sqrt(1.4);
    const auto cfl = parser.getArgument<double>("CFL");
//    const double cs = U / Ma;

    // chose scaling so that the right Ma-number is achieved
    const double reference_temperature = parser.getArgument<double>("ref-temp");
    const double gamma = 1.4;
    const double scaling = sqrt(3) * U / (Ma*sqrt(gamma*reference_temperature));
//    const double scaling = sqrt(3) * cs; // choose different? -> stencil larger/smaller -> from turb. channel

    // setup configuration
    boost::shared_ptr<SolverConfiguration> configuration = boost::make_shared<SolverConfiguration>();
    if (restart > 0) configuration->setRestartAtIteration(restart);
    configuration->setUserInteraction(false);
    configuration->setOutputCheckpointInterval(parser.getArgument<int>("ncheckpoint"));
    configuration->setOutputSolutionInterval(nout);
    configuration->setNoOutputInterval(parser.getArgument<int>("n-no-out"));
    configuration->setSimulationEndTime(time);
    configuration->setOutputGlobalTurbulenceStatistics(true);
    configuration->setOutputCompressibleTurbulenceStatistics(true);
    configuration->setOutputShearLayerStatistics(true);
    configuration->setOutputShearLayerInterval(parser.getArgument<int>("nstats"));
    configuration->setMachNumber(Ma);
    configuration->setStencilScaling(scaling);
    configuration->setStencil(Stencil_D3Q45);
    configuration->setAdvectionScheme(SEMI_LAGRANGIAN);
    configuration->setEquilibriumScheme(QUARTIC_EQUILIBRIUM);
    configuration->setHeatCapacityRatioGamma(gamma);
    configuration->setReferenceTemperature(reference_temperature);
    configuration->setPrandtlNumber(0.71);
    configuration->setSedgOrderOfFiniteElement(parser.getArgument<int>("order")); // TODO: set to 4
    configuration->setCFL(cfl); // TODO: should be 0.4<CFL<2
    configuration->setServerEndTime(parser.getArgument<int>("server-end"));
//    configuration->setInitializationScheme(COMPRESSIBLE_ITERATIVE);

    parser.applyToSolverConfiguration(*configuration);

    // standard output dir
    string m_dirname; //configuration->getOutputDirectory();
    if (not parser.hasArgument("output-dir")){
        std::stringstream dirName;
        dirName << getenv("NATRIUM_HOME");
        dirName << "/step-mixingLayer/Re" << Re
                << "-Ma" << floor(Ma*1000)/1000
                << "-ref" << ref_level
                << "-p" << configuration->getSedgOrderOfFiniteElement()
                << "-mesh" << meshname
                << "-randu" << randuname << "x" << floor(randuscaling*1000)/1000
                << "-uscale" << uscaling
                << "-refT" << reference_temperature << "_" << bc;
//        dirName << "-coll" << static_cast<int>(configuration->getCollisionScheme())
//                << "-sl" << static_cast<int>(configuration->getAdvectionScheme())
        if (configuration->getAdvectionScheme() != SEMI_LAGRANGIAN)
            dirName << "-int" << static_cast<int>(configuration->getTimeIntegrator()) << "_" << static_cast<int>(configuration->getDealIntegrator());
        dirName << "-CFL" << configuration->getCFL();
//        dirName << "-sten" << static_cast<int>(configuration->getStencil());
//        if (configuration->isFiltering()) (dirName << "-filt" << static_cast<int>(configuration->getFilteringScheme()) << "by_max_degree");
        if (configuration->getRegularizationScheme() != NO_REGULARIZATION)
            dirName << "-reg" << static_cast<int>(configuration->getRegularizationScheme());
//        if (configuration->getEquilibriumScheme()!= BGK_EQUILIBRIUM) (dirName << "-equili" << static_cast<int>(configuration->getEquilibriumScheme()));
        if (configuration->getCollisionScheme() == MRT_STANDARD) {
            dirName << "-mrt" << static_cast<int>(configuration->getMRTBasis());
        }
        if (configuration->getCollisionScheme() == MRT_STANDARD) {
            dirName << "-relax" << static_cast<int>(configuration->getMRTRelaxationTimes());
        }
        m_dirname = dirName.str();
    } else {
        m_dirname = parser.getArgument<string>("output-dir");
    }
    configuration->setOutputDirectory(m_dirname);

    double deltaTheta0 = 0.093;
    // Grid resolution
    double len_x = parser.getArgument<double>("lx") * deltaTheta0;
    double len_y = parser.getArgument<double>("ly") * deltaTheta0;
    double len_z = parser.getArgument<double>("lz") * deltaTheta0;
    double center = parser.getArgument<double>("center");
    double dy_scaling = parser.getArgument<double>("dy-scaling");
    boost::shared_ptr<MixingLayer3D> mixingLayer = boost::make_shared<MixingLayer3D>
            (viscosity, ref_level, repetitions, randuscaling, randuname, len_x, len_y, len_z, meshname, center, dy_scaling, deltaTheta0, U * uscaling, reference_temperature, bc);
    MixingLayer3D::UnstructuredGridFunc trafo(len_y);

    double ymin = trafo.trans(len_y / repetitions.at(1) / pow(2, ref_level));
    const double p = configuration->getSedgOrderOfFiniteElement();
    const double dt = configuration->getCFL() / (p * p) / (sqrt(2) * scaling) * ymin;

    if (is_MPI_rank_0()) {
        LOG(WELCOME) << "MIXING LAYER SETUP: " << endl
                     << "===================================================" << endl
                     << "Dimensions:    " << len_x << " x " << len_y << " x " << len_z
                     << endl << "Grid:          " << repetitions.at(0) << " x "
                     << repetitions.at(1) << " x " << repetitions.at(2)
                     << " blocks with 8^" << ref_level << " cells each " << endl
                     << "#Cells:        " << int(repetitions.at(0) * pow(2, ref_level))
                     << " x " << int(repetitions.at(1) * pow(2, ref_level)) << " x "
                     << int(repetitions.at(2) * pow(2, ref_level)) << " = "
                     << int(
                             repetitions.at(0) * repetitions.at(1) * repetitions.at(2)
                             * pow(2, 3 * ref_level)) << endl << "#Points:       "
                     << int(repetitions.at(0) * pow(2, ref_level) * p) << " x "
                     << int(repetitions.at(1) * pow(2, ref_level) * p) << " x "
                     << int(repetitions.at(2) * pow(2, ref_level) * p) << " = "
                     << int(
                             repetitions.at(0) * repetitions.at(1) * repetitions.at(2)
                             * pow(2, 3 * ref_level) * p * p * p) << endl
                     << endl << "               dt  = " << dt << endl
                     << "===================================================" << endl
                     << endl;
    }
//    MixingLayer3D::UnstructuredGridFunc trafo(mixingLayer->lx, mixingLayer->lx, mixingLayer->lx);
    /////////////////////////////////////////////////
    // run solver
    //////////////////////////////////////////////////
    CompressibleCFDSolver<3> solver(configuration, mixingLayer);
    const size_t table_output_lines_per_10s = 300;
    configuration->setOutputTableInterval(1 + 10.0 / solver.getTimeStepSize() / table_output_lines_per_10s);
    solver.appendDataProcessor(boost::make_shared<ShearLayerStats>(solver, configuration->getOutputDirectory(), shearLayerThickness, Re));
    solver.run();
    pout << "step-mixingLayer terminated." << endl;
    return 0;
}
