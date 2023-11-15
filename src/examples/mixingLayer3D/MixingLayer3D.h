/*
 * MixingLayer3D.h
 *
 *  Created on: Sep 18, 2014
 *      Author: dominik
 */

#ifndef MixingLayer3D_H_
#define MixingLayer3D_H_

/**
 * @file MixingLayer3D.h
 * @short Description of a simple Periodic Flow (in cubic domain).
 */

#include "deal.II/grid/tria.h"
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>
#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"
#include "deal.II/grid/grid_out.h"
#include <math.h>

#include <deal.II/lac/la_parallel_vector.h>
using namespace std;
using namespace natrium::DealIIExtensions;

namespace natrium {

/** @short Description of a simple Periodic Flow (flow in square domain).
 *  The domain is [0,1]^2. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 */
    class MixingLayer3D: public ProblemDescription<3> {
    public:
        /**
         * @short class to describe the x-component of the initial velocity
         * @note other are default (v0=w0=0, rho0=1)
         */
        class InitialVelocity: public dealii::Function<3> {
            private:
                MixingLayer3D* m_flow;
                vector<double> xvec;
                vector<double> yvec;
                vector<double> zvec;
//                DealIIExtensions::LinearAlgebra::distributed::Vector<
//                    DealIIExtensions::LinearAlgebra::distributed::Vector<double>> randomPsi;
                vector< vector< vector< vector<double> > > > randomPsi;
                vector< vector< vector< vector<double> > > > curlOfPsi;
                double minx, miny, minz;
                double maxx, maxy, maxz;
                double m_randu_scaling;
//                double dx, dy, dz;
//                double m_Ux;
                int nx, ny, nz;
                int kxmax, kymax, kzmax;
                bool m_print, m_recalculate;
                double InterpolateVelocities(double xq, double yq, double zq, const unsigned int dim) const;
            public:
                InitialVelocity(MixingLayer3D *flow, double randuscaling, string randuname);
                double value(const dealii::Point<3>& x, const unsigned int component = 0) const override;
        };
        class InitialDensity: public dealii::Function<3> {
            private:
                MixingLayer3D* m_flow;
            public:
                explicit InitialDensity(MixingLayer3D* flow) : m_flow(flow) { }
                virtual double value(const dealii::Point<3>& x, const unsigned int component = 0) const;
        };
        class InitialTemperature: public dealii::Function<3> {
        private:
            MixingLayer3D* m_flow;
        public:
            InitialTemperature(MixingLayer3D* flow) : m_flow(flow) { }
            virtual double value(const dealii::Point<3>& x, const unsigned int component = 0) const;
        };
        double m_initialT;

        /// constructor
        MixingLayer3D(double viscosity, size_t refinementLevel, vector<unsigned int> repetitions, double randu_scaling, string randuname,
                      double len_x, double len_y, double len_z, string meshname, double center, double scaling,
                      double deltaTheta0, double U = 1., double T = 1., string bc = "EQ_BC");
        /// destructor
        virtual ~MixingLayer3D();

        virtual void refine(Mesh<3>& mesh) {
            mesh.refine_global(m_refinementLevel);
        }
        virtual bool isCartesian() {return true;}

        /**
        * @short function to generate the unstructured mesh grid
        */
        struct UnstructuredGridFunc {
            double m_height;
            double m_center;
            double m_scaling;
            double m_a, m_b, m_c;
            UnstructuredGridFunc(double height, double center = 0.7, double scaling = 3) :
                    m_height(height / 2), m_center(center), m_scaling(scaling) {
                //// this one uses "m_scaling" to fix the maximum grid point distance (times the distance in the centre)
                m_a = (m_scaling - 1) / (2 - 2 * m_center);
                m_b = 1 - m_center * (m_scaling - 1) / (1 - m_center);
                m_c = (m_scaling - 1) / (2 - 2 * m_center) * m_center * m_center;
            }
            double trans(const double y) const {
                //// this one uses "m_scaling" to fix the maximum grid point distance (times the distance in the centre)
                double y_norm = y / m_height;
                if (std::abs(y_norm) < m_center) return y_norm * m_height;
                int sign = 1;
                if (y_norm < 0) {
                    sign = -1;
                    y_norm = std::abs(y_norm);
                }
                double y_new_norm = m_a * y_norm * y_norm + m_b * y_norm + m_c;
                return m_height * y_new_norm * sign;
            }
            dealii::Point<3> operator()(const dealii::Point<3> &in) const {
                return dealii::Point<3>(in(0), trans(in(1)), in(2));
            }
        };
        struct UnstructuredGridFunc2 {
            double m_height;
            double m_center = 0.7;
            double m_scaling = 1.5;
            double m_a, m_b, m_c;
            UnstructuredGridFunc2(double length, double height, double width, double gridDensity = 0.8) :
                    m_height(height / 2) {
                (void) length;
                (void) width;
                (void) gridDensity;
                //// this one uses "m_scaling" to fix the outer width of the domain
                m_a = (1 - m_scaling) / (-m_center*m_center + 2*m_center - 1);
                m_b = 1 - 2*m_a*m_center;
                m_c = m_scaling - m_a + 2*m_center*m_a - 1;
            }
            double trans(const double y) const {
                //// tangential transform
//                double y_norm = y / m_height;
//                double y_new_norm = std::tan(y_norm);
//                double y_new = m_height * y_new_norm; // std::tan(std::tan(std::tan(y_norm)));
//                if (is_MPI_rank_0()) cout << "y_norm = " << y << "/" << m_height << " = " << y_norm
//                                          << " -> y_new_norm = tan(tan(" << y_norm << ")) = " << y_new_norm
//                                          << " -> y_new = " << y_new_norm << "*" << m_height << " = " << y_new << endl;
//                return y_new;
                //// other trig. functions (asin, power)
//                double y_norm = y / m_height;
//                double y_new_norm = std::asin(y_norm); //; std::pow(y_norm, 5);
//                double y_new = m_height * y_new_norm;
//                return y_new;
                //// custom function (linear in 70% of domain, polynomail otherwise)
                //// this one uses "m_scaling" to fix the outer width of the domain
                double y_norm = y / m_height;
                if (std::abs(y_norm) < m_center) return y_norm * m_height;
                int sign = 1;
                if (y_norm < 0) {
                    sign = -1;
                    y_norm = std::abs(y_norm);
                }
                double y_new_norm = m_a*y_norm*y_norm + m_b*y_norm + m_c;
                return m_height * y_new_norm * sign;
            }
            dealii::Point<3> operator()(const dealii::Point<3> &in) const {
                return dealii::Point<3>(in(0), trans(in(1)), in(2));
            }
        };
        virtual void transform(Mesh<3>& mesh) {
            // transform grid to unstructured grid
            dealii::GridTools::transform(UnstructuredGridFunc(ly, m_center, m_scaling), mesh);
            // calculate ranges of dx, dy, dz
            vector<double> mindeltas(3,10), maxdeltas(3,0); // double mindx=0, maxdx=0, mindy=0, maxdy=0, mindz=0, maxdz=0;
            vector<double> deltas(3, 0); // double dx, dy, dz;
            vector<double> mincoords(3), maxcoords(3);
            //// get minimum and maximum coordinates
            for (typename Triangulation<3>::active_cell_iterator cell = mesh.begin_active(); cell != mesh.end(); ++cell) {
                for (size_t dim = 0; dim < 3; ++dim) {
                    mincoords.at(dim) = cell->vertex(0)[dim];
                    maxcoords.at(dim) = cell->vertex(0)[dim];
                }
                for (unsigned int f = 1; f < GeometryInfo<3>::vertices_per_cell; ++f) {
                    Point<3> x = cell->vertex(f);
                    for (size_t dim = 0; dim < 3; ++dim) {
                        mincoords.at(dim) = min(mincoords.at(0), x[dim]);
                        maxcoords.at(dim) = max(maxcoords.at(0), x[dim]);
                    }
                }
                for (size_t dim = 0; dim < 3; ++dim) {
                    deltas.at(dim) = maxcoords.at(dim) - mincoords.at(dim);
                    mindeltas.at(dim) = min(mindeltas.at(dim), deltas.at(dim));
                    maxdeltas.at(dim) = max(maxdeltas.at(dim), deltas.at(dim));
                }
            }

            // communicate
            for (size_t dim = 0; dim < 3; ++dim) {
                mindeltas.at(dim) = dealii::Utilities::MPI::min_max_avg(mindeltas.at(dim), MPI_COMM_WORLD).min;
                maxdeltas.at(dim) = dealii::Utilities::MPI::min_max_avg(maxdeltas.at(dim), MPI_COMM_WORLD).max;
            }//// calculate boundaries and set boundary ids
            if (is_MPI_rank_0()) LOG(DETAILED) << " dimensions: 3" << endl << " no. of cells: " << mesh.n_active_cells() << endl;
            double minx=0, maxx=0, miny=0, maxy=0, minz=0, maxz=0;
            //// get minimum and maximum coordinates
            for (typename Triangulation<3>::active_cell_iterator cell = mesh.begin_active(); cell != mesh.end(); ++cell) {
                for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f) {
                    if (cell->face(f)->at_boundary()) {
                        Point<3> x = cell->face(f)->center();
                        minx = min(minx, x[0]);
                        maxx = max(maxx, x[0]);
                        miny = min(miny, x[1]);
                        maxy = max(maxy, x[1]);
                        minz = min(minz, x[2]);
                        maxz = max(maxz, x[2]);
                    }
                }
            }
            // communicate
            minx = dealii::Utilities::MPI::min_max_avg(minx, MPI_COMM_WORLD).min;
            miny = dealii::Utilities::MPI::min_max_avg(miny, MPI_COMM_WORLD).min;
            minz = dealii::Utilities::MPI::min_max_avg(minz, MPI_COMM_WORLD).min;
            maxx = dealii::Utilities::MPI::min_max_avg(maxx, MPI_COMM_WORLD).max;
            maxy = dealii::Utilities::MPI::min_max_avg(maxy, MPI_COMM_WORLD).max;
            maxz = dealii::Utilities::MPI::min_max_avg(maxz, MPI_COMM_WORLD).max;

            lx = maxx-minx;
            ly = maxy-miny;
            lz = maxz-minz;

//            if (is_MPI_rank_0()) {
//                LOG(DETAILED) << "---------------------------------------" << endl
//                << "Mesh info before transform(): " << endl << "dx in [" << mindeltas.at(0) << "," << maxdeltas.at(0)
//                << "], dy in [" << mindeltas.at(1) << "," << maxdeltas.at(1)
//                << "], dz in [" << mindeltas.at(2) << "," << maxdeltas.at(2) << "]." << endl;
//                //// Print boundary indicators
//                map<types::boundary_id, unsigned int> boundary_count;
//                for (auto cell: mesh.active_cell_iterators()) {
//                    for (unsigned int face = 0; face < GeometryInfo<3>::faces_per_cell; ++face) {
//                        if (cell->face(face)->at_boundary()) boundary_count[cell->face(face)->boundary_id()]++;
//                    }
//                }
//                LOG(DETAILED) << " domain limits: x in[" << minx << "," << maxx
//                    << "], y in [" << miny << "," << maxy
//                    << "], z in [" << minz << "," << maxz << "]" << endl
//                    << " boundary indicators: ";
//                for (const pair<const types::boundary_id, unsigned int> &pair: boundary_count) {
//                    LOG(DETAILED) << pair.first << "(" << pair.second << " times) ";
//                }
//                LOG(DETAILED) << endl << "---------------------------------------" << endl;
//            }
	    }
        double lx, ly, lz, m_center, m_scaling, deltaTheta0;

    private:
        /// speed of sound
        double m_U;
        string m_bc;
        size_t m_refinementLevel;

        /**
         * @short create triangulation for couette flow
         * @return shared pointer to a triangulation instance
         */
        boost::shared_ptr<Mesh<3> > makeGrid(const string& meshname, double len_x, double len_y, double len_z, vector<unsigned int> repetitions);

        /**
         * @short create boundaries for couette flow
         * @return shared pointer to a vector of boundaries
         * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
         */
        boost::shared_ptr<BoundaryCollection<3> > makeBoundaries();
//    protected:
    };

} /* namespace natrium */

#endif /* MixingLayer3D_H_ */
