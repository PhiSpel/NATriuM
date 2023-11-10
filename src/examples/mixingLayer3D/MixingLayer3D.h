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
        MixingLayer3D(double viscosity, size_t refinementLevel, double randu_scaling, string randuname,
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
        boost::shared_ptr<Mesh<3> > makeGrid(const string& meshname, double len_x, double len_y, double len_z, vector<unsigned int> repetitions = {1, 1, 1});

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
