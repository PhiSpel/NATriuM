/**
 * @file DiamondObstacle2D.h
 * @short Flow around an obstacle
 * @date 14.12.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef DIAMONDOBSTACLE2D_H_
#define DIAMONDOBSTACLE2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

/**
 * @short Description of a flow around a diamond-shaped obstacle
 */
class DiamondObstacle2D: public ProblemDescription<2> {

private:
	double m_meanInflowVelocity;
	size_t m_refinementLevel;

public:
	class InitialVelocity: public dealii::Function<2> {
	private:
		DiamondObstacle2D* m_flow;
	public:
		InitialVelocity(DiamondObstacle2D* flow) : m_flow(flow) { }
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	class InitialDensity: public dealii::Function<2> {
	private:
		DiamondObstacle2D* m_flow;
	public:
		InitialDensity(DiamondObstacle2D* flow) : m_flow(flow) { }
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	class InitialTemperature: public dealii::Function<2> {
	private:
		DiamondObstacle2D* m_flow;
	public:
		InitialTemperature(DiamondObstacle2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	class InflowVelocity: public dealii::Function<2> {
	private:
		double m_averageU;
	public:
		InflowVelocity(double av_U) : m_averageU(av_U) {}
		virtual double value(const dealii::Point<2> &x, const unsigned int component) const {
//			double H = 4.1;
            (void) x;
			if (component == 0) {
				//return 4 * m_averageU * x(1) * (H - x(1)) / (H * H);
				return m_averageU;
			}
			return 0.0;
		}
		virtual void vector_value(const dealii::Point<2> &x, dealii::Vector<double> &return_value) const {
			return_value(0) = value(x, 0);
			return_value(1) = 0.0;
		}
	};
	/// constructor
	DiamondObstacle2D(double velocity, double viscosity, size_t refinementLevel, int aoa);

	/// destructor
	virtual ~DiamondObstacle2D();
	virtual double getCharacteristicVelocity() const {
		return m_meanInflowVelocity;
	}
	virtual void refine(Mesh<2>& mesh) {
		mesh.refine_global(m_refinementLevel);
	}
	virtual void transform(Mesh<2>& ){}
	virtual bool isCartesian(){
		return false;
	}

private:
	boost::shared_ptr<Mesh<2> > makeGrid(size_t refinementLevel, int aoa);
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();
};

} /* namespace natrium */
#endif /* LIDDRIVENCAVIT2D_H_ */
