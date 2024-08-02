/**
 * @file     test_2d_test_v7.cpp
 * @brief 	test
 * @details  test
 *
 * @author 	Sukang Peng
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file
//----------------------------------------------------------------------
std::string fixed_left_bottom = "./input/fixed_left_bottom.dat";
std::string fixed_right_bottom = "./input/fixed_right_bottom.dat";
std::string fixed_left_top_down = "./input/fixed_left_top_down.dat";
std::string fixed_right_top_down = "./input/fixed_right_top_down.dat";
std::string fixed_left_top_up = "./input/fixed_left_top_up.dat";
std::string fixed_right_top_up = "./input/fixed_right_top_up.dat";
//----------------------------------------------------------------------
//	Basic geometry parameters
//----------------------------------------------------------------------
Real domain_bound_right = 3.5;   /**< Domain Bound length right part. */
Real domain_bound_left = 3.5;    /**< Domain Bound length left part. */
Real domain_bound_height = 8;    /**< Domain Bound height. */
Real domain_bound_height1 = 0.5; /**< Domain Bound height1. */

Real channel_width = 0.6; /**< Channel width. */

Real resolution_ref = 0.07;            /**< Reference resolution. */
Real BW = resolution_ref * 4.0;        /**< Boundary width, determined by specific layer of boundary particles. */
Real DL_sponge = resolution_ref * 4.0; /**< Sponge region to impose inflow condition. */
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;    /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f;
Real Re = 100.0;                               /**< Reynolds number. */
Real mu_f = rho0_f * U_f * channel_width / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real rho0_s = 1265.0;
Real poisson = 0.45;
Real Youngs_modulus = 5e4;
Real physical_viscosity = 200.0;
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class Fixed : public MultiPolygonShape
{
  public:
    explicit Fixed(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygonFromFile(fixed_left_bottom, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(fixed_right_bottom, ShapeBooleanOps::add);
		multi_polygon_.addAPolygonFromFile(fixed_left_top_down, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(fixed_right_top_down, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(fixed_left_top_up, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(fixed_right_top_up, ShapeBooleanOps::add);
    }
};

//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-domain_bound_left, -domain_bound_height1), Vec2d(domain_bound_right, domain_bound_height));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(true);  // tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(false);         // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments

    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody fixed(sph_system, makeShared<Fixed>("Fixed"));
    fixed.defineMaterial<Solid>();
    fixed.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_fixed_to_vtp(fixed);
    ReloadParticleIO write_fixed_reload_files(fixed);
    write_fixed_to_vtp.writeToFile(0);
    write_fixed_reload_files.writeToFile(0);
    return 0;
}
