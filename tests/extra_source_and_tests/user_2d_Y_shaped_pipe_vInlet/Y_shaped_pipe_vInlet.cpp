/**
 * @file 	T_shaped_pipe.cpp
 * @brief 	This is the benchmark test of multi-inlet and multi-outlet.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string wall_inner_path = "./input/Yshape_inner_boundary_for_code.dat";
std::string wall_outer_path = "./input/Yshape_outer_boundary_for_code.dat";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(Vec2d(-3.0, -12.0), Vecd(44.0, 16.0));
Real resolution_ref = 0.125;           /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
//Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */
//-------------------------------------------------------
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;    /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f * SMAX(Real(1), Real(7.0) / (Real(1.6) + Real(2.55)));
Real Re = 100.0;                    /**< Reynolds number. */
Real mu_f = rho0_f * U_f * Real(7.0) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygonFromFile(wall_inner_path, ShapeBooleanOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygonFromFile(wall_outer_path, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(wall_inner_path, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        //Real run_time = GlobalStaticVariables::physical_time_;
        //Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        //target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        target_velocity[0] = 1;
        target_velocity[1] = 0;

        return target_velocity;
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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);       // Tag for computation with save particles distribution
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment(); // handle command line arguments
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();

    // test buffer location
    /*Vec2d test_half = Vec2d(0.5 * 3.0, 0.5 * BW);
    Vec2d test_translation = Vec2d(32.716, 13.854);
    Real test_rotation = -0.8506;
    SolidBody test_up(
        sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(test_rotation), Vec2d(test_translation)), test_half, "TestBody"));
    test_up.defineParticlesAndMaterial<SolidParticles, Solid>();
    test_up.generateParticles<Lattice>();*/
    //----------------------------------------------------------------------
    //	SPH Particle relaxation section
    //----------------------------------------------------------------------
    /** check whether run particle relaxation for body fitted particle distribution. */
    if (sph_system.RunParticleRelaxation())
    {
        InnerRelation wall_inner(wall_boundary);
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_particles(wall_boundary);
        RelaxationStepInner relaxation_step_inner(wall_inner);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_wall_state_to_vtp({wall_boundary});
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        //----------------------------------------------------------------------
        //	Physics relaxation starts here.
        //----------------------------------------------------------------------
        random_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_wall_state_to_vtp.writeToFile(0.0);
        //----------------------------------------------------------------------
        // From here the time stepping begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            relaxation_step_inner.exec();
            ite++;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_wall_state_to_vtp.writeToFile(ite);
            }
        }

        std::cout << "The physics relaxation process of wall particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);

        return 0;
    }

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_wall_contact);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * 7.0);
    Vec2d emitter_translation = Vec2d(-0.5 * BW, 0.0);
    BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, inlet_particle_buffer);

    Vec2d inlet_flow_buffer_halfsize = emitter_halfsize;
    Vec2d inlet_flow_buffer_translation = emitter_translation;
    BodyAlignedBoxByCell inlet_flow_buffer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(inlet_flow_buffer_translation)), inlet_flow_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_condition(inlet_flow_buffer);
   
    Vec2d disposer_up_halfsize = Vec2d(0.5 * 3.0, 0.5 * BW);
    Vec2d disposer_up_translation = Vec2d(33.0, 14.2);
    Real disposer_up_rotation = -0.8506;
    BodyAlignedBoxByCell disposer_up(
        water_block, makeShared<AlignedBoxShape>(yAxis, Transform(Rotation2d(disposer_up_rotation), Vec2d(disposer_up_translation)), disposer_up_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_up_outflow_deletion(disposer_up);

    Vec2d disposer_down_halfsize = Vec2d( 0.5 * 4.0, 0.5 * BW);
    Vec2d disposer_down_translation = Vec2d(42.0, -9.7);
    Real disposer_down_rotation = 4.3807;
    BodyAlignedBoxByCell disposer_down(
        water_block, makeShared<AlignedBoxShape>(yAxis, Transform(Rotation2d(disposer_down_rotation), Vec2d(disposer_down_translation)), disposer_down_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_down_outflow_deletion(disposer_down);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    write_body_states.addToWrite<Real>(water_block, "Pressure");
    write_body_states.addToWrite<int>(water_block, "Indicator");
    write_body_states.addToWrite<Real>(water_block, "Density");

    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>>
        write_water_kinetic_energy(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 100.0;
    Real output_interval = end_time / 200.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                           /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";
                write_water_kinetic_energy.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** inflow injection*/
            emitter_inflow_injection.exec();
            disposer_up_outflow_deletion.exec();
            disposer_down_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        write_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_water_kinetic_energy.generateDataBase(1.0e-3);
    }
    else
    {
        write_water_kinetic_energy.testResult();
    }

    return 0;
}
