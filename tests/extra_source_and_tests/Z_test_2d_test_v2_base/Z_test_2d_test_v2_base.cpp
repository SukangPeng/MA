/**
 * @file     test_2d_test_v2_base.cpp
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
std::string vessel_inner = "./input/vessel_inner.dat";
std::string vessel_outer = "./input/vessel_outer.dat";
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
Real rho0_s = 10.0; /**< Reference density.*/
Real poisson = 0.4; /**< Poisson ratio.*/
Real Ae = 1.4e3;    /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** create a fixed shape */
std::vector<Vecd> createFixedShape()
{
    // geometry
    std::vector<Vecd> fixed_shape;
    fixed_shape.push_back(Vecd(0.15, 1.5));
    fixed_shape.push_back(Vecd(-0.15, 1.5));
    fixed_shape.push_back(Vecd(-0.15, 1));
    fixed_shape.push_back(Vecd(0.15, 1));
    fixed_shape.push_back(Vecd(0.15, 1.5));

    return fixed_shape;
}
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &import_model_name) : MultiPolygonShape(import_model_name)
    {
        multi_polygon_.addAPolygonFromFile(vessel_inner, ShapeBooleanOps::add);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygonFromFile(vessel_outer, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(vessel_inner, ShapeBooleanOps::sub);
    }
};
class Fixed : public MultiPolygonShape
{
  public:
    explicit Fixed(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createFixedShape(), ShapeBooleanOps::add);
    }
};
/** create the emitter shape. */
MultiPolygon createFixedBaseShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(createFixedShape(), ShapeBooleanOps::add);
    return multi_polygon;
}
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
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        target_velocity[1] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[0] * position[0] / halfsize_[0] / halfsize_[0]);
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
    BoundingBox system_domain_bounds(Vec2d(-domain_bound_left, -domain_bound_height1), Vec2d(domain_bound_right, domain_bound_height));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false);  // tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);         // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments

    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineAdaptationRatios(1.15, 2.0);
    wall_boundary.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    wall_boundary.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation vessel_wall_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_vessel_wall_particles(wall_boundary);
        RelaxationStepInner relaxation_step_inner(vessel_wall_inner);
        BodyStatesRecordingToVtp write_vessel_wall_to_vtp(wall_boundary);
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_vessel_wall_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_vessel_wall_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the vessel wall N = " << ite_p << "\n";
                write_vessel_wall_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
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
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    /*
    TimeDependentAcceleration time_dependent_acceleration(Vec2d::Zero());
    SimpleDynamics<GravityForce> apply_gravity_force(water_block, time_dependent_acceleration);
    */

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_wall_contact);

    /** Calculate the convective time step of the fluid to ensure the numerical stability of the simulation.
    Calculate the acoustic time step of the fluid to ensure numerical stability of the simulation. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    /**Emitter & Disposer. */
    /*Emitter*/
    Vec2d emitter_halfsize = Vec2d(0.5 * channel_width, 0.5 * DL_sponge);
    Vec2d emitter_translation = Vec2d(-0.5 * channel_width, 0.0) + emitter_halfsize;
    BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(yAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, inlet_particle_buffer);
    /*Buffer*/
    Vec2d inlet_flow_buffer_halfsize = Vec2d(0.5 * channel_width, 0.5 * DL_sponge);
    Vec2d inlet_flow_buffer_translation = Vec2d(0.0, 0.5 * DL_sponge);
    BodyAlignedBoxByCell inlet_flow_buffer(water_block, makeShared<AlignedBoxShape>(yAxis, Transform(Vec2d(inlet_flow_buffer_translation)), inlet_flow_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_condition(inlet_flow_buffer);
    /*Disposer_right*/
    Vec2d disposer_right_halfsize = Vec2d(0.5 * channel_width, 0.5 * DL_sponge);
    Vec2d disposer_right_translation = Vec2d(3.159845357, 7.113012702);
    Real disposer_right_rotation = -1.0472;
    BodyAlignedBoxByCell disposer_right(
        water_block, makeShared<AlignedBoxShape>(yAxis, Transform(Rotation2d(disposer_right_rotation), Vec2d(disposer_right_translation)), disposer_right_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_right_outflow_deletion(disposer_right);
    /*Disposer_left*/
    Vec2d disposer_left_halfsize = Vec2d(0.5 * channel_width, 0.5 * DL_sponge);
    Vec2d disposer_left_translation = Vec2d(-3.159845357, 7.113012702);
    Real disposer_left_rotation = 1.0472;
    BodyAlignedBoxByCell disposer_left(
        water_block, makeShared<AlignedBoxShape>(yAxis, Transform(Rotation2d(disposer_left_rotation), Vec2d(disposer_left_translation)), disposer_left_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_left_outflow_deletion(disposer_left);

    /** Test Shape
        // water block
    BodyStatesRecordingToVtp write_water_to_vtp(water_block);
    ReloadParticleIO write_water_reload_files(water_block);
    write_water_to_vtp.writeToFile(0);
    write_water_reload_files.writeToFile(0); //

    // WallBoundary
    BodyStatesRecordingToVtp write_wall_to_vtp(wall_boundary);
    ReloadParticleIO write_wall_reload_files(wall_boundary);
    write_wall_to_vtp.writeToFile(0);
    write_wall_reload_files.writeToFile(0); //

   // emitter
    SolidBody emitter_(sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize), "Emitter");
    emitter_.defineMaterial<Solid>();
    emitter_.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_emitter_to_vtp(emitter_);
    ReloadParticleIO write_emitter_reload_files(emitter_);

    write_emitter_to_vtp.writeToFile(0);
    write_emitter_reload_files.writeToFile(0);

    // buffer
    SolidBody buffer(sph_system, makeShared<AlignedBoxShape>(yAxis, Transform(Vec2d(inlet_flow_buffer_translation)), inlet_flow_buffer_halfsize), "Buffer");
    buffer.defineMaterial<Solid>();
    buffer.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_buffer_to_vtp(buffer);
    ReloadParticleIO write_buffer_reload_files(buffer);

    write_buffer_to_vtp.writeToFile(0);
    write_buffer_reload_files.writeToFile(0);

    // disposer_right
    SolidBody disposer_right_(sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(disposer_right_rotation), Vec2d(disposer_right_translation)), disposer_right_halfsize), "Disposer_right");
    disposer_right_.defineMaterial<Solid>();
    disposer_right_.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_disposer_right_to_vtp(disposer_right_);
    ReloadParticleIO write_disposer_right_reload_files(disposer_right_);

    write_disposer_right_to_vtp.writeToFile(0);
    write_disposer_right_reload_files.writeToFile(0);

    // disposer_left
    SolidBody disposer_left_(sph_system, makeShared<AlignedBoxShape>(xAxis, Transform(Rotation2d(disposer_left_rotation), Vec2d(disposer_left_translation)), disposer_left_halfsize), "Disposer_left");
    disposer_left_.defineMaterial<Solid>();
    disposer_left_.generateParticles<BaseParticles, Lattice>();

    BodyStatesRecordingToVtp write_disposer_left_to_vtp(disposer_left_);
    ReloadParticleIO write_disposer_left_reload_files(disposer_left_);

    write_disposer_left_to_vtp.writeToFile(0);
    write_disposer_left_reload_files.writeToFile(0);
    */

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    write_body_states.addToWrite<Real>(water_block, "Pressure"); // output for debug
    write_body_states.addToWrite<int>(water_block, "Indicator"); // output for debug
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_water_kinetic_energy(water_block);
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
            disposer_right_outflow_deletion.exec();
            disposer_left_outflow_deletion.exec();

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
