/**
* @file     test_2d_test_v4.cpp
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
Real domain_bound_height = 8;  /**< Domain Bound height. */
Real domain_bound_height1 = 0.5; /**< Domain Bound height1. */

Real channel_width = 0.6; /**< Channel width. */

Real resolution_ref = 0.07;           /**< Reference resolution. */
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

class VesselWall : public MultiPolygonShape
{
  public:
    explicit VesselWall(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygonFromFile(vessel_outer, ShapeBooleanOps::add);
        multi_polygon_.addAPolygonFromFile(vessel_inner, ShapeBooleanOps::sub);
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
    sph_system.setRunParticleRelaxation(false);   // tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(true);        // Tag for computation with save particles distribution
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments

    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);

    SolidBody vessel_wall(sph_system, makeShared<VesselWall>("VesselWall"));
    vessel_wall.defineAdaptationRatios(1.15, 2.0);
    vessel_wall.defineBodyLevelSetShape()->writeLevelSet(sph_system);
    vessel_wall.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? vessel_wall.generateParticles<BaseParticles, Reload>(vessel_wall.getName())
        : vessel_wall.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation vessel_wall_inner(vessel_wall);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_vessel_wall_particles(vessel_wall);
        RelaxationStepInner relaxation_step_inner(vessel_wall_inner);
        BodyStatesRecordingToVtp write_vessel_wall_to_vtp(vessel_wall);
        ReloadParticleIO write_particle_reload_files(vessel_wall);
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
    InnerRelation vessel_wall_inner(vessel_wall);
    ContactRelation water_block_contact(water_block, {&vessel_wall, &wall_boundary});
    ContactRelation vessel_wall_contact(vessel_wall, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> vessel_wall_normal_direction(vessel_wall);
  
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_block_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> vessel_wall_corrected_configuration(vessel_wall_inner);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> vessel_wall_stress_relaxation_first_half(vessel_wall_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> vessel_wall_stress_relaxation_second_half(vessel_wall_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> vessel_wall_computing_time_step_size(vessel_wall);

    //Constrain
    /*   
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> insert_body_computing_time_step_size(insert_body);
    BodyRegionByParticle beam_base(insert_body, makeShared<MultiPolygonShape>(createBeamBaseShape()));
    SimpleDynamics<FixBodyPartConstraint> constraint_beam_base(beam_base);
    */

    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_block_contact);

    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_block_contact);

    /** Calculate the convective time step of the fluid to ensure the numerical stability of the simulation.
    Calculate the acoustic time step of the fluid to ensure numerical stability of the simulation. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    /*Emitter & Disposer. */
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

    /* compute vorticity in the fluid field
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
    */

    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(vessel_wall);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> vessel_wall_update_normal(vessel_wall);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(vessel_wall_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(vessel_wall_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_water_kinetic_energy(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing surface normal direction for the wall. */
    wall_boundary_normal_direction.exec();
    /** computing surface normal direction for the vessel wall. */
    vessel_wall_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    vessel_wall_corrected_configuration.exec();
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
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
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

            /** FSI for viscous force. */
            viscous_force_from_fluid.exec();
            /** Update normal direction on elastic body.*/
            vessel_wall_update_normal.exec();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                inflow_condition.exec();
                /** FSI for pressure force. */
                pressure_force_from_fluid.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);
           
                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(vessel_wall_computing_time_step_size.exec(), dt - dt_s_sum);
                    vessel_wall_stress_relaxation_first_half.exec(dt_s);
                    vessel_wall_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
                write_water_kinetic_energy.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
            /** inflow injection*/
            emitter_inflow_injection.exec();
            disposer_left_outflow_deletion.exec();
            disposer_right_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            water_block_complex.updateConfiguration();

            /** one need update configuration after periodic condition. */
            vessel_wall.updateCellLinkedList();
            vessel_wall_contact.updateConfiguration();
        }
        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

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
