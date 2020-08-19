add_executable(milk_carton_optimizer
               mod_particle.f90
               mod_potential.f90
               mod_particle_object.f90
               mod_carton.f90
               milk_carton_optimizer.f90)
install(TARGETS milk_carton_optimizer
        DESTINATION bin)
