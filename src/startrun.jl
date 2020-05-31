function runsimulation(input::JuSwalbe.Inputconstants)
    lx = input.lx
    ly = input.ly
    maxruntime = input.maxruntime
    
    if ly == 0
        println("Starting a one dimensional simulation")
        # Define an initial fluid state
        mom = JuSwalbe.Macroquant(height = 2 .+ sin.(4*π.*[1:lx;])./lx, 
                                  velocity = zeros(lx),
                                  pressure = zeros(lx),
                                  energy = zeros(lx))

        # Define an struct containing all forces
        forces = JuSwalbe.Forces(slip = zeros(lx), 
                                 h∇p = zeros(lx),
                                 bathymetry = zeros(lx),
                                 thermal = zeros(lx))

        # Set up an equilibrium distribtution from the chosen fluid state
        equilibrium = calc_equilibrium_distribution(mom)
        temp_dist = deepcopy(equilibrium)

        # Next we iterate until maxruntime is reached
        for time in 1:maxruntime
            println(mom.height[5], time)
            # Compute all forces
            computeslip(mom, forces, input)
            pressure(mom; γ=input.γ, θ=ones(Float64,1)*1/9)
            ∇p(mom, forces)
            computethermalcapillarywaves(mom, forces, input)

            # Compute the updated moments
            mom = calculatemoments(temp_dist, forces)

            # Add boundary conditions here

            # Performe the collision
            out_dist = collisionBGK(mom, forces, temp_dist, input)

            # Add obstacles and bounce back here

            # Stream the distribtutions
            temp_dist = streamdistperiodic(out_dist)
            
            if time%input.dumping == 0
                plot(mom.height, fill = (0, 0.5, :blue), label="h(x)")
                println("Hello $time")
                println(mom.height)
            end
        end
    else 
        println("Starting a two dimensional simulation")
    end

end
