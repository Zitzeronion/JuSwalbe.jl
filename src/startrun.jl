using Plots

function runsimulation_1D(input::JuSwalbe.Inputconstants)
    lx = input.lx
    ly = input.ly
    maxruntime = input.maxruntime

    h0 = 2 .+ sin.(4*π.*[1:lx;]./lx)
    if ly == 0
        println("Starting a one dimensional simulation")
        # Define an initial fluid state
        mom = JuSwalbe.Macroquant(height = h0, 
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
            # Compute the updated moments
            mom = calculatemoments(temp_dist)

            # Compute all forces
            computeslip(mom, forces, input)
            pressure(mom, γ=input.γ, θ=ones(Float64,1)*1/9)
            ∇p(mom, forces)
            computethermalcapillarywaves(mom, forces, input)


            # Add boundary conditions here

            # Performe the collision
            out_dist = collisionBGK(mom, forces, temp_dist, input)

            # Add obstacles and bounce back here

            # Stream the distribtutions
            temp_dist = streamdistperiodic(out_dist)
            
            if time%100 == 0
                println("Heigth: ", mom.height[5], " Velocity: ", mom.velocity[5], " Pressure: ", mom.pressure[5], " Time step: ", time)
                println("Forces -> Pressure-grad: ", forces.h∇p[5], " slippage: ", forces.slip[5])
            end
        end
        plot(h0)
        plot!(mom.height)
        gui()
        return mom
    else 
        println("Please use `runsimulation_2D`")
    end
end

function runsimulation_2D(input::JuSwalbe.Inputconstants)
    lx = input.lx
    ly = input.ly
    maxruntime = input.maxruntime

    h0 = [2 + (sin(4*π*i/lx) * sin(4*π*j/ly)) for i in 1:lx, j in 1:ly] 
    if ly == 0
        println("Please use `runsimulation_1D`")
    else 
        println("Starting Simulation in two dimensions")
        # Define an initial fluid state
        mom = JuSwalbe.Macroquant(height = h0, 
                                  velocity = zeroTwovector(lx, ly),
                                  pressure = zeros(lx, ly),
                                  energy = zeros(lx, ly))

        # Define an struct containing all forces
        forces = JuSwalbe.Forces(slip = zeroTwovector(lx, ly), 
                                 h∇p = zeroTwovector(lx, ly),
                                 bathymetry = zeroTwovector(lx, ly),
                                 thermal = zeroTwovector(lx, ly))

        # Set up an equilibrium distribtution from the chosen fluid state
        equilibrium = calc_equilibrium_distribution(mom)
        temp_dist = deepcopy(equilibrium)

        # Next we iterate until maxruntime is reached
        for time in 1:maxruntime
            # Compute the updated moments
            mom = @spawn calculatemoments(temp_dist)

            # Compute all forces
            computeslip(fetch(mom), forces, input)
            pressure(fetch(mom), γ=input.γ, θ=ones(Float64,(1,1))*1/9)
            ∇p(fetch(mom), forces)
            computethermalcapillarywaves(fetch(mom), forces, input)


            # Add boundary conditions here

            # Performe the collision
            out_dist = @spawn collisionBGK(fetch(mom), forces, temp_dist, input)

            # Add obstacles and bounce back here

            # Stream the distribtutions
            temp_dist = streamdistperiodic(fetch(out_dist))
            if time % 1000 == 0
                println("Time step: $time")
            end

        end
        s1 = heatmap(h0)
        s2 = heatmap(fetch(mom).height)
        plot(s1, s2)
        gui()
        return mom
    end

end

function runsimulation_2DGPU(input::JuSwalbe.Inputconstants)
    lx = input.lx
    ly = input.ly
    maxruntime = input.maxruntime

    h0 = [2 + (sin(4*π*i/lx) * sin(4*π*j/ly)) for i in 1:lx, j in 1:ly]
    h0gpu = Float32.(h0) 
    hnew = CUDA.CuArray(h0gpu)
    if ly == 0
        println("Please use `runsimulation_1D`")
    else 
        println("Starting Simulation in two dimensions")
        # Define an initial fluid state
        mom = JuSwalbe.Macroquant(height = hnew, 
                                  velocity = zeroTwovectorCU(lx, ly),
                                  pressure = CUDA.zeros(lx, ly),
                                  energy = CUDA.zeros(lx, ly))

        # Define an struct containing all forces
        forces = JuSwalbe.Forces(slip = zeroTwovectorCU(lx, ly), 
                                 h∇p = zeroTwovectorCU(lx, ly),
                                 bathymetry = zeroTwovectorCU(lx, ly),
                                 thermal = zeroTwovectorCU(lx, ly))

        # Set up an equilibrium distribtution from the chosen fluid state
        equilibrium = calc_equilibrium_distribution(mom)
        temp_dist = deepcopy(equilibrium)

        # Next we iterate until maxruntime is reached
        for time in 1:maxruntime
            # Compute the updated moments
            mom = calculatemoments(temp_dist)

            # Compute all forces
            computeslip(mom, forces, input)
            pressure(mom, γ=input.γ, θ=ones(Float32,(1,1))*1/9)
            ∇p(mom, forces)
            computethermalcapillarywaves(mom, forces, input)


            # Add boundary conditions here

            # Performe the collision
            out_dist = collisionBGK(mom, forces, temp_dist, input)

            # Add obstacles and bounce back here

            # Stream the distribtutions
            temp_dist = streamdistperiodic(out_dist)
            if time % 1000 == 0
                println("Time step: $time")
            end

        end
        s1 = heatmap(h0)
        s2 = heatmap(mom.height)
        plot(s1, s2)
        gui()
        return mom
    end

end