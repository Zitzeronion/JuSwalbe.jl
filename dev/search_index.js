var documenterSearchIndex = {"docs":
[{"location":"manual/#Manual-1","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"manual/#","page":"Manual","title":"Manual","text":"Short and clean instructions how to use the Package. Functions and methods can be found here with usage examples.","category":"page"},{"location":"manual/#","page":"Manual","title":"Manual","text":"Modules = [JuSwalbe]\nOrder = [:function, :type, :struct] ","category":"page"},{"location":"manual/#JuSwalbe.calc_equilibrium_distribution-Union{Tuple{JuSwalbe.Macroquant{Array{T,1},Array{T,1}}}, Tuple{T}} where T<:Number","page":"Manual","title":"JuSwalbe.calc_equilibrium_distribution","text":"calc_equilibrium_distribution(mom::JuSwalbe.Macroquant; gravity=0)\n\nCalculates the equilibrium distributions based on the macroscopic quantities mom and gravitational acceleration gravity.\n\nThe equilibrium distribtutions are at the heart of the lattice Boltzmann method. As the expansion is made around the equilibrium, lattice Boltzmann always assumes that the flow field is close to equilibrium. Therefore the equilibrium is calculated from macroscopic quantities, i.e. height .height and velocity .velocity.\n\nMath:\n\nOne spatial dimension\n\nDefining equations in one spatial dimensions are expressed in the following way\n\nf_0^eq =  h - frac12v^2gh^2 - frac1v^2hu^2 \n\nf_1^eq =  frac14v^2gh^2 + frac12vhu + frac12v^2hu^2\n\nf_2^eq =  frac14v^2gh^2 - frac12vhu + frac12v^2hu^2\n\nTwo spatial dimension\n\nIn two spatial dimensions the velocity because a vector with a x and y component. Such the equations look a little more complex\n\nf_0^eq =  h - frac49h(frac152gh - frac32u^2)\n\nf_i^eq = w_i h(frac32gh + 3 mathbfc_icdotmathbfu + frac92(mathbfc_icdotmathbfu)^2 - frac32u^2))\n\nwith i are the number of lattice speeds and wi and ci are the weights and sets of lattice velocities. \n\nExample\n\nTBD\n\nReferences\n\nOne spatial dimension\n\nThere are plenty of ways to calculate them I mainly the ones derived in Eq.(13) of Study of the 1D lattice Boltzmann shallow water equation and its coupling to build a canal network.\n\nTwo spatial dimensions\n\nHere I stick to the paper written by Paul Dellar, Eq.(26) of Nonhydrodynamic modes and a priori construction of shallow water lattice Boltzmann equations Originally these equilibria have been worked out by Paul Salmon.\n\n\n\n\n\n","category":"method"},{"location":"manual/#JuSwalbe.pressure-Union{Tuple{JuSwalbe.Macroquant{Array{T,2},JuSwalbe.Twovector{Array{T,2}}}}, Tuple{T}} where T<:Number","page":"Manual","title":"JuSwalbe.pressure","text":"pressure(mom::JuSwalbe.Macroquant; γ = 0.01, θ = 1/9)\n\nFilm pressure of the thin film equation.\n\nThis is just the summation of the laplacian of the height Δh and the disjoining potential Π. Since this function uses Π one has to be careful of θ as it should be an array. However it can be either a single element array or as large as the whole domain. See also: Δh, Π\n\nMath\n\nThe film pressure is simply\n\np_film = -gammaDelta h + Pi(h)\n\nThe sign is not 100% fixed and can change from paper to paper. Here I choose to be in agreement with Thiele et al., which is known to be a good theoretician. \n\nExample\n\njulia> using JuSwalbe\n\njulia> height = reshape([i for i in 1.0:1.0:16.0],4,4)\n4×4 Array{Float64,2}:\n 1.0  5.0   9.0  13.0\n 2.0  6.0  10.0  14.0\n 3.0  7.0  11.0  15.0\n 4.0  8.0  12.0  16.0\n\njulia> p = pressure(height)\n4×4 Array{Float64,2}:\n -0.200016  -0.0400001   -0.04        0.12\n -0.160002  -7.44536e-8  -1.6082e-8   0.16\n -0.160001  -4.68862e-8  -1.20826e-8  0.16\n -0.12       0.04         0.04        0.2\n\n\nSo what does this mean? A positive pressure is force that drives the film down in height. While a negative pressure generates a flux towards that location. In the example a linear increasing height field was used, an equilibrium though is reached when the film is flat.\n\nReferences\n\nRecent (short)\n\nSignatures of slip in dewetting polymer films\n\nReview\n\nDynamics and stability of thin liquid films\n\n\n\n\n\n","category":"method"},{"location":"manual/#JuSwalbe.readinput-Tuple{Any}","page":"Manual","title":"JuSwalbe.readinput","text":"readinput(file)\n\nReads input parameters from a file.\n\nThe expected amount of parameters can be addressed with inputconstants. For now it expects seven values for different runtime constants.\n\nExample\n\njulia> using JuSwalbe, DelimitedFiles\n\njulia> args = [\"Lattice_points_x\" 10; \"Lattice_points_y\" 5; \"Max_run_time\" 1000; \"Output_dump\" 100; \"gravity\" 0.0; \"surface_tension\" 0.01; \"slippage\" 1.0] # Generate a text file with input\n7×2 Array{Any,2}:\n \"Lattice_points_x\"    10\n \"Lattice_points_y\"     5\n \"Max_run_time\"      1000\n \"Output_dump\"        100\n \"gravity\"              0.0\n \"surface_tension\"      0.01\n \"slippage\"             1.0\n\njulia> writedlm(\"test.txt\", args)\n\njulia> test = readinput(\"test.txt\")\nJuSwalbe.inputconstants(10, 5, 1000, 100, 0.0, 0.01, 1.0)\n\njulia> test.lx\n10\n\njulia> test.γ\n0.01\n\njulia> test.γ + test.δ\n1.01\n\njulia> isa(test.lx + test.gravity, Int32)\nfalse\n\njulia> rm(\"test.txt\")\n\n\n\n\n\n","category":"method"},{"location":"manual/#JuSwalbe.Δh-Union{Tuple{JuSwalbe.Macroquant{Array{T,2},JuSwalbe.Twovector{Array{T,2}}}}, Tuple{T}} where T<:Number","page":"Manual","title":"JuSwalbe.Δh","text":"Δh(mom::JuSwalbe.Macroquant)\n\nCalculates the laplacian of the height field with periodic boundaries.\n\nThe calculation of the laplacian is central for the thin film evolution. Only the pressure gradient will induce a flow, at least for the non-fluctuating version. Therefore it is fairly important to have an accurate computation of the laplacian.\n\nMath\n\nThe Laplace equation is given by \n\nDelta rho = 0\n\nThe laplace operator is simply Δ = ∇ . ∇.  Since the film pressure has a contribution Δh we need a discreticed laplace operator.\n\nTwo spatial dimensions\n\nIn two dimensions I follow the paper from Santosh and Succi. A nine-point stencil is used and the equation goes as follow\n\nDelta h_ij = frac164sum_nnh_ij + sum_diagh_ij - 20 h_ij\n\nWhere i and j are the x and y coordinates. The index nn relates to the nearest neighbors, all four elements which are exactly Δx away from (i,j). On the other hand the diagonal elements are those four which are have a distance √2Δx to (i,j).  Periodicity is taken care of by a circular padding of the height array. \n\nOne spatial dimension\n\nFor a single spatial dimension it is mostly sufficient to use the central difference approach. This approach is given by\n\n`` \\partialx^2 h{i} = h{i-1} - 2h{i} + h_{i+1}\n\nwith i being the index along the spatial dimension.\n\nExample\n\njulia> using JuSwalbe\n\njulia> height = reshape(collect(1.0:16.0), (4,4))\n4×4 Array{Float64,2}:\n 1.0  5.0   9.0  13.0\n 2.0  6.0  10.0  14.0\n 3.0  7.0  11.0  15.0\n 4.0  8.0  12.0  16.0\n\njulia> velx = zeros(Float64, (4,4))\n4×4 Array{Float64,2}:\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n\njulia> vely = zeros(Float64, (4,4))\n4×4 Array{Float64,2}:\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n \njulia> vel = JuSwalbe.Twovector(velx, vely)\nJuSwalbe.Twovector{Array{Float64,2}}([0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])\n\njulia> moment = JuSwalbe.Macroquant(height, vel, zeros(4,4), zeros(4,4))\nJuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}([1.0 5.0 9.0 13.0; 2.0 6.0 10.0 14.0; 3.0 7.0 11.0 15.0; 4.0 8.0 12.0 16.0], JuSwalbe.Twovector{Array{Float64,2}}([0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]), [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])\n\njulia> moment.height\n4×4 Array{Float64,2}:\n 1.0  5.0   9.0  13.0\n 2.0  6.0  10.0  14.0\n 3.0  7.0  11.0  15.0\n 4.0  8.0  12.0  16.0\n\njulia> laplace = Δh(moment)\n4×4 Array{Float64,2}:\n 20.0   4.0   4.0  -12.0\n 16.0   0.0   0.0  -16.0\n 16.0   0.0   0.0  -16.0\n 12.0  -4.0  -4.0  -20.0\n\njulia> height = collect(1.0:16.0)\n16-element Array{Float64,1}:\n  1.0\n  2.0\n  3.0\n  4.0\n  5.0\n  6.0\n  7.0\n  8.0\n  9.0\n 10.0\n 11.0\n 12.0\n 13.0\n 14.0\n 15.0\n 16.0\n\njulia> Δh(height)\n16-element Array{Float64,1}:\n  16.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n   0.0\n -16.0\n\n\nReferences\n\nTwo spatial dimensions\n\nThere are plenty of papers concerning discret differentail operators. Many of them are good, although definitly not the best I go here with:\n\nIsotropic discrete Laplacian operators from lattice hydrodynamics\n\nOne spatial dimension\n\nAlmost any reference on discrete differentiation is good enough.\n\n\n\n\n\n","category":"method"},{"location":"manual/#JuSwalbe.Π-Union{Tuple{JuSwalbe.Macroquant{Array{T,2},JuSwalbe.Twovector{Array{T,2}}}}, Tuple{T}} where T<:Number","page":"Manual","title":"JuSwalbe.Π","text":"Π(mom::JuSwalbe.Macroquant; h_star = 0.1, exponents = [9,3], γ = 0.01, θ = 1.0/9.0)\n\nCalculates the disjoining pressure for a given surface tension and contact angle at every lattice point.\n\nThe disjoing pressure potential enables the simulation to study phenomena like dewetting of thin liquid films. Only with this potential it is possible to \"dewett\" without using highly sophisticated boundary conditions at hrightarrow 0. The functional shape of this potential can vary, ranging from the powerlaw used here to exponentials to a mixture of both. Using the powerlaw shape the exponents (9,3) mimic the Lenard-Jones potential for the liquid substrate interaction.\n\nArguments\n\nmom::JuSwalbe.Macroquant : Macroscopic moments of the simulation, important here .height which contains the height field. Allows for array input as well.\nh_star::T : Minimum of the pontential. T needs to be a subtype of Number! \nexponents::Array{Int64,1} : Exponents for the powerlaw potential\nγ::T : Surface tension\nθ::T : Equilibrium contact angle in multiple of π (for accuarcy reasons). Needs to be an array when doing patterning. For simple substrate use a one entry array with correct dimension.\n\nMath\n\nThe potential can be derived from the assumption that a given surface energy demands an equilibrium contact angle for the fluid. For the exact derivation take a look at the references, in principle it is the derivative of the interfacial potential with respect to h.\n\nPhi(h) = Pi(h) \n\nPi(h) = kappa f(h) = (1-cos(theta))frac(n-1)(m-1)(n-m)h_astBiggBigg(\frach_asthBigg)^n - Bigg(\frach_asthBigg)^mBigg\n\nPi(h_ast) = 0 \n\nExample\n\njulia> using JuSwalbe\n\njulia> n,m = (4,4)\n(4, 4)\n\njulia> moment = JuSwalbe.Macroquant{Matrix{Float64}, JuSwalbe.Twovector{Matrix{Float64}}}(ones(n,m), JuSwalbe.Twovector(zeros(n,m),zeros(n,m)), zeros(n,m),zeros(n,m))\nJuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}([1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0], JuSwalbe.Twovector{Array{Float64,2}}([0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]), [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])\n\njulia> p = Π(moment)\n4×4 Array{Float64,2}:\n -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5\n -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5\n -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5\n -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5\n\njulia> p = Π(moment, θ=zeros(1,1)) # Fully wetting substrate\n4×4 Array{Float64,2}:\n -0.0  -0.0  -0.0  -0.0\n -0.0  -0.0  -0.0  -0.0\n -0.0  -0.0  -0.0  -0.0\n -0.0  -0.0  -0.0  -0.0\n\njulia> h = reshape([i for i in 1.0:1.0:(n*m)],n,m)\n4×4 Array{Float64,2}:\n 1.0  5.0   9.0  13.0\n 2.0  6.0  10.0  14.0\n 3.0  7.0  11.0  15.0\n 4.0  8.0  12.0  16.0\n\njulia> p = Π(h)\n4×4 Array{Float64,2}:\n -1.6082e-5   -1.28656e-7  -2.20603e-8  -7.31997e-9\n -2.01025e-6  -7.44536e-8  -1.6082e-8   -5.86078e-9\n -5.95628e-7  -4.68862e-8  -1.20826e-8  -4.76503e-9\n -2.51281e-7  -3.14101e-8  -9.30669e-9  -3.92626e-9\n\n\nReferences\n\nTo get a good understanding:\n\nLong-scale evolution of thin liquid films\nWetting and spreading\nDynamics and stability of thin liquid films\n\nA rather recent new setup for the shape of Π can be found in \n\nSignatures of slip in dewetting polymer films \n\n\n\n\n\n","category":"method"},{"location":"manual/#JuSwalbe.velocitysquared-Union{Tuple{JuSwalbe.Macroquant{Array{T,2},JuSwalbe.Twovector{Array{T,2}}}}, Tuple{T}} where T<:Number","page":"Manual","title":"JuSwalbe.velocitysquared","text":"velocitysquared(mom::Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}})\n\nComputes the square of the velocity vector (ux, uy) at every lattice point.\n\nThe magnitude of the velocity is needed to calculate the equilibrium distribution in the dimensional case. In the one dimensional case the velocity is just a vector and therefore has no .x and .y component. See also calc_equilibrium_distribution\n\nMath\n\nThe velocity squared u^2(xy) is computed according to \n\nu^2(xy) = (u_x u_y)^2(xy) = u_x^2(xy) + u_y^2(xy)\n\nWith lower case x and y the respective component of the velocity vector is addressed.\n\nExample\n\njulia> using JuSwalbe\n\njulia> velocities = JuSwalbe.Twovector{Matrix{Float64}}(fill(0.1, (4,4)),fill(0.2, (4,4)))\nJuSwalbe.Twovector{Array{Float64,2}}([0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1], [0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2])\n\njulia> moment = JuSwalbe.Macroquant{Matrix{Float64},JuSwalbe.Twovector{Matrix{Float64}}}(ones(Float64, (4,4)), velocities, ones(Float64, (4,4)), ones(Float64, (4,4)))\nJuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}([1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0], JuSwalbe.Twovector{Array{Float64,2}}([0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1], [0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2]), [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0], [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0])\n\njulia> moment.velocity.x\n4×4 Array{Float64,2}:\n 0.1  0.1  0.1  0.1\n 0.1  0.1  0.1  0.1\n 0.1  0.1  0.1  0.1\n 0.1  0.1  0.1  0.1\n\njulia> JuSwalbe.velocitysquared(moment)\n4×4 Array{Float64,2}:\n 0.05  0.05  0.05  0.05\n 0.05  0.05  0.05  0.05\n 0.05  0.05  0.05  0.05\n 0.05  0.05  0.05  0.05\n\n\n\n\n\n","category":"method"},{"location":"manual/#JuSwalbe.inputconstants","page":"Manual","title":"JuSwalbe.inputconstants","text":"inputconstants = new(lx, ly, maxruntime, dumping, gravity, γ, δ)\n\nStruct containing input parameters.\n\nContains .lx lattice points in x-direction, .ly lattice points in y-direction.  Other fields are .maxruntime for the maximal number of time steps and .dumping to limit the number of output files. On top of these there are physical quantities such as .gravity, .γ and .δ  for the values of gravitational acceleration, fluids surface tension and the slip length. The example relates to an quadratic lattice 20 times 20 lattice units in area.  Run for 100 lattice Boltzmann time steps only printing output every 10 time steps. Having no gravity and a surface tension of 0.01 and a slip length of 1.  \n\nExample\n\njulia> using JuSwalbe\n\njulia> new_input = JuSwalbe.inputconstants(20, 20, 100, 10, 0.0, 0.01, 1.0)\nJuSwalbe.inputconstants(20, 20, 100, 10, 0.0, 0.01, 1.0)\n\njulia> new_input.γ\n0.01\n\n\n\n\n\n","category":"type"},{"location":"manual/#JuSwalbe.inputfile","page":"Manual","title":"JuSwalbe.inputfile","text":"inputfile\n\nAbstract type for all kinds of input files\n\n\n\n\n\n","category":"type"},{"location":"reference/#References-1","page":"References","title":"References","text":"","category":"section"},{"location":"#Home-1","page":"Home","title":"Home","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"JuSwalbe is a package to simulate thin film flows with a lattice Boltzmann approach. The method is introduced in the paper: Lattice Boltzmann method for thin-liquid-film hydrodynamics.","category":"page"},{"location":"devnotes/#Developer-notes-1","page":"Developer notes","title":"Developer notes","text":"","category":"section"}]
}
