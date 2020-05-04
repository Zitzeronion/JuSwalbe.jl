var documenterSearchIndex = {"docs":
[{"location":"manual/#Manual-1","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"manual/#","page":"Manual","title":"Manual","text":"Short and clean instructions how to use the Package. Functions and methods can be found here with usage examples.","category":"page"},{"location":"manual/#","page":"Manual","title":"Manual","text":"Modules = [JuSwalbe]\nOrder = [:function, :type, :struct] ","category":"page"},{"location":"manual/#JuSwalbe.calc_equilibrium_distribution","page":"Manual","title":"JuSwalbe.calc_equilibrium_distribution","text":"calc_equilibrium_distribution(mom:macroquant64_1d, gravity=0)\n\nCalculates the equilibrium distributions based on the macroscopic quantities mom and gravitational acceleration gravity.\n\nThe equilibrium distribtutions are at the heart of the lattice Boltzmann method. As the expansion is made around the equilibrium, lattice Boltzmann always assumes that the flow field is close to equilibrium. Therefore the equilibrium is calculated from macroscopic quantities, i.e. height .height and velocity .velocity. There are plenty of ways to calculate them I mainly the ones derived in Eq.(13) of Study of the 1D lattice Boltzmann shallow water equation and its coupling to build a canal network.\n\nMath\n\nDefining equations in one spatial dimensions are expressed in the following way ``\\begin{align}  f0^{eq} &=  h - \\frac{1}{2v^2}gh^2 - \\frac{1}{v^2}hu^2 \\  f1^{eq} &=  \\frac{1}{4v^2}gh^2 + \\frac{1}{2v}hu + \\frac{1}{2v^2}hu^2 \nf_2^{eq} &=  \\frac{1}{4v^2}gh^2 - \\frac{1}{2v}hu + \\frac{1}{2v^2}hu^2  \\end{align}``\n\n\n\n\n\n","category":"function"},{"location":"manual/#JuSwalbe.readinput-Tuple{Any}","page":"Manual","title":"JuSwalbe.readinput","text":"readinput(file)\n\nReads input parameters from a file.\n\nThe expected amount of parameters can be addressed with inputconstants. For now it expects seven values for different runtime constants.\n\nExample\n\njulia> using JuSwalbe, DelimitedFiles\n\njulia> args = [\"Lattice_points_x\" 10; \"Lattice_points_y\" 5; \"Max_run_time\" 1000; \"Output_dump\" 100; \"gravity\" 0.0; \"surface_tension\" 0.01; \"slippage\" 1.0] # Generate a text file with input\n7×2 Array{Any,2}:\n \"Lattice_points_x\"    10\n \"Lattice_points_y\"     5\n \"Max_run_time\"      1000\n \"Output_dump\"        100\n \"gravity\"              0.0\n \"surface_tension\"      0.01\n \"slippage\"             1.0\n\njulia> writedlm(\"test.txt\", args)\n\njulia> test = readinput(\"test.txt\")\nJuSwalbe.inputconstants(10, 5, 1000, 100, 0.0, 0.01, 1.0)\n\njulia> test.lx\n10\n\njulia> test.γ\n0.01\n\njulia> test.γ + test.δ\n1.01\n\njulia> isa(test.lx + test.gravity, Int32)\nfalse\n\njulia> rm(\"test.txt\")\n\n\n\n\n\n","category":"method"},{"location":"manual/#JuSwalbe.inputconstants","page":"Manual","title":"JuSwalbe.inputconstants","text":"inputconstants = new(lx, ly, maxruntime, dumping, gravity, γ, δ)\n\nStruct containing input parameters.\n\nContains .lx lattice points in x-direction, .ly lattice points in y-direction.  Other fields are .maxruntime for the maximal number of time steps and .dumping to limit the number of output files. On top of these there are physical quantities such as .gravity, .γ and .δ  for the values of gravitational acceleration, fluids surface tension and the slip length. The example relates to an quadratic lattice 20 times 20 lattice units in area.  Run for 100 lattice Boltzmann time steps only printing output every 10 time steps. Having no gravity and a surface tension of 0.01 and a slip length of 1.  \n\nExample\n\njulia> using JuSwalbe\n\njulia> new_input = JuSwalbe.inputconstants(20, 20, 100, 10, 0.0, 0.01, 1.0)\nJuSwalbe.inputconstants(20, 20, 100, 10, 0.0, 0.01, 1.0)\n\njulia> new_input.γ\n0.01\n\n\n\n\n\n","category":"type"},{"location":"manual/#JuSwalbe.inputfile","page":"Manual","title":"JuSwalbe.inputfile","text":"inputfile\n\nAbstract type for all kinds of input files\n\n\n\n\n\n","category":"type"},{"location":"reference/#References-1","page":"References","title":"References","text":"","category":"section"},{"location":"#Home-1","page":"Home","title":"Home","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"JuSwalbe is a package to simulate thin film flows with a lattice Boltzmann approach. The method is introduced in the paper: Lattice Boltzmann method for thin-liquid-film hydrodynamics.","category":"page"},{"location":"devnotes/#Developer-notes-1","page":"Developer notes","title":"Developer notes","text":"","category":"section"}]
}
