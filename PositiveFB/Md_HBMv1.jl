# Hill + Basal Model (v01)

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	myODE = @ode_def begin
		dXi = (mY * (aX + ((1 - aX) * (((H * Xi)^nO)/(((H * Xi)^nO) + (KO^nO)))))) + (bX * Xi) - (gY * Xi)
		dXa =  0
		dH  =  0
	end gY mY aX nO KO bX;
end