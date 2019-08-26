# Full Saturating Model (v01)

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	myODE = @ode_def begin
		dXi = (iM * (1 - (mY * exp(-kB * (((Xa + (aX * Xi))^nO)/(((Xa + (aX * Xi))^nO)+(KO^nO))) / iM)))) - (gY * Xi) + (kM * Xa) - (kP * H * Xi)
		dXa =  - (gY * Xa) - (kM * Xa) + (kP * H * Xi)
		dH  =    (gY * Xa) + (kM * Xa) - (kP * H * Xi)
	end gY kM kP iM mY kB aX nO KO;
	
	# dHT/dt = 0;
	# KX = (gY+kM)/kP
end