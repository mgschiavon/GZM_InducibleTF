# Steady state calculation functions

# Julia v.1.1.1

module fn
	# Required libraries
	using DifferentialEquations
	using Statistics

	# Steady state function for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	# OUPUT: ss   - Vector of steady state of the ODE system
	function SS(syst, p, x0, rtol)
		pV = [p[i] for i in syst.params];
		#ss = solve(SteadyStateProblem(syst, x0, pV), SSRootfind());
		#return ss.u
		ss = solve(SteadyStateProblem(syst, x0, pV), DynamicSS(Rodas5(); reltol=rtol));
		return last(ss.u)
	end;
end