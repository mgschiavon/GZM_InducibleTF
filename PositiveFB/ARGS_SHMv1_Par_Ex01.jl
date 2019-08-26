# Kinetic parameters
p = Dict([
    :gY  => 0.01,      # Degradation/dilution rate of the output ([1/min])
    :mY  => 0.38,      # Maximum synthesis rate ([nM/min])
	:aX  => 0.0052,     # Related to basal synthesis of the naked promoter ([0,1])
	:nO  => 1.5,       # Promoter occupancy nonlinearity (Hill coefficient)
	:KO  => 20,        # Related to TF-promoter dissociation rate ([nM])
]);