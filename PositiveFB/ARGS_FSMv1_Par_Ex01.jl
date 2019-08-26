# Kinetic parameters
p = Dict([
    :gY  => 0.01,      # Degradation/dilution rate of the output ([1/min])
    :kM  => 36.1395,   # X:H unbinding rate
    :kP  => 1,     	   # X:H binding rate
    :iM  => 0.0383,    # Maximum synthesis rate given the gene & translocation rate
	:mY  => 0.9946,    # Related to basal synthesis of the naked promoter ([0,1])
	:kB  => 0.6,       # Related to efficiency rate of the transcription factor
	:aX  => 3.5563e-4, # Related to basal activity of free (non-active) TF ([0,1])
	:nO  => 1.5757,    # Promoter occupancy nonlinearity (Hill coefficient)
	:KO  => 2.6626,    # Related to TF-promoter dissociation rate ([nM])
]);