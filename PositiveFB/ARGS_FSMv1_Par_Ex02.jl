# Kinetic parameters
p = Dict([
    :gY  => 0.01,      # Degradation/dilution rate of the output ([1/min])
    :kM  => 23.49,     # X:H unbinding rate
    :kP  => 1,     	   # X:H binding rate
    :iM  => 0.4,       # Maximum synthesis rate given the gene & translocation rate
	:mY  => 0.99,      # Related to basal synthesis of the naked promoter ([0,1])
	:kB  => 0.58,      # Related to efficiency rate of the transcription factor
	:aX  => 0.039,     # Related to basal activity of free (non-active) TF ([0,1])
	:nO  => 1.5,       # Promoter occupancy nonlinearity (Hill coefficient)
	:KO  => 4.8,       # Related to TF-promoter dissociation rate ([nM])
]);