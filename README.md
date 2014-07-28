irtools
=======

Some tools for IR and 2DIR spectroscopy modelling.

ir2dspectral2.c is highly recommended as it contains all the features to obtain a realistic 2D IR plot,
including the lifetime broadening for the 0->1 and 1->2 transitions.

The lifetime broadening is not provided by ir2dspectra.c. ir2dspectral.c only provides it for the 0->1 transition.

ir2dspectral_many.c works like ir2dspectral2.c but averages the response functions over many molecules.
irspec.cu is the single precision GPU version of ir2dspectral_many.c

Please not that all the previous code requires the frequency fluctuations correlation function (FFCF) to be modelled
by a tri-exponential:
g(t)=c_1*exp(-t/tau_1)+c_2*exp(-t/tau_1t2)+c_3*exp(-t/tau_3)

ir2dspectra_num.c is supposed to avoid any assumption in the functional form of the FFCF by using numerical integration.
So far, the integration scheme used is too trivial to provide any acceptable result. It will be updated soon with a more
accurate integration scheme. Nevertheless, this approach requires a highly converged FFCF.
