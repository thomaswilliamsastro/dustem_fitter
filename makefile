all:
	f2py -L/usr/lib -llapack -c fortran_funcs.f95 -m fortran_funcs
