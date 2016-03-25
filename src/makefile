FC = gfortran $(FFLAGS)  #-O3 -static #-lacml -llapack

p1 = run_magnus
p2 = calc_TBME
p3 = run_CI
p4 = gs_decouple
p5 = ex_decouple

FFLAGS =  -O3 -fbounds-check -fopenmp  #-CB #-check all
TFLAGS =  -g -fbounds-check -fopenmp

LIBS =  -L/user/local/lib/ -llapack -lblas -larpack

gs: ME_general.o commutators.o adams_ode.o IMSRG_tools.o HFQDmod.o operators.o EOM_IMSRG.o gs_IMSRG.o 
	${FC} ME_general.o commutators.o IMSRG_tools.o HFQDmod.o adams_ode.o operators.o EOM_IMSRG.o gs_IMSRG.o -o ${p4} ${LIBS}

magnus: ME_general.o commutators.o IMSRG_tools.o HFQDmod.o operators.o magnus_IMSRG.o 
	${FC} ME_general.o commutators.o IMSRG_tools.o HFQDmod.o operators.o magnus_IMSRG.o -o ${p1} ${LIBS}

tda: ME_general.o commutators.o adams_ode.o IMSRG_tools.o HFQDmod.o operators.o TDA_IMSRG.o 
	${FC} ME_general.o commutators.o IMSRG_tools.o HFQDmod.o adams_ode.o operators.o TDA_IMSRG.o -o ${p5} ${LIBS}

TBME: ME_general.o get_TBME.o
	${FC} ME_general.o get_TBME.o -o ${p2} ${LIBS}

ME_general.o: ME_general.f90
	${FC} -c ME_general.f90 ${LIBS}

commutators.o:  commutators.f90
	${FC} -c  commutators.f90 ${LIBS}

IMSRG_tools.o: IMSRG_tools.f90
	${FC} -c IMSRG_tools.f90 ${LIBS}

HFQDmod.o: HFQDmod.f90
	${FC} -c HFQDmod.f90 ${LIBS}

operators.o: operators.f90
	${FC} -c operators.f90 ${LIBS}

magnus_IMSRG.o: magnus_IMSRG.f90
	${FC} -c magnus_IMSRG.f90 ${LIBS}

flow_IMSRG.o: flow_IMSRG.f90
	${FC} -c flow_IMSRG.f90 ${LIBS}

gs_IMSRG.o: gs_IMSRG.f90
	${FC} -c gs_IMSRG.f90 ${LIBS}

TDA_IMSRG.o: TDA_IMSRG.f90
	${FC} -c TDA_IMSRG.f90 ${LIBS}

flow_IMSRG3.o: flow_IMSRG3.f90
	${FC} -c flow_IMSRG3.f90 ${LIBS}

adams_ode.o: adams_ode.f90
	${FC} -c adams_ode.f90  ${LIBS}

EOM_IMSRG.o: EOM_IMSRG.f90
	${FC} -c EOM_IMSRG.f90  ${LIBS}

get_TBME.o: get_TBME.f90
	${FC} -c get_TBME.f90  ${LIBS}


CI: simple_shell_model.f90 
	${FC} simple_shell_model.f90 -o ${p3} ${LIBS}

clean:
	rm -f ${p1} $
	rm -f ${p2} $
	rm -f ${p3} $
	rm -f ${p4} $
	rm -f ${p5} $
	rm -f *.o
	rm -f *.mod
	rm -f *~
	rm -f *.su

