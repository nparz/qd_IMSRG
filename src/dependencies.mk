bin/adams_ode.o: bin/ME_general.o
bin/commutators.o: bin/ME_general.o
bin/HFQDmod.o: bin/ME_general.o bin/ME_general.o bin/ME_general.o bin/ME_general.o bin/ME_general.o
bin/IMSRG_tools.o: bin/ME_general.o bin/commutators.o
bin/magnus_IMSRG.o: bin/operators.o bin/commutators.o bin/ME_general.o bin/IMSRG_tools.o bin/commutators.o bin/ME_general.o bin/IMSRG_tools.o bin/ME_general.o bin/IMSRG_tools.o bin/ME_general.o bin/IMSRG_tools.o bin/ME_general.o bin/IMSRG_tools.o
bin/ME_general.o:
bin/operators.o: bin/ME_general.o
