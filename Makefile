SRC= dgp.cpp  main.cpp  vec.cpp
INC= dgp.h  fobj_sph.h  option_handler.h  sph.h  vec.h

# release
sph.bin: $(SRC) $(INC)
	g++ -O3 -std=c++17 -Wall -Werror $(SRC) -I./ -o sph.bin -fopenmp -lgsl -lm -lgslcblas
	# pgc++ -Werror -std=c++11 -fast -ta=multicore -Minfo=all $(SRC) -I./ -o sph.bin -lgsl -lm -lgslcblas

# debug 
sph.dbg: $(SRC) $(INC)
	g++ -O0 -std=c++17 -Wall -Werror -g  $(SRC) -I./ -o sph.dbg -fopenmp -lgsl -lm -lgslcblas
	# pgc++ -Werror -std=c++11 -fast -ta=multicore -Minfo=all $(SRC) -I./ -o sph.bin -lgsl -lm -lgslcblas

run: sph.bin
	./sph.bin -dgp DATA_EPSD_01_DMAX_60/1bdk.nmr -ftol 1E-2
	
valgrind: sph.bin
	valgrind -v --tool=memcheck --leak-check=full --show-leak-kinds=all ./sph.bin -dgp DATA_EPSD_01_DMAX_60/1bdk.nmr
	
clean:
	rm -rf sph.bin
	rm -rf sph.dbg
