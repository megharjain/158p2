all: mmm_omp heat_distribution_serial heat_distribution_omp

mmm_omp: mmm_omp.c matrix_checksum.c
	gcc -o mmm_omp mmm_omp.c matrix_checksum.c -Wall -Werror -fopenmp -O2

heat_distribution_serial: heat_distribution_serial.c matrix_checksum.c
	gcc heat_distribution_serial.c matrix_checksum.c -o heat_distribution_serial -Wall -Werror -O2
	
heat_distribution_omp: heat_distribution_omp.c matrix_checksum.c
	gcc heat_distribution_omp.c matrix_checksum.c -o heat_distribution_omp -Wall -Werror -fopenmp -O2
	
	
clean:
	rm -f mmm_omp heat_distribution_serial heat_distribution_omp
