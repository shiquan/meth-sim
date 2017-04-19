all:
	gcc -g -lm -lhts -lz -o methsim number.c sequence.c stats.c methsim.c
