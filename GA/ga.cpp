/*
 *	Student Name: Chih-Ang Wang |  AndrewID: chihangw
 *	[To Run]
 *	1. Enter "g++ -o ea gen_ea.cpp" to terminal for compiling.
 *	2. Enter "./ea" to run the program.
 *	3. a file called "report.txt" will be generated with detailed information
 *	   about best, worst and mean fitness in each generation, each run.
 *	4. time taken for each run will be displayed on the terminal screen.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <ctime>
using namespace std;

#define GEN_NUM       100
#define GENE_LEN      25
#define RUN           10
#define SEAD          1234
#define OPT_PHENOTYPE 25
#define OPT_GENOTYPE  0x1ffffff// unsigned integer representing 25 bits, all 1

string chromosomes[GEN_NUM];
string mating_pool[GEN_NUM];
int fitness_arr[GEN_NUM];
ofstream report;

/* get random double number between 0 - 1 */
double get_double_rand(void) {
	return rand() / (double)(RAND_MAX);
}

/* get random integer between bound_l and bound_h */
int get_int_rand(int bound_l, int bound_h) {
	int r_val = bound_l + rand() % bound_h;
	return r_val;
}

/* initialize the GA algorithm's pool */
void init(void) {
	for(int i=0; i<GEN_NUM; i++) {
		chromosomes[i] = to_string(get_int_rand(0, OPT_GENOTYPE));
	}
}

/* sum all the 25 bits in the gene */
int calculate_fitness(int gene) {

	if(gene > OPT_GENOTYPE){
		cout << "GENE out of range !!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	}

	int i = 0;
	int mask = 0x1;
	int fitness = 0;
	while(i < 25) {
		int lsb = (gene >> i) & mask;
		if(lsb) {
			fitness += 1;
		}
		i += 1;
	}
	return fitness;
}

/* roulette wheel algorithm */
void roulette_wheel_select(void) {
	int fitness = 0;
	int fitness_sum = 0;
	double fitness_prop[GEN_NUM];
	double ai_array[GEN_NUM];

	for(int i=0; i<GEN_NUM; i++) {
		fitness = calculate_fitness(atoi(chromosomes[i].c_str()));
		fitness_sum += fitness;
		fitness_arr[i] = fitness;
	}

	for(int i=0; i<GEN_NUM; i++) {
		fitness_prop[i] = fitness_arr[i] / (double)fitness_sum;
	}

	for(int i=0; i<GEN_NUM; i++) {
		double sigma_fitness = 0;
		for(int j=0; j<=i; j++) {
			sigma_fitness += fitness_prop[j];
		}
		ai_array[i] = sigma_fitness;
	}

	int curr_member = 0;

	while(curr_member < GEN_NUM) {

		int j = 0;
		double rand_val = get_double_rand();

		while(ai_array[j] < rand_val) {
			j += 1;
		}

		mating_pool[curr_member] = chromosomes[j];
		curr_member += 1;
	}
}

/* pc: probability of crossover */
void crossover(double pc) {

	for(int i=0; i<GEN_NUM; i+=2) {
		/* crossover occurs */
		if(pc > get_double_rand()) {
			/* integer between 1 to 24 when considering L = 25 */
			int cross_pt = get_int_rand(1, GENE_LEN-1);
			unsigned int parent1 = atoi(mating_pool[i].c_str());
			unsigned int parent2 = atoi(mating_pool[i+1].c_str());

			/* mask = 0000.......01111111111111111111
			 *        ^           \___ cross_pt ___/^
			 *       MSB                           LSB
			 */
			unsigned int mask = ~(~0 >> cross_pt << cross_pt);

			unsigned int tail1 = parent1 & mask;
			unsigned int tail2 = parent2 & mask;

			/* crossover the tail here */
			parent1 = parent1 - tail1 + tail2;
			parent2 = parent2 - tail2 + tail1;

			mating_pool[i]   = to_string(parent1);
			mating_pool[i+1] = to_string(parent2);
		}
		else {

		}
	}
}

/* invert offset-bit in gene g */
unsigned int inverse_bit(unsigned int g, int offset) {
	/* mask points to the bit that needs to be inverted */
	unsigned int mask = 1 << offset;
	return (g ^ mask);
}

/* pm: probability of mutation (bit-flip) */
void mutate(float pm) {
	/* traverse each gene in the mating pool */
	for(int i=0; i<GEN_NUM; i++) {

		unsigned int gene = atoi(mating_pool[i].c_str());

		/* traverse each bit in the genotype */
		for(int j=0; j<GENE_LEN; j++) {

			/* mutation occurs, inverse the current j-bit */
			if(pm > get_double_rand()) {
				gene = inverse_bit(gene, j);
			}
		}

		mating_pool[i] = to_string(gene);
	}
}

/* replace the next generation by those in the mating pool */
void replace(void) {
	for(int i=0; i<GEN_NUM; i++) {
		chromosomes[i] = mating_pool[i];
	}
}

/* search if the optimum phenotype has existed */
bool search_opt(void) {
	for(int i=0; i<GEN_NUM; i++) {
		if(fitness_arr[i] == OPT_PHENOTYPE) {
			return true;
		}
	}
	return false;
}

/* print out min, max, mean fitness of the current generation */
void print_report(int gen) {
	int max=0, min=0, sum=0;
	double mean;

	int curr_val = fitness_arr[0];

	max = min = sum = curr_val;

	/* find the mean, max, min in gene pool */
	for(int i=1; i<GEN_NUM; i++) {

		curr_val = fitness_arr[i];

	  	if(curr_val >= max) {
	  		max = curr_val;
	  	}
	  	if(curr_val <= min) {
	  		min = curr_val;
	  	}
	  	sum += curr_val;
	}
	mean = (double)sum / GEN_NUM;

	/* print the max, min, mean to output.txt */
	report << "Generation[" << gen << "] best: " << max << " | worst: " << min\
	       << " | mean: " << mean << "\n";
}

int main(int argc, char *argv[]) {

	int generation, hit = 0;
	double duration = 0;
	clock_t start;
	report.open("report.txt");
	/* seeding the random number generator */
	srand(SEAD);

	/* 10 run to calculate the standard deviation and the mean
	   of time taken for each round */
	for(int i=0; i<RUN; i++) {

		/* reinitialize generation to first generation */
		init();
		generation = 1;
		hit = 0;

		report << "RUN[" << i << "] Begins Here ----------------------" <<
				"------" << endl;
		/* start counting time */
		start = clock();

		/* start EA algorithm */
		while(generation <= 100) {

			roulette_wheel_select();
			crossover(0.7);
			mutate(1/GENE_LEN);
			replace();

			print_report(generation);

			if(search_opt() == true) {
				cout << "Optimum found at " << generation << "th generation!"
				     << endl;
				hit = 1;
				break;
			}

			generation += 1;
		}
		if(!hit) {
			cout << "EA reaches the end of generation limit." << endl;
		}

		/* stop the timer and calculate the duration */
		duration = (clock() - start) / (double)CLOCKS_PER_SEC;
		cout << "Time taken: " << duration << " seconds.\n" << endl;

		report << "RUN[" << i << "] Ends Here ------------------------" <<
		    "------\n\n" << endl;
	}

	report.close();
	return 0;
}