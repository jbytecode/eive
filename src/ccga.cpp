#include <Rcpp.h>
using namespace Rcpp;

//' @name cga_generate_chromosome
//' @title Generate Chromosome
//' @description Generate a binary vector using a probability vector
//' 	This function is not directly called by user. CGAs (Compact genetic algorithms)
//' 	sample chromosomes using this probability vector. A probability vector
//' 	contains[P1, P2, ..., PN] and the function generates and returns a chromosome[B1, B2, ..., BN].
//' 	The probability of BK having the value of 1 is PK. So, it has more chance to have
//' 	[1, 1, 1, 0, 0] rather than [0, 0, 0, 1, 1] when the probability vector is
//'  	[0.9, 0.9, 0.9, 0.1, 0.1]. 
//' @param proc_vec Vector of probabilities
//' @param vect Vector of bits.
//' @return Mutates the vect. Returns null.
//' @export
// [[Rcpp::export]]
void cga_generate_chromosome(NumericVector prob_vec, NumericVector vect)
{
	int len_prob_vec = prob_vec.length();
	NumericVector unifs = runif(len_prob_vec);
	int i;
	for (i = 0; i < len_prob_vec; i++)
	{
		if (unifs[i] < prob_vec[i])
		{
			vect[i] = 1;
		}else{
			vect[i] = 0;
		}
	}
}

//' @name cga
//' @title Compact Genetic Algorithm
//' @description Performs a Compact Genetic Algorithm (CGA) search
//' 	for a given chromosome size, population size (mutation rate), 
//'     and an objective function. 
//' @param chsize Number of bits.
//' @param popsize Size of population. The value is used for mutating 
//'		the probability vector by 1/popsize. 
//' @param evalFunc Objective function.
//' @return Binary vector of size chsize.
//' @export
// [[Rcpp::export]]
NumericVector cga(int chsize, int popsize, Function evalFunc)
{
	NumericVector prob_vec = rep(0.5, chsize);
	NumericVector chromosome1(chsize);
	NumericVector chromosome2(chsize);
	NumericVector winner, loser;
	NumericVector cost1, cost2;
	int i, t;
	double mutation = 1.0 / (double)popsize;
	while (1)
	{
		cga_generate_chromosome(prob_vec, chromosome1);
		cga_generate_chromosome(prob_vec, chromosome2);
		cost1 = evalFunc(chromosome1);
		cost2 = evalFunc(chromosome2);
		winner = chromosome1;
		loser = chromosome2;
		if (cost2[0] < cost1[0])
		{
			winner = chromosome2;
			loser = chromosome1;
		}

		for (i = 0; i < chsize; i++)
		{
			if (winner[i] != loser[i])
			{
				if (winner[i] == 1)
				{
					prob_vec[i] = prob_vec[i] + mutation;
				}
				else
				{
					prob_vec[i] = prob_vec[i] - mutation;
				}
			}
		}

		t = 0;
		for (i = 0; i < chsize; i++)
		{
			if (prob_vec[i] <= 0.001 || prob_vec[i] >= 0.999)
			{
				t++;
			}
		}
		if (t == chsize)
		{
			break;
		}
	}
	return (prob_vec);
}
