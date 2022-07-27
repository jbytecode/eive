\name{cga_generate_chromosome}
\alias{cga_generate_chromosome}
\title{
cga_generate_chromosome
}
\description{
Generates vector of zeros and ones for a given probability vector. 
}
\usage{
cga_generate_chromosome(prob_vec, vect)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{prob_vec}{
Vector of probabilities. 
}
\item{vect}{
Integer vector with the same number elements of probabilities. The returned 0-1 vector is held by this variable.
}
}
\details{
This function is not directly called by user. CGA (Compact genetic algorithms) sample chromosomes using this probability vector.
A probability vector contains [P1,P2,...,PN] and the function generates and returns a chromosome [B1,B2,...,BN]. The probability
of BK having the value of 1 is PK. So, it has more chance to have [1,1,1,0,0] than [0,0,0,1,1] when the probability vector is
[0.9,0.9,0.9,0.1,0.1]. 
}
\value{
Returns NULL. The wrapped C++ function returns void. The result is held by the `vect` argument.
}
\author{
Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
Erkin Diyarbakirlioglu <ediyarbakirlioglu@gmail.com>
}


