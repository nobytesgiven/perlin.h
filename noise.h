#include <math.h>
#include <stdlib.h>

// Based on GDC 2017 talk "Math for Game Programmers: Noise-Based RNG"
// SquirrelNoise5 - Squirrel's Raw Noise utilities (version 5)
// Code available under CC-BY-3.0 US
double SquirrelNoise5( int positionX, unsigned int seed )
{
	unsigned int SQ5_BIT_NOISE1 = 0xd2a80a3f;	// 11010010101010000000101000111111
	unsigned int SQ5_BIT_NOISE2 = 0xa884f197;	// 10101000100001001111000110010111
	unsigned int SQ5_BIT_NOISE3 = 0x6C736F4B; // 01101100011100110110111101001011
	unsigned int SQ5_BIT_NOISE4 = 0xB79F3ABB;	// 10110111100111110011101010111011
	unsigned int SQ5_BIT_NOISE5 = 0x1b56c4f5;	// 00011011010101101100010011110101

	unsigned int mangledBits = (unsigned int) positionX;
	mangledBits *= SQ5_BIT_NOISE1;
	mangledBits += seed;
	mangledBits ^= (mangledBits >> 9);
	mangledBits += SQ5_BIT_NOISE2;
	mangledBits ^= (mangledBits >> 11);
	mangledBits *= SQ5_BIT_NOISE3;
	mangledBits ^= (mangledBits >> 13);
	mangledBits += SQ5_BIT_NOISE4;
	mangledBits ^= (mangledBits >> 15);
	mangledBits *= SQ5_BIT_NOISE5;
	mangledBits ^= (mangledBits >> 17);
	return mangledBits/2147483648.0f-1.0f; // 2^32
}

double interpolate(double x) {
	return x*x*x*(x*(6*x-15)+10);
}


void pseudo_gradient(double** grad, int *coords, unsigned int n) { // random n dimensional gradient vector
	const int primes[] = {
		498493943,
		398491637,
		198491317,
		89851219,
		29849987,
		6542989,
		544139,
		15971,
		9679,
		823
	};
	
	double norm = 0.0f;
	int uniq = coords[0];
	
	for (int i=1; i<(int)n; i++) {
		uniq += primes[i%3]*coords[i];
	}
	
	double pseudo_noise_unused;
	double pseudo_noise;
	for (int i=0; i<(int)n; i++) {
		if ((i&1) == 0) {
			double U1 = 0.5f*(1.0f+SquirrelNoise5(coords[i], uniq+2*n*i)); // Techincally uniform (0,1) random values
			double U2 = 0.5f*(1.0f+SquirrelNoise5(coords[i], uniq+2*n*i+1));
			
			// Box-Muller transform to generate normally distributed random numbers Xi
			// The resulting vectors (X0, X1,...,Xn) are uniformly distributed on the n-sphere surface (After normalization)
			pseudo_noise = sqrt(-2*log(U1))*cos(2*PI*U2);			
			pseudo_noise_unused = sqrt(-2*log(U2))*cos(2*PI*U1);
		} else {
			pseudo_noise = pseudo_noise_unused;
		}

		(*grad)[i] = pseudo_noise;
		norm += pseudo_noise*pseudo_noise;
	}
	norm = sqrt(norm);

	for (int i=0; i<(int)n; i++) (*grad)[i] /= norm; // Normalize vector
	
	return;
}

double noise(const double* coords, const unsigned int n) { // n = dimension
	double ret = 0.0f;
	double* distance_vec = malloc(n*sizeof(double));
	double* grad = malloc(n*sizeof(double)); // Gradient vector
	int* hyperpoint = malloc(n*sizeof(double)); // Hypercube point

	for (int i=0; i<pow(2,n); i++) {
		double pos = 1.0f;
		double dot_product = 0.0f;
		
		/*
				  N-DIMENSIONAL INTERPOLATION
							n
						   2
		f(X0,X1,...,Xn) =  Σ  An*Q0(X0)*Q1(X1)*...*Qn(Xn)
						  j=0

		An is the gradient*distance vector dot product
		of the hypercube point n
			 __
			 |1-p(x) , i=1  Where i is i-th bit of the
		Qi = |              hypercube point n (starting
			 |p(x)   , i=0  from the right)
			 +-

		p(x) is the interpolation function. In our case
		it is:
							5    4    3
				   p(x) = 6x -15x +10x
				   
		*/
		
		for (int k=0; k<(int)n; k++) {
			int isset = (i>>(n-k-1))&1; // Is k-th bit set? (starting from the left)
			
			if (isset) {
				hyperpoint[k] = (int)floor(coords[k]) + 1;
				pos *= interpolate(coords[k]-(int)coords[k]);
			} else {
				hyperpoint[k] = (int)floor(coords[k]);
				pos *= 1-interpolate(coords[k]-(int)coords[k]);
			}
			distance_vec[k] = coords[k] - hyperpoint[k];
		}

		pseudo_gradient(&grad, hyperpoint, n);
		for (int j=0; j<(int)n; j++) {
			dot_product += grad[j]*distance_vec[j];
		}

		ret += dot_product*pos;
	}

	free(distance_vec);
	free(grad);
	free(hyperpoint);
	
	// https://digitalfreepen.com/2017/06/20/range-perlin-noise.html
	// Maximum possible noise value is sqrt(N/4), where N=dimension
	return ret*sqrt((double)4/n); // Map to (-1,1) range. 
}

