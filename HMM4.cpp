// g++ -std=c++11

#include <iostream>
#include <float.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

using namespace std;

// Create a matrix which is a vector of vectors
typedef vector< vector<double> > matrix;

// Class definition
class Lamda {
	public:
		matrix A, B, P;
		matrix alpha, beta, gamma;
		vector<int> Obs;
		vector<double> c;
		vector<matrix> digamma;
		double logprob, oldlogprob;
		int T, N, M, iters, maxiters;

		Lamda(matrix AA, matrix BB, matrix PP, vector<int> ObsObs);
		void alpha_pass();
		void beta_pass();
		void gamma_pass();
		void re_estimate();
		void logprob_func();
		void iterate();
};


// Function Declaration
string makestr (matrix);
matrix makemat (string);
matrix mult (matrix, matrix);
vector<int> makeseq (string input);

// Main
int main(void)
{
	string input;
	matrix A_; matrix B_; matrix P_; vector<int> Obs_;

	// Get input
	getline(cin, input); A_ = makemat(input);
	getline(cin, input); B_ = makemat(input);
	getline(cin, input); P_ = makemat(input);
  getline(cin, input); Obs_ = makeseq(input);

	// Pass input to class and iterate
	Lamda system(A_, B_, P_, Obs_);
	system.iterate();

	// Output
	cout << makestr(system.A) << "\n" << makestr(system.B);


}

// Functions

// Convert matrix to string
string makestr (matrix input)
{
	int row = input.size();
	int col = input[0].size();
	string output;
	string space = " ";
	output += to_string(row) + space + to_string(col);

	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			stringstream ss;
			ss << input[i][j];
			string tmp = ss.str();
			if (tmp == "0")
			{
				tmp = "0.0";
			}
			output += space + tmp;
		}
	}
	return output;
}

// Convert string to matrix
matrix makemat (string input)
{
	istringstream ss(input);
	vector <string> record;

	while(ss)
	{
		string s;
		{
			if (!getline( ss, s, ' ' )) break;
  			record.push_back( s );
		}
	}
	int row = stoi(record[0]);
	int col = stoi(record[1]);
	matrix output (row, vector<double> (col));

	for(int i = 0; i < row; i++)
	{
		for(int j = 0; j < col; j++)
		{
			output[i][j] = stod(record[i*col + j + 2]);
		}
	}
	return output;
}

// Multiply matrix
matrix mult (matrix M1, matrix M2)
{
	matrix Mout (M1.size(),vector<double>(M2[0].size()));

	for(int i = 0; i < M1.size(); i++)
	{
		for(int j = 0; j < M2[0].size(); j++)
		{
			for(int k = 0; k < M1[0].size(); k++)
			{
				Mout[i][j] += M1[i][k]*M2[k][j];
			}
		}

	}
	return Mout;
}

// Make sequence of string
vector<int> makeseq (string input)
{
	istringstream ss(input);
	vector <string> record;

	while(ss)
	{
		string s;
		{
			if (!getline( ss, s, ' ' )) break;
  			record.push_back( s );
		}
	}
	int length = stoi(record[0]);
	vector<int> output(length);

	for (int i = 0; i < length; i++)
	{
		output[i] = stoi(record[i+1]);
	}
	return output;
}

// Get a specified column from matrix
matrix getcol (matrix input, int col)
{
	matrix output (input.size(), vector<double> (1));

	for (int i = 0; i < input.size(); i++)
	{
		output[i][0] = input[i][col];
	}
	return output;
}

// Perform dot multiplication
matrix dotmult (matrix M1, matrix M2)
{
	matrix output (1, vector<double> (M2.size()) );
	for(int i = 0; i < M2.size(); i++)
	{
		output[0][i] = M1[0][i]*M2[i][0];
	}
	return output;
}

// Define Class Functions

// Initialize parameters
Lamda::Lamda(matrix AA, matrix BB, matrix PP, vector<int> ObsObs)
{
	maxiters = 100;
	iters = 0;
	oldlogprob = -DBL_MAX;
	logprob = -DBL_MAX;
	T = ObsObs.size();
	N = AA.size();
	M = BB[0].size();
	A = AA;
	B = BB;
	P = PP;
	Obs = ObsObs;
}

void Lamda::alpha_pass()
{
	// Define alpha and c
	matrix alpha_tot (T, vector<double>(N));
	vector<double> c_tot(T);
	alpha = alpha_tot;
	c = c_tot;

	// Compute alpha_0
	for (int i = 0; i < N; i++)
	{
		alpha[0][i] = P[0][i]*B[i][Obs[0]];
		c[0] += alpha[0][i];
	}

	// Scale alpha_0
	c[0] = 1.0/c[0];
	for (int i = 0; i < N; i ++)
	{
		alpha[0][i] = c[0]*alpha[0][i];
	}

	// Compute alpha_t
	for (int t = 1; t < T; t++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				alpha[t][i] += alpha[t-1][j]*A[j][i];
			}
			alpha[t][i] *= B[i][Obs[t]];
			c[t] += alpha[t][i];
		}

		// Scale alpha_t
		c[t] = 1.0/c[t];
		for (int i = 0; i < N; i++)
		{
			alpha[t][i] *= c[t];
		}
	}
}



// Beta pass
void Lamda::beta_pass()
{
	// Define beta
	matrix beta_tot (T, vector<double>(N));
	beta = beta_tot;

	// Compute beta_T-1
	for (int i = 0; i < N; i ++) beta[T-1][i] = c[T-1];

	// Compute and scale beta_t
	for (int t = T-2; t >= 0; t--)
	{
		for (int i = 0; i < N; i ++)
		{
			for (int j = 0; j < N; j++)
			{
				beta[t][i] += A[i][j]*B[j][Obs[t+1]]*beta[t+1][j];
			}
			beta[t][i] *= c[t]; // Scale
		}
	}
}

// Gamma
void Lamda::gamma_pass()
{
	// Define gamma and digamme
	matrix gamma_tot (T, vector<double>(N));
	vector<matrix> digamma_tot (T-1, matrix(N, vector<double>(N)));
	gamma = gamma_tot;
	digamma = digamma_tot;

	for (int t = 0; t < T-1; t++)
	{
		// Compute denominator
		double denom = 0;
		for (int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				denom += alpha[t][i]*A[i][j]*B[j][Obs[t+1]]*beta[t+1][j];
			}
		}
		// Compute digamma and gamma
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				digamma[t][i][j] = (alpha[t][i]*A[i][j]*B[j][Obs[t+1]]*beta[t+1][j])/denom;
				gamma[t][i] += digamma[t][i][j];
			}
		}
	}
	// Special case for gamma_T-1
	double denom = 0;
	for (int i = 0; i < N; i++)
	{
		denom += alpha[T-1][i];
	}
	for (int i = 0; i < N; i++)
	{
		gamma[T-1][i] = alpha[T-1][i]/denom;
	}
}

// Recalculate
void Lamda::re_estimate()
{
	// Define
	matrix A_tot (N, vector<double>(N));
	matrix B_tot (N, vector<double>(M));
	matrix P_tot (1, vector<double>(N));
	A = A_tot;
	B = B_tot;
	P = P_tot;

	// Reestimate Pi
	for(int i = 0; i < N; i++) P[0][i] = gamma[0][i];

	// Reestimate A
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			double numer = 0, denom = 0;
			for(int t = 0; t < T-1; t++)
			{
				numer+= digamma[t][i][j];
				denom+= gamma[t][i];
			}
			A[i][j] = numer/denom;
		}
	}

	// Reestimate B
	for (int i = 0; i < N; i++)
	{
		for(int j = 0; j < M; j++)
		{
			double numer = 0, denom = 0;
			for(int t = 0; t < T; t++)
			{
				if(Obs[t] == j) {numer+= gamma[t][i];}
				denom += gamma[t][i];
			}
			B[i][j] = numer/denom;
		}
	}
}

// Compute logprob
void Lamda::logprob_func()
{
	logprob = 0;
	for(int i = 0; i < T; i++) logprob += log(c[i]);
	logprob = -logprob;
}

// To iterate or not to iterate
void Lamda::iterate()
{
	while((iters < maxiters) && (logprob >= oldlogprob))
	{
		oldlogprob = logprob;
		alpha_pass();
		beta_pass();
		gamma_pass();
		re_estimate();
		logprob_func();
		iters ++;
	}
}
