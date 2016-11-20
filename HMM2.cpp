// g++ -std=c++11

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

typedef vector< vector<double> > matrix;

// Fnction Declaration
string makestr (matrix);
matrix makemat (string);
matrix mult (matrix, matrix);
vector<int> makeseq (string input);
matrix alpha1 (matrix B, matrix P, vector<int> Obs);
matrix alpha_pass(matrix A, matrix B, matrix P, vector<int> Obs);
matrix alpha_step (matrix A, matrix B, matrix prev_alpha, int Obs_);
double sum_alpha (matrix M);
matrix getcol (matrix input, int col);
matrix dotmult (matrix M1, matrix M2);

// Main
int main(void)
{
	string input;
	matrix A; matrix B; matrix P; matrix out; vector<int> Obs;

	getline(cin, input); A = makemat(input);
	getline(cin, input); B = makemat(input);
	getline(cin, input); P = makemat(input);
  getline(cin, input); Obs = makeseq(input);

  matrix output = alpha_pass(A, B, P, Obs);
  cout << output[0][Obs.size()-1];

  return 0;
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

// Initialize alpha
matrix alpha1 (matrix B, matrix P, vector<int> Obs)
{
	matrix output (1, vector<double> (P[0].size()) );
	output = dotmult(P,getcol(B,Obs[0]));
	return output;
}

// Alpha pass
matrix alpha_pass(matrix A, matrix B, matrix P, vector<int> Obs)
{
	matrix output (1, vector<double> (Obs.size()) );
	matrix prev_alpha = alpha1(B, P, Obs);
	output[0][0] = sum_alpha(prev_alpha);

	for (int i = 1; i < Obs.size(); i++)
	{
		prev_alpha = alpha_step(A, B, prev_alpha, Obs[i]);
		output[0][i] = sum_alpha(prev_alpha);
	}
	return output;
}

// Perform a step in alpha pass
matrix alpha_step (matrix A, matrix B, matrix prev_alpha, int Obs_)
{
	matrix output (1, vector<double> (prev_alpha[0].size()) );
	matrix alpha_ = mult(prev_alpha, A);
	matrix Obs_b = getcol(B, Obs_);
	output = dotmult(alpha_, Obs_b);
	return output;
}

//  Sum alpha matrix
double sum_alpha (matrix M)
{
	double sum = 0;
	for (int i = 0; i < M[0].size(); i++)
	{
		sum += M[0][i];
	}
	return sum;
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
