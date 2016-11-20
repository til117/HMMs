// g++ -std=c++11

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>

using namespace std;

typedef vector< vector<double> > matrix;

struct the_struct
{
	matrix value;
	matrix index;
};

// Function Declaration
string makestr (matrix);
matrix makemat (string);
matrix mult (matrix, matrix);
vector<int> makeseq (string input);
matrix getcol (matrix input, int col);
matrix dotmult (matrix M1, matrix M2);
matrix repmat(double value, int row, int col);
int max_index(matrix input);
double max_value(matrix input);
matrix delta1 (matrix B, matrix P, vector<int> Obs);
the_struct viterbi(matrix A, matrix B, matrix P, vector<int> Obs);
vector<int> viterbi_path(matrix value, matrix index);
the_struct viterbi_step(matrix A, matrix B, matrix prev_delta, int Obs_);

// Main
int main(void)
{
	string input;
	matrix A; matrix B; matrix P; matrix out; vector<int> Obs;

	getline(cin, input); A = makemat(input);
	getline(cin, input); B = makemat(input);
	getline(cin, input); P = makemat(input);
  getline(cin, input); Obs = makeseq(input);

  the_struct output1;
  output1 = viterbi(A, B, P, Obs);
  vector<int> output2 = viterbi_path(output1.value, output1.index);

  for (int i = 1; i < output2.size()+1; i++)
  {
    cout << output2[output2.size()-i] << " ";
  }

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

// Return max index of vector
int max_index(matrix input)
{
	return(distance(input[0].begin(), max_element(input[0].begin(), input[0].end())));
}

// Return max value of vector
double max_value(matrix input)
{
	return *max_element(input[0].begin(), input[0].end());
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

// Creates a matrix with double value
matrix repmat(double value, int row, int col)
{
	matrix output (row, vector<double> (col) );
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			output[i][j] = value;
		}
	}
	return output;
}

// Initialize delta
matrix delta1 (matrix B, matrix P, vector<int> Obs)
{
	matrix output (1, vector<double> (P[0].size()) );
	output = dotmult(P,getcol(B,Obs[0]));
	return output;
}

// Obtain final deltas and indexes for every timestep
the_struct viterbi(matrix A, matrix B, matrix P, vector<int> Obs)
{
	matrix prev_delta = delta1 (B, P, Obs);
	the_struct v_out;
	the_struct output;
	matrix index_size (Obs.size()-1, vector<double> (prev_delta.size()) );
	matrix value_size (Obs.size(), vector<double> (prev_delta.size()) );
	output.value = value_size;
	output.index = index_size;

	output.value[0] = prev_delta[0];

	for (int i = 1; i < Obs.size(); i++)
	{
		v_out = viterbi_step(A, B, prev_delta, Obs[i]);
		prev_delta = v_out.value;
		output.value[i] = prev_delta[0];
		output.index[i-1] = v_out.index[0];
	}

	return output;
}


// Performs a step in viterbi algorithm returning deltas and indexes for a timestep as a struct
the_struct viterbi_step(matrix A, matrix B, matrix prev_delta, int Obs_)
{
	the_struct output;
	matrix value (1, vector<double> (prev_delta[0].size()) );
	matrix index (1, vector<double> (prev_delta[0].size()) );
	matrix tmp (1, vector<double> (prev_delta[0].size()) );
	matrix col_b;
	for (int i = 0; i < prev_delta[0].size(); i++)
	{
		col_b = repmat(B[i][Obs_],prev_delta[0].size(),1);
		tmp = dotmult(prev_delta,getcol(A,i));
		tmp = dotmult(tmp, col_b);
		value[0][i] = max_value(tmp);
		index[0][i] = max_index(tmp);
	}
	output.value = value;
	output.index = index;
	return output;
}

// Returns path
vector<int> viterbi_path(matrix value, matrix index)
{
	vector<int> output;
	int value_size = value[0].size();
	matrix tmp (1, vector<double> (value_size) );

	for (int i = 0; i < value_size; i++)
	{
		tmp[0][i] = value[value_size-1][i];
	}
	int T_end = max_index(tmp);
	output.push_back(T_end);
	int next_state = T_end;

	for(int i = 1; i < index.size() + 1; i++)
	{
		next_state = index[index.size()-i][next_state];
		output.push_back(next_state);
	}
	return output;
}
