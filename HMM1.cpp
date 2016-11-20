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


// Main
int main(void)
{
	string line;
	matrix A; matrix B; matrix P; matrix out;

	getline(cin,line); A = makemat(line);
	getline(cin, line); B = makemat(line);
	getline(cin, line); P = makemat(line);

	out = mult(mult(P,A),B);

	cout << makestr(out);

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
