/*
	Written by Jelle Mulyadi, 2021
*/

#include "Tools.h"
#include "Entry.h"
#include "Scoring.h"
#include<vector>
#include<algorithm>
#include<bitset>

//Based on: Li, Guoliang, et al. "Pass-join: A partition-based method for similarity joins." arXiv preprint arXiv:1111.7171 (2011).
int LengthAwareED(std::string query, std::string candidate, int startQ, int lengthQ, int startC, int lengthC, int threshold)
{
	//Base cases
	if (lengthQ == 0 && lengthC == 0)
		return 0;
	else if (lengthQ == 0)
		return lengthC;
	else if (lengthC == 0)
		return lengthQ;
	//String A has to be longer than string B, so if A < B swap
	if (lengthQ < lengthC)
	{
		std::string temp = candidate;
		candidate = query;
		query = temp;
		int tempI = startC;
		startC = startQ;
		startQ = tempI;
		tempI = lengthC;
		lengthC = lengthQ;
		lengthQ = tempI;
	}
	std::vector<std::string> result;
	std::vector<std::vector<int>> matrix;
	matrix.push_back({});
	matrix.push_back({});
	//Initialize matrix
	for (int i = 0; i <= lengthQ; i++)
	{
		matrix[0].push_back(i);
		matrix[1].push_back(INT_MAX - 1);
	}

	//Calculate values
	int delta = abs(lengthQ - lengthC);
	for (int i = 1; i < lengthC + 1; i++)
	{
		//Length-aware verification
		int lowerJ = std::max(0, (int)(i - floor((threshold - delta) / 2)));
		int upperJ = std::min(lengthQ, (int)(i + floor((threshold + delta) / 2)));
		int ETcount = 0;
		int ETmax = upperJ - lowerJ + 1;
		//Set first j value
		if (lowerJ == 0)
		{
			matrix[1][lowerJ] = matrix[0][lowerJ] + 1;
		}
		else
		{
			int val = (candidate[startC + i - 1] != query[startQ + lowerJ - 1]);
			matrix[1][lowerJ] = std::min({
				(matrix[0][lowerJ - 1] + val),
				(matrix[0][lowerJ] + 1) });
		}
		if (matrix[1][lowerJ] > threshold)
			ETcount++;
		lowerJ++;
		//Set remaining j values
		for (int j = lowerJ; j <= upperJ; j++)
		{
			int val = (candidate[startC + i - 1] != query[startQ + j - 1]);
			matrix[1][j] = std::min({
				(matrix[0][j - 1] + val),
				(matrix[0][j] + 1),
				(matrix[1][j - 1] + 1) });
			if (matrix[1][j] > threshold)
				ETcount++;
		}
		if (ETcount == ETmax) //Early termination if all cells have edit distance above threshold
			return threshold + 1;
		matrix[0] = matrix[1];
	}
	return matrix[1][lengthQ];
}

int SubstringEditDistance(std::string query, std::string candidate, int startQ, int lengthQ, int startC, int lengthC)
{
	//Initialize matrix
	std::vector<std::vector<int>> matrix;
	matrix.push_back({});
	matrix.push_back({});
	//Initialize matrix
	for (int i = 0; i <= lengthC; i++)
	{
		matrix[0].push_back(i);
		matrix[1].push_back(0);
	}

	//Calculate values
	int min = INT_MAX;
	for (int i = 1; i < lengthQ + 1; i++)
	{
		for (int j = 1; j < lengthC + 1; j++)
		{
			int val = (candidate[startC + j - 1] != query[startQ + i - 1]);
			matrix[1][j] = std::min({
				(matrix[0][j] + 1),
				(matrix[0][j - 1] + val),
				(matrix[1][j - 1] + 1) });
		}
		if (matrix[1][lengthC] < min)
			min = matrix[1][lengthC];
		matrix[0] = matrix[1];
	}
	return min;
}

double SubstringHammingDistance(std::string query, std::string candidate, int startQ, int lengthQ, int startC, int lengthC)
{
	//Generate bitvector substring
	std::bitset<53> bsSub1; //TODO: allow for variable alphabet size
	for (int i = 0; i < lengthC; i++)
	{
		int pos = CharacterIndex(candidate[startC + i]);
		bsSub1[pos] = true;
	}
	//Generate bitvectors substrings in range
	int min = INT_MAX;
	for (int i = 1; i < lengthQ - lengthC + 1; i++)
	{
		std::bitset<53> bsSub2; //TODO: allow for variable alphabet size
		for (int j = 0; j < lengthC; j++)
		{
			int pos = CharacterIndex(query[startQ + i + j]);
			bsSub2[pos] = true;
		}
		std::bitset<53> difference = bsSub1 ^ bsSub2; //TODO: allow for variable alphabet size
		//Popcount
		int t = difference.count();
		if (t < min)
			min = t;
	}
	//Return lowerbound on substring edit distance
	return (double)min / (double)2;
}

bool CompareLength(Entry i, Entry j)
{
	return (i.sequence.length() < j.sequence.length());
}