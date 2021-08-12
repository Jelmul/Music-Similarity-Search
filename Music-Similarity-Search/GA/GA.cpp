/*
	Written by Jelle Mulyadi, 2021
*/

#include "GA.h"
#include "../General/Scoring.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>

//Algorithms and code based on:
// - Needleman, Saul B., and Christian D. Wunsch. "A general method applicable to the search for similarities in the amino acid sequence of two proteins." Journal of molecular biology 48.3 (1970): 443-453.
// - Gotoh, Osamu. "An improved algorithm for matching biological sequences." Journal of molecular biology 162.3 (1982): 705-708.
// - Hirschberg, Daniel S. "A linear space algorithm for computing maximal common subsequences." Communications of the ACM 18.6 (1975): 341-343.

GA::GA(Database* db, int k)
{
	this->db = db;
	this->k = k;

	//Construct substitution matrix
	for (int i = 0; i < 53; i++)
		for (int j = 0; j < 53; j++)
			substitutionMatrix[i][j] = IntervalScore(db->alphabet[i], db->alphabet[j]);
}

//Based on Needleman-Wunsch, Gotoh and Hirschberg algorithms
int GA::Gotoh(std::string query, std::string candidate)
{
	std::vector<std::vector<int>> Ix;
	std::vector<std::vector<int>> Iy;
	std::vector<std::vector<int>> B;

	//Initialize matrices
	for (int i = 0; i < 2; i++)
	{
		Ix.push_back(std::vector<int>(query.size() + 1, 0));
		Iy.push_back(std::vector<int>(query.size() + 1, 0));
		B.push_back(std::vector<int>(query.size() + 1, 0));
	}

	int d = 2; //Gap opening cost
	int e = 1; //Gap extension cost
	int value = -d;
	for (int i = 1; i <= query.size(); i++)
	{
		B[0][i] = value;
		Ix[0][i] = -INFINITY;
		Iy[0][i] = -INFINITY;
		value -= e;
	}

	//Calculate values
	value = -d;
	for (int i = 0; i < candidate.size(); i++)
	{
		for (int j = 0; j <= query.size(); j++)
		{
			if (j == 0)
			{
				B[1][j] = value;
				Ix[1][j] = -INFINITY;
				Iy[1][j] = -INFINITY;
				value -= e;
			}
			else
			{
				int M = B[0][j - 1] + substitutionMatrix[CharacterIndex(candidate[i])][CharacterIndex(query[j - 1])];
				B[1][j] = std::max({ Ix[0][j], Iy[1][j - 1],  M });
				Ix[1][j] = std::max((M - d), (Ix[0][j] - e));
				Iy[1][j] = std::max((M - d), (Iy[1][j - 1] - e));
			}
		}
		B[0] = B[1];
		Ix[0] = Ix[1];
		Iy[0] = Iy[1];
	}
	return B[1][query.size()];
}

std::vector<Result> GA::SearchSequence(std::string query)
{
#pragma region Searching
	std::vector<Result> result(db->db.size());
	std::clock_t start = std::clock();
#pragma omp parallel for
	for (int i = 0; i < db->db.size(); i++)
	{
		Entry dbentry = db->db[i];
		std::string entry = dbentry.sequence;
		double score = Gotoh(query, entry);
		result[i] = Result(dbentry.index, score, true);
	}
	double duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << ";;"; //2 x empty
	std::cout << duration << ";"; //Searching
#pragma endregion

#pragma region Sorting results
	start = std::clock();
	std::sort(result.begin(), result.end());
	duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";\n"; //Sorting
#pragma endregion
	return result;
}
