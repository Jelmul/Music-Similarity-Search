/*
	Written by Jelle Mulyadi, 2021
*/

#include "ED.h"
#include "../General/Tools.h";
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>

//Algorithms and code based on:
// - Levenshtein, Vladimir I. "Binary codes capable of correcting deletions, insertions, and reversals." Soviet physics doklady. Vol. 10. No. 8. 1966.
// - Wagner, Robert A., and Michael J. Fischer. "The string-to-string correction problem." Journal of the ACM (JACM) 21.1 (1974): 168-173.
// - Hirschberg, Daniel S. "A linear space algorithm for computing maximal common subsequences." Communications of the ACM 18.6 (1975): 341-343.

ED::ED(Database* db, int k)
{
	this->db = db;
	this->k = k;
}

//Based on Levenshtein, Wagner-Fischer and Hirschberg algorithms
int ED::WagnerFischer(std::string query, std::string candidate)
{
	std::vector<std::string> result;
	std::vector<std::vector<int>> matrix;
	matrix.push_back({});
	matrix.push_back({});
	//Initialize matrix
	for (int i = 0; i <= query.size(); i++)
	{
		matrix[0].push_back(i);
		matrix[1].push_back(0);
	}
	//Calculate values
	for (int i = 0; i < candidate.size(); i++)
	{
		for (int j = 0; j <= query.size(); j++)
		{
			if (j == 0)
				matrix[1][j] = matrix[0][j] + 1;
			else
			{
				int val = 0;
				if (candidate[i] != query[j - 1])
					val = 1;
				matrix[1][j] = std::min({
					(matrix[0][j - 1] + val),
					(matrix[0][j] + 1),
					(matrix[1][j - 1] + 1) });
			}
		}
		matrix[0] = matrix[1];
	}
	return matrix[1][query.size()];
}

std::vector<Result> ED::SearchSequence(std::string query)
{
#pragma region Searching
	std::vector<Result> result(db->db.size());
	std::clock_t start = std::clock();
#pragma omp parallel for
	for (int i = 0; i < db->db.size(); i++)
	{
		Entry dbentry = db->db[i];
		std::string entry = dbentry.sequence;
		double score = WagnerFischer(query, entry);
		result[i] = Result(dbentry.index, score, true);
	}
	double duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << ";;"; //2 x empty
	std::cout << duration << ";"; //Searching
#pragma endregion

#pragma region Sorting results
	start = std::clock();
	std::sort(result.rbegin(), result.rend());
	duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";\n"; //Sorting
#pragma endregion
	return result;
}
