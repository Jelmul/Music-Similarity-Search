/*
	Written by Jelle Mulyadi, 2021
*/

#include "BLAST.h"
#include "../General/Scoring.h"
#include <algorithm>
#include <ctime>
#include <iostream>

//Algorithms and code based on:
// - Altschul, Stephen F., et al. "Basic local alignment search tool." Journal of molecular biology 215.3 (1990): 403-410.
// - Altschul, Stephen F., et al. "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic acids research 25.17 (1997): 3389-3402.

BLAST::BLAST(Database* db, int T, int A, int X, int q)
{
	this->db = db;
	this->T = T; //T high scoring word threshold
	this->A = A; //Window size
	this->X = X; //Drop off score
	this->q = q; //Word size

	//Construct substitution matrix
	for (int i = 0; i < 53; i++)
		for (int j = 0; j < 53; j++)
			substitutionMatrix[i][j] = IntervalScore(db->alphabet[i], db->alphabet[j]);
}

template <typename Mat>
void FindHSWs(std::string currentString, int score, int T, int pos, std::unordered_map<std::string, std::vector<HSW>>& index, std::string alphabet, int HSWpos, Mat& substitutionMatrix)
{
	for (int i = pos; i < currentString.length(); i++)
	{
		std::string n = currentString;
		char originalC = n[i];
		for (int j = 0; j < alphabet.length(); j++)
		{
			char newC = alphabet[j];
			if (newC == originalC)
				continue;
			int scoreDifference = substitutionMatrix[CharacterIndex(originalC)][CharacterIndex(newC)] - 1;
			int newScore = score + scoreDifference;
			if (newScore >= T) //also for = FindHSWs because octaves are also score difference 0!
			{
				n[i] = newC;
				index[n].push_back(HSW(newScore, HSWpos));
				FindHSWs(n, newScore, T, i + 1, index, alphabet, HSWpos, substitutionMatrix);
			}
		}
	}
}

std::unordered_map<std::string, std::vector<HSW>> BLAST::GenerateHSWIndex(std::string query)
{
	//Generate high scoring words & insert into index
	std::unordered_map<std::string, std::vector<HSW>> index;
	std::vector<bool> p(q, false);
	for (int i = 0; i < query.length() - q + 1; i++)
	{
		std::string qgram = query.substr(i, q);
		int score = qgram.length();
		index[qgram].push_back(HSW(score, i));
		FindHSWs(qgram, score, T, 0, index, db->alphabet, i, substitutionMatrix);
	}
	return index;
}

void BLAST::UngappedExtension(int qpos, int epos, std::string query, std::string entry, int& score, int& ql, int& qr, int& el, int& er)
{
	//Extend to left
	ql = qpos;
	el = epos;
	int best = T;
	while (ql > 0 && el > 0 && (score + substitutionMatrix[CharacterIndex(query[ql - 1])][CharacterIndex(entry[el - 1])] >= (best - X)))
	{
		ql--;
		el--;
		score += substitutionMatrix[CharacterIndex(query[ql])][CharacterIndex(entry[el])];
		if (score > best)
			best = score;
	}
	//Extend to right
	qr = qpos + q - 1;
	er = epos + q - 1;
	while (qr < query.length() - 1 && er < entry.length() - 1 && (score + substitutionMatrix[CharacterIndex(query[qr + 1])][CharacterIndex(entry[er + 1])]) >= (best - X))
	{
		qr++;
		er++;
		score += substitutionMatrix[CharacterIndex(query[qr])][CharacterIndex(entry[er])];
		if (score > best)
			best = score;
	}
}

//Searching
double ProbSGreaterOrEqual(double m, int n, int S)
{
	//			K									Y						//TODO: calculate parameters in code
	double y = 0.27231929561537233 * m * n * exp(-1.0499992370605469 * S);	//Parameters K and Y are specifically calculated for EMO database
	return 1 - exp(-y);
}

std::vector<Result> BLAST::SearchSequence(std::string query)
{
#pragma region Initialization
	std::vector<Result> result(db->db.size());
	std::clock_t start = std::clock();
	//Generate high scoring word index from query
	std::unordered_map<std::string, std::vector<HSW>> indexHSW = GenerateHSWIndex(query);
	double duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";;"; //Indexing + empty
#pragma endregion

#pragma region Searching
	start = std::clock();
	//Scan database
#pragma omp parallel for
	for (int i = 0; i < db->db.size(); i++)
	{
		Entry dbentry = db->db[i];
		std::string entry = dbentry.sequence;
		std::vector<int> diagonals(entry.size() + query.size() - 1, -1);

		double highest = 0;
		bool significant = false;
		//For every word in entry search index
		for (int j = 0; j < entry.length() - q + 1; j++)
		{
			std::string qgram = entry.substr(j, q);
			std::vector<HSW> HSWs;
			if (indexHSW.find(qgram) != indexHSW.end())
				HSWs = indexHSW[qgram]; //High scoring words
			for (int x = 0; x < HSWs.size(); x++)
			{
				//Check if two non-overlapping hits on same diagonal within distance A
				int qpos = HSWs[x].pos;
				int epos = j;
				int previous = diagonals[qpos - epos + entry.size()];
				if (qpos > previous + q - 1) //Found hit doesn't overlap with previous one.
				{
					diagonals[qpos - epos + entry.size()] = qpos; //Update most recent found hit
					if (previous != -1 && qpos - previous <= A) //Two non-overlappings hits found within distance A -> ungapped extension
					{
						int score = HSWs[x].score;
						int ql; int qr; int el; int er;
						UngappedExtension(qpos, epos, query, entry, score, ql, qr, el, er);

						double p = ProbSGreaterOrEqual(db->sumLengths, query.length(), score);
						if (p < 0.01)
							significant = true;
						if (score > highest)
							highest = score;
					}
				}
			}
		}
		result[i] = Result(dbentry.index, highest, significant);
	}
	duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";"; //Searching
#pragma endregion

#pragma region Sorting results
	start = std::clock();
	//Sort results && write to file
	std::sort(result.begin(), result.end());
	duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";\n"; //Sorting
#pragma endregion
	return result;
}