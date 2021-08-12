/*
	Written by Jelle Mulyadi, 2021
*/

#include "SSMAW.h"
#include "Suffix.h"
#include "Stack.h"
#include "../General/Result.h"
#include "../General/Tools.h"
#include <set>
#include <iostream>
#include <ctime>
#include <vector>
#include <algorithm>
#include <bitset>

//Algorithms and code based on:
// - Barton, Carl, et al. "Linear-time computation of minimal absent words using suffix array." BMC bioinformatics 15.1 (2014): 1-10.
// - Crawford, Tim, Golnaz Badkobeh, and David Lewis. "Searching page-images of early music scanned with OMR: A scalable solution using minimal absent words." (2018): 233-239.

SSMAW::SSMAW(Database* db, int min, int max, int nrOfResults)
{
	this->db = db;
	this->min = min;
	this->max = max;
	this->nrOfResults = nrOfResults;

	//Indexing
	std::clock_t start = std::clock();
	Indexing();
	double duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	indexTime = (int)duration;
	std::cout << "Time elapsed indexing: " << duration << "\n";
}

//Indexing
int LongestCommonPrefix(std::string entry, int startA, int startB)
{
	if (startA < startB)
	{
		int temp = startA;
		startA = startB;
		startB = temp;
	}
	int result = 0;
	for (int i = 0; i < entry.length() - startA; i++)
	{
		if (entry[startA + i] == entry[startB + i])
			result++;
		else
			break;
	}
	return result;
}

void TopDownPass(std::string entry, std::vector<int>& SA, std::vector<int>& LCP, std::vector<std::bitset<53>>& B1, std::vector<std::bitset<53>>& B2)
{
	//Initialize interval array: saves all left neighbors of all factors sized <= max LCP
	std::vector<std::bitset<53>> interval;
	for (int i = 0; i <= *std::max_element(LCP.begin(), LCP.end()); i++)
		interval.push_back(std::bitset<53>());
	Stack LIFOLCP = Stack();
	LIFOLCP.Push(0);
	for (int i = 0; i < entry.size(); i++)
	{
		if (i > 0 && LCP[i] < LCP[i - 1]) //Prefix is unrelated to some of the previous prefixes -> reset interval array for those lengths
		{
			std::bitset<53> save; //Contains the left neighbors of the LCP above current LCP
			while (LIFOLCP.Count() != 0 && LIFOLCP.Top()->value > LCP[i]) //Reset all intervals of unrelated prefixes = previously encountered bigger prefixes
			{
				int length = LIFOLCP.Pop();
				save = interval[length];
				interval[length] = std::bitset<53>();
			}
			if (LIFOLCP.Top()->value < LCP[i]) //Current LCP is bigger then top (not equal) -> interval doesn't have values yet -> set interval to last found interval
				interval[LCP[i]] = save;
			B1[2 * i - 1] = save; //Complete factor -> set to last found interval
			B2[2 * i - 1] = interval[LCP[i]]; //Longest prefix factor, same prefix as current -> set to previously found interval
		}
		if (SA[i] > 0) //Index > 0 -> add left neighbor (Index 0 doesn't have left neighbour)
		{
			int ln = CharacterIndex(entry[SA[i] - 1]); //Left neighbour of suffix array
			Node* current = LIFOLCP.Top();
			while (current != nullptr && !interval[current->value][ln]) //Add left neighbor to all earlier encountered prefixes of suffix too
			{															//When encounter already set to one, all smaller also already set, so stop
				interval[current->value][ln] = 1;
				current = current->next;
			}
			interval[LCP[i]][ln] = 1; //Add left neighbor to interval of current length & B arrays
			B1[2 * i][ln] = 1;
			B2[2 * i][ln] = 1; //TODO: Can be removed?
			B1[2 * i + 1][ln] = 1;
			B2[2 * i + 1][ln] = 1;
		}
		if (i > 0 && LCP[i] > 0 && SA[i - 1] > 0) //LCP with previous suffix array is > 0 and previous suffix array is not index 0 -> add its left neighbor to interval too
		{
			int ln = CharacterIndex(entry[SA[i - 1] - 1]);
			interval[LCP[i]][ln] = 1;
		}
		B2[2 * i] = interval[LCP[i]]; //Add previously known left neighbors
		if (LIFOLCP.Top()->value != LCP[i]) //New LCP -> push to stack
			LIFOLCP.Push(LCP[i]);
	}
}

void BottomUpPass(std::string entry, std::vector<int>& LCP, std::vector<std::bitset<53>>& B1, std::vector<std::bitset<53>>& B2)
{
	//Initialize interval array: saves all left neighbors of all factors sized <= max LCP
	std::vector<std::bitset<53>> interval;
	for (int i = 0; i <= (*std::max_element(LCP.begin(), LCP.end())) + 1; i++)
		interval.push_back(std::bitset<53>());
	Stack LIFOLCP = Stack();
	Stack LIFORem = Stack();
	LIFOLCP.Push(0);
	int saveB;
	for (int i = entry.size() - 1; i >= 0; i--)
	{
		int saveA = LCP[i] + 1; //saveA is the lcp that is one above LCP[i]
		if (i < entry.size() - 1 && LCP[i] < LCP[i + 1]) //If current LCP is smaller then previous -> save all higher LCP values to be removed later
		{
			while (LIFOLCP.Count() != 0 && LIFOLCP.Top()->value > LCP[i]) //Only smaller/equal LCP remain
			{
				saveA = LIFOLCP.Pop();
				LIFORem.Push(saveA); //Larger are saved
			}
			if (LIFOLCP.Top()->value < LCP[i]) //Current LCP is bigger then top (not equal) -> interval doesn't have values yet -> set interval to last found interval
				interval[LCP[i]] = interval[saveA];
		}
		for (int ln = 0; ln < 53; ln++) //For every character in alphabet check if it is left neighbor of complete factor
		{							//If set, add to all previous unset intervals < current & current interval
			if (B1[2 * i][ln])
			{
				Node* current = LIFOLCP.Top();
				while (current != nullptr && !interval[current->value][ln]) //When encounter already set to one, all smaller also already set, so stop
				{
					interval[current->value][ln] = 1;
					current = current->next;
				}
				interval[LCP[i]][ln] = 1;
			}
		}
		B2[2 * i] = B2[2 * i] | interval[LCP[i]]; //Current prefix -> current LCP[i]
		if (i < entry.size() - 1)
		{
			B2[2 * i + 1] = B2[2 * i + 1] | interval[LCP[i + 1]]; //Current prefix -> Previous LCP[i + 1]
			B1[2 * i + 1] = B1[2 * i + 1] | interval[saveB]; //Current complete factor -> saveB
		}													//saveB is the lcp that is one above LCP[i + 1], so >= Current complete factor

		saveB = saveA; //Save for next loop (current becomes previous in next loop)
		if (i < entry.size() - 1 && LCP[i] < LCP[i + 1])
		{
			B1[2 * i] = B1[2 * i] | interval[saveA]; //saveA is one above LCP[i], so >= complete factor
			while (LIFORem.Count() != 0) //Remove all saved LCPs
			{
				int length = LIFORem.Pop();
				interval[length] = std::bitset<53>();
			}
		}
		if (LIFOLCP.Top()->value != LCP[i]) //New LCP -> push to stack
			LIFOLCP.Push(LCP[i]);
	}
}

void CalculateArrays(std::string entry, std::vector<int>& SA, std::vector<int>& LCP, std::vector<std::bitset<53>>& B1, std::vector<std::bitset<53>>& B2)
{
	std::vector<Suffix> suffixes;
	//Find suffixes
	for (int i = 0; i < entry.size(); i++)
		suffixes.push_back(Suffix(entry.substr(i), i));
	//Create SA
	std::sort(suffixes.begin(), suffixes.end());
	for (int i = 0; i < entry.size(); i++)
		SA.push_back(suffixes[i].pos);
	//Create LCP -> LCP[i] = lcp(SA[i-1], SA[i])
	LCP.push_back(0);
	for (int i = 1; i < SA.size(); i++)
	{
		int startA = SA[i - 1];
		int startB = SA[i];
		LCP.push_back(LongestCommonPrefix(entry, startA, startB));
	}
	//Calculate B1 and B2 (B1 = left neighbors factor, B2 = left neigbors its longest proper prefix
	TopDownPass(entry, SA, LCP, B1, B2); //Top down pass through arrays
	BottomUpPass(entry, LCP, B1, B2); //Bottom up pass through arrays
}

void SSMAW::CalculateMAWs(std::string entry, std::vector<int>& SA, std::vector<int>& LCP,
	std::vector<std::bitset<53>>& B1, std::vector<std::bitset<53>>& B2, std::set<std::string>& MAWs)
{
	for (int j = 0; j < entry.size() * 2 - 1; j++)
	{
		std::bitset<53> difference = B1[j] ^ B2[j]; //Difference between set of neigbors of factor and longest proper prefix of that factor gives all the MAWs
		int index = j / 2;
		int plus = (j % 2);
		if ((SA[index] + LCP[index + plus]) > entry.size() - 1) //Index out of bounds -> continue
			continue;
		for (int z = 0; z < difference.size(); z++)
		{
			if (difference[z])
			{
				char l = db->alphabet[z];
				int length = LCP[index + plus] + 2;
				//Take MAWs of sizes between min and max
				if (min <= length && length <= max)
				{
					std::string MAW = l + entry.substr(SA[index], LCP[index + plus] + 1);
					MAWs.insert(MAW);
				}
			}
		}
	}
}

void SSMAW::Indexing()
{
	MAWcounts = std::vector<int>(db->db.size(), 0);
#pragma omp parallel for
	for (int i = 0; i < db->db.size(); i++)
	{
		std::string entry = db->db[i].sequence;
		//Calculate all arrays
		std::vector<int> SA; //Suffix array
		std::vector<int> LCP; //Longest common prefix array
		std::vector<std::bitset<53>> B1(entry.size() * 2, std::bitset<53>()); //TODO: allow for variable alphabet size
		std::vector<std::bitset<53>> B2(entry.size() * 2, std::bitset<53>()); //TODO: allow for variable alphabet size
		CalculateArrays(entry, SA, LCP, B1, B2);
		//Calculate MAWs
		std::set<std::string> MAWs;
		CalculateMAWs(entry, SA, LCP, B1, B2, MAWs);
		//Save number of MAWs for each database entry
		MAWcounts[i] = MAWs.size();
		//Build Trie for MAWs for quick search
#pragma omp critical
		{
			for (std::string w : MAWs)
				MAWsTrie.Insert(w, i);
		}
	}
}

//Searching
std::vector<Result> SSMAW::SearchSequence(std::string query)
{
#pragma region Initialization
	std::clock_t start = std::clock();
	//Calculate all arrays
	std::vector<int> SA; //Suffix array
	std::vector<int> LCP; //Longest common prefix array
	std::vector<std::bitset<53>> B1(query.size() * 2, std::bitset<53>()); //TODO: allow for variable alphabet size
	std::vector<std::bitset<53>> B2(query.size() * 2, std::bitset<53>()); //TODO: allow for variable alphabet size
	CalculateArrays(query, SA, LCP, B1, B2);
	//Calculate MAWs
	std::set<std::string> MAWs;
	CalculateMAWs(query, SA, LCP, B1, B2, MAWs);
	double duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";"; //Indexing query
#pragma endregion

#pragma region Searching
	start = std::clock();
	//Calculate scores = number of common MAWs between query and database entry
	auto it = MAWs.begin();
	std::vector<int> score = std::vector<int>(db->db.size(), 0);
	while (it != MAWs.end())
	{
		std::vector<int> indicesEntries = MAWsTrie.Search((*it));
		for (int i = 0; i < indicesEntries.size(); i++)
		{
			score[indicesEntries[i]]++;
		}
		it++;
	}
	duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";"; //Scoring

	start = std::clock();
	//Calculate results -> jaccard distance = 1-(intersect/union)
	std::vector<Result> result(db->db.size());
//#pragma omp parallel for
	for (int i = 0; i < score.size(); i++)
	{
		double jaccard = (double)1 - ((double)score[i] / (double)(MAWcounts[i] + MAWs.size() - score[i]));
		result[i] = Result(db->rIndex[i], jaccard, true); //TODO: Fix rIndex! If PassJoin sorts, rIndex is not correct (right now SSMAW first in testing order, so not a problem)
	}
	duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
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