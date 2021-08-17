/*
	Written by Jelle Mulyadi, 2021
*/

#include "PassJoin.h"
#include "..\General\Tools.h"
#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>

//Algorithms and code based on:
// - Li, Guoliang, et al. "Pass-join: A partition-based method for similarity joins." arXiv preprint arXiv:1111.7171 (2011).
// - Qin, Jianbin, and Chuan Xiao. "Pigeonring: A principle for faster thresholded similarity search." arXiv preprint arXiv:1804.01614 (2018).

PassJoin::PassJoin(Database* db, int threshold, int chainLength)
{
	this->db = db;
	this->threshold = threshold;
	this->chainLength = chainLength;

	//Indexing
	std::clock_t start = std::clock();
	Indexing();
	double duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	indexTime = (int)duration;
	std::cout << "Time elapsed indexing: " << duration << "\n";
}

//Indexing
void PassJoin::Indexing()
{
	std::sort(db->db.begin(), db->db.end(), CompareLength); //Sort on length
	//Sort strings of same length on alphabetical order
	auto it = db->db.begin();
	auto it2 = db->db.begin() + 1;
	while (it != db->db.end())
	{
		while (it2 != db->db.end() && (*it).sequence.length() == (*it2).sequence.length())
			std::advance(it2, 1);
		std::sort(it, it2);
		it = it2;
	}
	//Iterate over database & build index
	int nrOfSegments = threshold + 1;
	for (int i = 0; i < db->db.size(); i++)
	{
		std::string entry = db->db[i].sequence;
		int f = floor((double)entry.length() / nrOfSegments);
		int c = ceil((double)entry.length() / nrOfSegments);
		int k = entry.length() - f * nrOfSegments;
		int pos = 0;
		int segLength;
		if (index.find(entry.length()) == index.end()) //Initialize index if current length doesn't exist yet
			for (int j = 0; j < nrOfSegments; j++)
				index[entry.length()].push_back(std::unordered_map<std::string, std::vector<int>>());
		for (int j = 0; j < nrOfSegments; j++)
		{
			if (j < nrOfSegments - k)
				segLength = f;
			else
				segLength = c;
			index[entry.length()][j][entry.substr(pos, segLength)].push_back(i);
			pos += segLength;
		}
	}
}

//Searching
bool PassJoin::PigeonRing(std::string candidate, int segmentNr, int segmentPos, std::string query, std::vector<double>& thresholds)
{
	if (chainLength == 0) //Chain length 0 = off, chain length threshold + 1 = equal to alignment filter
		return true;

	int nrOfSegments = threshold + 1;
	int f = floor((double)candidate.length() / nrOfSegments);
	int c = ceil((double)candidate.length() / nrOfSegments);
	int k = candidate.length() - f * nrOfSegments;
	int pos = segmentPos;
	int nr = segmentNr + 1;
	//Find prefix-viable chain
	double errors = 0; //double because hamming distance returns double
	for (int j = 1; j < chainLength; j++)
	{
		pos = pos % candidate.length();
		int index = nr % nrOfSegments;
		int segLength;
		if (index < nrOfSegments - k)
			segLength = f;
		else
			segLength = c;
		int startQ = std::max(0, pos - threshold);
		int lengthQ = std::min((int)query.length(), (pos + segLength + threshold)) - startQ;
		int startC = pos;
		int lengthC = segLength;
		errors += SubstringEditDistance(query, candidate, startQ, lengthQ, startC, lengthC);
		if (errors > thresholds[j - 1])
			return false;
		pos += segLength;
		nr++;
	}
	//Prefix-viable chain found -> pigeonring filter passed
	return true;
}

bool PassJoin::AlignmentFilter(std::string candidate, std::string query)
{
	int nrOfSegments = threshold + 1;
	int f = floor((double)candidate.length() / nrOfSegments);
	int c = ceil((double)candidate.length() / nrOfSegments);
	int k = candidate.length() - f * nrOfSegments;
	int pos = 0;
	std::string segment;
	int errors = 0;
	for (int i = 0; i < nrOfSegments; i++)
	{
		int segLength;
		if (i < nrOfSegments - k)
			segLength = f;
		else
			segLength = c;
		int startQ = std::max(0, pos - threshold);
		int lengthQ = std::min((int)query.length(), (pos + segLength + threshold)) - startQ;
		int startC = pos;
		int lengthC = segLength;
		errors += SubstringEditDistance(query, candidate, startQ, lengthQ, startC, lengthC);
		if (errors > threshold)
			return false;
		pos += segLength;
	}
	//Alignment filter passed
	return true;
}

void PassJoin::Verification(std::string query, std::vector<int>& candidates, int entryPos, int entrySeg, int segLength, std::vector<double>& thresholds, std::vector<Result>& result, std::vector<bool>& checked)
{
	for (int i = 0; i < candidates.size(); i++)
	{
		if (!checked[candidates[i]]) //Only verify entries that haven't been verified yet
		{
			Entry dbCandidate = db->db[candidates[i]];
			std::string candidate = dbCandidate.sequence;
			if (PigeonRing(candidate, entrySeg, entryPos + segLength, query, thresholds)) //Pigeonring filter -> try to find a prefix viable chain
			{
				int score = LengthAwareED(query, candidate, 0, query.length(), 0, candidate.length(), threshold); //Verify candidate, using expensive ED calculation
				if (score <= threshold)
					result.push_back(Result(dbCandidate.index, score, true));
				checked[candidates[i]] = true;
			}
		}
	}
}

void PassJoin::SubstringSelection(std::string query, int pos, int segment, int segLength, int entryLength, int& start, int& end)
{
	//Multi-match aware method
	int delta = abs((int)query.length() - entryLength);
	int minL = std::max(0, pos - segment);
	int maxL = std::min((int)query.length() - segLength, pos + segment);
	int minR = std::max(0, pos + delta - (threshold - segment));
	int maxR = std::min((int)query.length() - segLength, pos + delta + (threshold - segment));
	start = std::max(minL, minR);
	end = std::min(maxL, maxR);
}

std::vector<Result> PassJoin::SearchSequence(std::string query)
{
#pragma region Initialization
	std::clock_t startC = std::clock();
	//Variables
	std::vector <Result> result;
	std::vector<bool> checked = std::vector<bool>(db->db.size(), false); //Inserted entries
	int minL = (*index.begin()).first; //Minimum length
	int maxL = (*index.rbegin()).first; //Maximum length
	//Pre calculate thresholds for pigeonring
	int nrOfSegments = threshold + 1;
	double single = (double)threshold / (double)nrOfSegments;
	std::vector<double> thresholds;
	for (int i = 1; i < chainLength; i++)
		thresholds.push_back((double)(i + 1) * single);
	double duration = (std::clock() - startC) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";"; //Initialization
#pragma endregion

#pragma region Searching
	startC = std::clock();
	auto startLength = index.lower_bound(std::max((int)query.length() - threshold, minL));
	auto endLength = index.upper_bound(std::min((int)query.length() + threshold, maxL));
	for (auto lengthIT = startLength; lengthIT != endLength; lengthIT++) //Length-based filter: only consider entries with possible lengths
	{
		int pos = 0;
		for (int i = 0; i < nrOfSegments; i++) //For every segment find matches -> candidates (Pigeonhole filter)
		{
			int currentLength = (*lengthIT).first;
			int segLength = (*(*lengthIT).second[i].begin()).first.length();
			int start; int end;
			SubstringSelection(query, pos, i, segLength, currentLength, start, end); //Select substrings for pigeonhole filter
			for (int j = start; j <= end; j++)
			{
				std::string substring = query.substr(j, segLength);
				auto substringIT = (*lengthIT).second[i].find(substring); //Look for exact matches of substring
				if (substringIT != (*lengthIT).second[i].end()) //Pigeonhole filter: if there is an exact match, entry is a candidate
					Verification(query, (*substringIT).second, pos, i, segLength, thresholds, result, checked); //For candidates: verify
			}
			pos += segLength;
		}
	}
	duration = (std::clock() - startC) / (CLOCKS_PER_SEC / 1000);
	std::cout << ";" << duration << ";"; //Searching
#pragma endregion

#pragma region Sorting results
	startC = std::clock();
	std::sort(result.rbegin(), result.rend());
	duration = (std::clock() - startC) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";\n"; //Sorting
#pragma endregion
	return result;
}
