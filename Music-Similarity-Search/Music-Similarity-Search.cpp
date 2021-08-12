/*
	Written by Jelle Mulyadi, 2021
*/

#include "Music-Similarity-Search.h"
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <ctime>
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "windows.h"
#include "psapi.h"
//TODO: add k to algorithms

//VARIABLES
int k = 100;
int EDthreshold;
//MAW
int MAWmin = 4;
int MAWmax = 8;
//BLAST
int BLT = 4;
int BLA = 10;
int BLX = 4;
int BLq = 4;
//Pass-Join
int PASSchain;
//PIVOTAL
int PIVq;
int PIVchain;

//FUNCTIONS
struct GroundTruth
{
	std::string query;
	std::string result;
	GroundTruth() {};
	GroundTruth(std::string query, std::string result)
	{
		this->query = query;
		this->result = result;
	}
};

std::vector<GroundTruth> QueryList(std::string name)
{
	std::vector<GroundTruth> list;
	std::fstream file;
	file.open(name, std::ios::in);
	if (file.is_open())
	{
		std::string line;
		while (std::getline(file, line))
		{
			int j = 0;
			while (line[j] != '\t')
				j++;
			int j2 = j + 1;
			while (line[j2] != '\t')
				j2++;
			j2 = j2 - (j + 1);
			std::string query = line.substr(0, j);
			std::string result = line.substr(j + 1, j2);
			list.push_back(GroundTruth(query, result));
		}
	}
	return list;
}

bool getBool(char c)
{
	if (c == '0')
		return false;
	return true;
}

double StandardDeviation(std::vector<double> array, double mean, int nr)
{
	double variance = 0;
	for (int n = 0; n < nr; n++)
	{
		variance += ((array[n] - mean) * (array[n] - mean));
	}
	variance /= nr;
	return sqrt(variance);
}

//RETRIEVAL EXPERIMENTS
void QueryRetrieval(SimilaritySearch* ss, std::string algorithm, std::vector<std::vector<GroundTruth>>& querylists)
{
	for (int i = 0; i < querylists.size(); i++)
	{
		//Initialization
		std::string task;
		switch (i) {
		case 0:
			task = "duplicates";
			break;
		case 1:
			task = "same";
			break;
		case 2:
			task = "relevant";
			break;
		}
		std::vector<GroundTruth> list = querylists[i];
		std::vector<double> ranks;
		std::vector<double> scores;
		int sumRank = 0;
		int minRank = INT_MAX;
		int maxRank = INT_MIN;
		double sumScore = 0;
		double minScore = DBL_MAX;
		double maxScore = DBL_MIN;

		//Retrieval
		std::clock_t start = std::clock();
		int nrOfQueries = list.size();
		std::ofstream file;
		file.open("Retrieval performance " + algorithm + " " + task + ".csv");
		if (file.is_open())
		{
			for (int j = 0; j < nrOfQueries; j++)
			{
				GroundTruth current = list[j];
				std::string query = current.query;
				std::string truth = current.result;
				std::vector<Result> result = ss->SearchSequenceID(query);
				auto it = std::find(result.begin(), result.end(), truth); //Find truth in result list
				int rank = std::distance(result.begin(), it) + 1;
				if (rank < minRank)
					minRank = rank;
				if (rank > maxRank)
					maxRank = rank;
				double score = (*it).score;
				if (score < minScore)
					minScore = score;
				if (score > maxScore)
					maxScore = score;
				sumRank += rank;
				sumScore += score;
				ranks.push_back(rank);
				scores.push_back(score);
				file << query << ";" << truth << ";" << rank << ";" << score << ";\n";
			}

			//Statistics query time
			double duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
			double avgQueryTime = duration / nrOfQueries;
			file << "Mean query time;\n";
			file << avgQueryTime << ";\n";
			//Statistics rank
			double avgRank = sumRank / nrOfQueries;
			double stdevRank = StandardDeviation(ranks, avgRank, nrOfQueries);
			file << "Mean rank; Min rank; Max rank; Standard deviation rank;\n";
			file << avgRank << ";" << minRank << ";" << maxRank << ";" << stdevRank << ";\n";
			//Statistics score
			double avgScore = sumScore / nrOfQueries;
			double stdevScore = StandardDeviation(scores, avgScore, nrOfQueries);
			file << "Mean score; Min score; Max score; Standard deviation score;\n";
			file << avgScore << ";" << minScore << ";" << maxScore << ";" << stdevScore << ";\n";
		}
		file.close();
	}
}

void RetrievalPerformance(std::string algorithmBools, std::string lists)
{
	Database* database = new Database("emo.txt");
	std::vector<std::vector<GroundTruth>> querylists;
	if (getBool(lists[0]))
		querylists.push_back(QueryList("duplicates.txt"));
	if (getBool(lists[1]))
		querylists.push_back(QueryList("same_music.txt"));
	if (getBool(lists[2]))
		querylists.push_back(QueryList("relevant.txt"));

	if (getBool(algorithmBools[0]))
	{
		std::cout << "Starting SSMAW, database size= " << database->db.size() << "\n";
		SSMAW* MinimalAbsentWordSS = new SSMAW(database, MAWmin, MAWmax, k); //min, max
		QueryRetrieval(MinimalAbsentWordSS, "MAW", querylists);
		delete MinimalAbsentWordSS;
	}

	if (getBool(algorithmBools[1]))
	{
		std::cout << "Starting ED, database size= " << database->db.size() << "\n";
		ED* EditDistanceSS = new ED(database, k);
		QueryRetrieval(EditDistanceSS, "ED", querylists);
		delete EditDistanceSS;
	}

	if (getBool(algorithmBools[2]))
	{
		std::cout << "Starting GA, database size= " << database->db.size() << "\n";
		GA* GlobalAlignmentSS = new GA(database, k);
		QueryRetrieval(GlobalAlignmentSS, "GA", querylists);
		delete GlobalAlignmentSS;
	}

	if (getBool(algorithmBools[3]))
	{
		std::cout << "Starting LA, database size= " << database->db.size() << "\n";
		LA* LocalAlignmentSS = new LA(database, k);
		QueryRetrieval(LocalAlignmentSS, "LA", querylists);
		delete LocalAlignmentSS;
	}

	if (getBool(algorithmBools[4]))
	{
		std::cout << "Starting BLAST, database size= " << database->db.size() << "\n";
		BLAST* BasicLocalAlignmentSS = new BLAST(database, BLT, BLA, BLX, BLq); //T, A, X, q
		QueryRetrieval(BasicLocalAlignmentSS, "BLAST", querylists);
		delete BasicLocalAlignmentSS;
	}

	if (getBool(algorithmBools[5]))
	{
		std::cout << "Starting PassJoin, database size= " << database->db.size() << "\n";
		PassJoin* PassJoinSS = new PassJoin(database, EDthreshold, PASSchain);
		QueryRetrieval(PassJoinSS, "PassJoin", querylists);
		delete PassJoinSS;
	}

	if (getBool(algorithmBools[6]))
	{
		std::cout << "Starting PIVOTAL, database size= " << database->db.size() << "\n";
		PivotalSearch* PivotalSS = new PivotalSearch(database, PIVq, EDthreshold, PIVchain);
		QueryRetrieval(PivotalSS, "PIVOTAL", querylists);
		delete PivotalSS;
	}
	delete database;
}

//SCALABILITY EXPERIMENTS
void ScalabilityTest(SimilaritySearch* ss, std::vector<GroundTruth>& queryList, std::string& line)
{
	std::clock_t start = std::clock();

	for (int i = 0; i < queryList.size(); i++)
		ss->SearchSequenceID(queryList[i].query);
	double duration = ((std::clock() - start) / (CLOCKS_PER_SEC / 1000)) / queryList.size();
	line += std::to_string((int)duration) + ";";
}

void Scalability(std::string algorithmBools, std::string queryFile)
{
	std::ofstream file;
	file.open("Scalability.csv");
	if (file.is_open())
	{
		std::clock_t startToFinish = std::clock();
		double duration;
		file << "Number of sequences;";
		std::vector<std::string> labels{ "SSMAW IndexSize;SSMAW;SSMAW IndexTime", "ED", "GA", "LA", "BLAST", "PassJoin IndexSize;PassJoin;PassJoin IndexTime", "Pivotal IndexSize;Pivotal;Pivotal IndexTime" };
		for (int i = 0; i < labels.size(); i++)
		{
			if (getBool(algorithmBools[i]))
				file << labels[i] << ";";
		}
		file << "\n";

		std::vector<GroundTruth> queryList = QueryList(queryFile);
		for (int i = 0; i < 10; i++)
		{
			int nr = i + 1;
			Database* database = new Database("emo_" + std::to_string(nr) + "a.txt");
			std::string line = std::to_string(database->db.size()) + ";";

			if (getBool(algorithmBools[0]))
			{
				std::cout << "Starting SSMAW, database size= " << database->db.size() << "\n";
				SSMAW* MinimalAbsentWordSS = new SSMAW(database, MAWmin, MAWmax, k); //min, max
				PROCESS_MEMORY_COUNTERS_EX pmc;
				GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
				line += std::to_string(pmc.PrivateUsage / 1000000) + ";";
				ScalabilityTest(MinimalAbsentWordSS, queryList, line);
				line += std::to_string(MinimalAbsentWordSS->indexTime) + ";";
				delete MinimalAbsentWordSS;
			}

			if (getBool(algorithmBools[1]))
			{
				std::cout << "Starting ED, database size= " << database->db.size() << "\n";
				ED* EditDistanceSS = new ED(database, k);
				ScalabilityTest(EditDistanceSS, queryList, line);
				delete EditDistanceSS;
			}

			if (getBool(algorithmBools[2]))
			{
				std::cout << "Starting GA, database size= " << database->db.size() << "\n";
				GA* GlobalAlignmentSS = new GA(database, k);
				ScalabilityTest(GlobalAlignmentSS, queryList, line);
				delete GlobalAlignmentSS;
			}

			if (getBool(algorithmBools[3]))
			{
				std::cout << "Starting LA, database size= " << database->db.size() << "\n";
				LA* LocalAlignmentSS = new LA(database, k);
				ScalabilityTest(LocalAlignmentSS, queryList, line);
				delete LocalAlignmentSS;
			}

			if (getBool(algorithmBools[4]))
			{
				std::cout << "Starting BLAST, database size= " << database->db.size() << "\n";
				BLAST* BasicLocalAlignmentSS = new BLAST(database, BLT, BLA, BLT, BLq); //T, A, X, q
				ScalabilityTest(BasicLocalAlignmentSS, queryList, line);
				delete BasicLocalAlignmentSS;
			}

			if (getBool(algorithmBools[5]))
			{
				std::cout << "Starting PassJoin, database size= " << database->db.size() << "\n";
				PassJoin* PassJoinSS = new PassJoin(database, EDthreshold, PASSchain);
				PROCESS_MEMORY_COUNTERS_EX pmc;
				GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
				line += std::to_string(pmc.PrivateUsage / 1000000) + ";";
				ScalabilityTest(PassJoinSS, queryList, line);
				line += std::to_string(PassJoinSS->indexTime) + ";";
				delete PassJoinSS;
			}

			if (getBool(algorithmBools[6]))
			{
				std::cout << "Starting PIVOTAL, database size= " << database->db.size() << "\n";
				PivotalSearch* PivotalSS = new PivotalSearch(database, PIVq, EDthreshold, PIVchain);
				PROCESS_MEMORY_COUNTERS_EX pmc;
				GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
				line += std::to_string(pmc.PrivateUsage / 1000000) + ";";
				ScalabilityTest(PivotalSS, queryList, line);
				line += std::to_string(PivotalSS->indexTime) + ";";
				delete PivotalSS;
			}

			file << line << "\n";
			delete database;
		}
		duration = (std::clock() - startToFinish) / (CLOCKS_PER_SEC / 1000);
		file << "Total time elapsed: " + std::to_string((int)duration) + "\n";
	}
	file.close();
}

//MAIN
void Experiments()
{
	int mode;
	std::cout << "Choose experiment music similarity search:\n";
	std::cout << "0: Retrieval Performance, 1: Scalability, 2: Both\n";
	std::cin >> mode;
	std::string algorithms;
	std::cout << "Choose which algorithms to test:\n";
	std::cout << "Enter string of 7 zeroes and ones (e.g. 1111111 for all algorithms)\n";
	std::cout << "In order: SSMAW, ED, GA, LA, BLAST, PassJoin, PIVOTAL\n";
	std::cin >> algorithms;
	if (getBool(algorithms[5]) || getBool(algorithms[6]))
	{
		std::cout << "Enter edit distance threshold:\n";
		std::cin >> EDthreshold;
	}
	if (getBool(algorithms[5]))
	{
		std::cout << "Enter PassJoin chain length:\n";
		std::cin >> PASSchain;
	}
	if (getBool(algorithms[6]))
	{
		std::cout << "Enter PIVOTAL q size:\n";
		std::cin >> PIVq;
		std::cout << "Enter PIVOTAL chain length:\n";
		std::cin >> PIVchain;
	}

	//Experiment 1 & 3 - Retrieval performance
	if (mode == 0 || mode == 2)
	{
		std::string lists;
		std::cout << "Choose which query lists to test:\n";
		std::cout << "Enter string of 3 zeroes and ones (e.g. 111 for all query lists)\n";
		std::cout << "In order: dupl, same, relv\n";
		std::cin >> lists;
		RetrievalPerformance(algorithms, lists);
	}

	//Experiment 2 - Scalability
	if (mode == 1 || mode == 2)
		Scalability(algorithms, "duplicates.txt");

	std::cout << "Finished\n";
	int end = 0;
	std::cin >> end;
}

int main()
{
	Experiments();
}
