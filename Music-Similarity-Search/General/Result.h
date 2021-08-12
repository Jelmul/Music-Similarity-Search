/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include <string>
struct Result
{
	std::string s;
	double score;
	bool significant;
	Result() {}
	Result(std::string s, double score, bool significant)
	{
		this->s = s;
		this->score = score;
		this->significant = significant;
	}

	bool operator<(const Result& rhs) const noexcept
	{
		if (this->significant != rhs.significant)
			return this->significant;
		return this->score > rhs.score;
	}

	bool operator==(const std::string& rhs) const noexcept
	{
		return this->s == rhs;
	}
};