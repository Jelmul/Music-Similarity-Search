/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include <string>

struct Entry
{
	std::string index;
	std::string sequence;
	Entry() {};
	Entry(std::string index, std::string sequence)
	{
		this->index = index;
		this->sequence = sequence;
	}

	bool operator<(const Entry& rhs) const noexcept
	{
		return this->sequence < rhs.sequence;
	}
};