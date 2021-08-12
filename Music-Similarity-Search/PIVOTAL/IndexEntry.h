/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
struct PrefixEntry
{
	int index;
	int pos;
	PrefixEntry() {};
	PrefixEntry(int index, int pos)
	{
		this->index = index;
		this->pos = pos;
	}
};

struct PivotalEntry
{
	int index;
	int pos;
	int pivotalNr;
	PivotalEntry() {};
	PivotalEntry(int index, int pos, int pivotalNr)
	{
		this->index = index;
		this->pos = pos;
		this->pivotalNr = pivotalNr;
	}
};