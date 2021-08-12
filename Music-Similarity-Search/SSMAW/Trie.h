/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "../General/Scoring.h"
#include <bitset>
#include <vector>

//Code based on:
// - https://www.geeksforgeeks.org/trie-insert-and-search/

struct TrieNode
{
	struct TrieNode* children[53];
	std::vector<int> entries;
};

struct TrieNode* makeNode();

void TrieDestructor(TrieNode* node);

class Trie
{
private:
	//Variables
	TrieNode* root;
public:
	Trie()
	{
		this->root = makeNode();
	}
	~Trie()
	{
		TrieDestructor(root);
	}
	void Insert(std::string word, int entry)
	{
		struct TrieNode* current = root;
		for (int i = 0; i < word.length(); i++)
		{
			int index = CharacterIndex(word[i]);
			if (current->children[index] == NULL)
				current->children[index] = makeNode();
			current = current->children[index];
		}
		current->entries.push_back(entry);
	}
	std::vector<int> Search(std::string word)
	{
		TrieNode* current = root;
		for(int i = 0; i < word.length(); i++)
		{
			int index = CharacterIndex(word[i]);
			if (!current->children[index])
				return std::vector<int>();
			current = current->children[index];
		}
		if (current != NULL)
			return current->entries;
		else
			return std::vector<int>();
	}
};

