/*
	Written by Jelle Mulyadi, 2021
*/

#include "Trie.h"
struct TrieNode* makeNode()
{
	struct TrieNode* node = new TrieNode;
	for (int i = 0; i < 53; i++)
	{
		node->children[i] = NULL;
	}
	return node;
}

void TrieDestructor(TrieNode* node)
{
	if (node != nullptr)
	{
		int count = 0;
		for (int i = 0; i < 53; i++)
			if (node->children[i] != NULL)
				count++;
		if (count == 0)
		{
			delete(node);
			return;
		}

		for (TrieNode* child : node->children)
		{
			if (child != nullptr)
				TrieDestructor(child);
		}
		delete(node);
		return;
	}
}