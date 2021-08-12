/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include <algorithm>
class Node
{
public:
	//Variables
	int value;
	Node* next;
	//Methods
	Node()
	{
		this->value = INT_MAX;
		next = nullptr;
	}
	Node(int value)
	{
		this->value = value;
		next = nullptr;
	}
};

class Stack
{
private:
	Node* top;
	int count = 0;
public:
	Stack()
	{
		top = nullptr;
	}

	~Stack()
	{
		Node* current = top;
		while (current != nullptr)
		{
			Node* save = current;
			current = current->next;
			delete save;
		}
	}

	Node* Top()
	{
		return top;
	}

	int Pop()
	{
		if (count != 0)
		{
			int value = top->value;
			Node* temp = top;
			top = top->next;
			delete temp;
			count--;
			return value;
		}
		else
			return INT_MIN;
	}

	void Push(int value)
	{
		Node* newTop = new Node(value);
		newTop->next = top;
		top = newTop;
		count++;
	}
	
	int Count()
	{
		return count;
	}
};