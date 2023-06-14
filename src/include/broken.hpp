#ifndef __DV_BROKEN__
#define __DV_BROKEN__


#include <vector>
#include <iostream>
#include <cmath>

#define LSB(x) ((x) & -(x))

template <typename NodeVal>
class FenwickNode {
	public:
		// Constructor
		FenwickNode() : branch(0.0), node(0.0), parent(nullptr), firstchild(nullptr), sibling(nullptr) {std::cout << "Constructor\n";};

		// Functions
		void add(const double value);

		void update(const double value);

		void zero();

		void print();

		// getters
		NodeVal & operator * ();

		double get_branch() const;

		double get_p() const;

		template <typename X>
		friend class reBrokenStick;
	private:
		NodeVal content; // stored value
		double branch; // the weight of the branch
		double node; // the weight of the node
		FenwickNode<NodeVal> * parent;
		FenwickNode<NodeVal> * firstchild;
		FenwickNode<NodeVal> * sibling;

		void _add(const double value);

		// to find largest node with with prefix_sum(i) <= value.
		FenwickNode<NodeVal> * find(double rn);

		double dummyCumsum();

};

template <typename NodeVal>
class reBrokenStick {
	public:
		reBrokenStick(unsigned int _size);


		void push_back(double p, int val); 

		// Find the largest i with prefix_sum(i) <= value.
		int draw(double rn);

		double cumsum() const;

		void print();
		    /*int range_sum(int i, int j) {
			int sum = 0;
			for (; j > i; j -= LSB(j))
				sum += A[j];
			for (; i > j; i -= LSB(i))
				sum -= A[i];
			return sum;
		    }
		    
		    // Returns the sum of the first i elements (indices 0 to i)
		    // Equivalent to range_sum(0, i)
		    int prefix_sum(int i) {
			int sum = A[0];
			for (; i != 0; i -= LSB(i))
				sum += A[i];
			return sum;
		    }*/
		FenwickNode<NodeVal> * operator[] (int i) ;
	private:
		unsigned int used; // number of nodes used
		std::vector<FenwickNode<int>> tree;

		void init(unsigned int newsize);

		void updateRelation(unsigned int index);

		void updateRelations() ;
};


#endif
