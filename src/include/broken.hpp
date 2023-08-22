#ifndef __DV_BROKEN__
#define __DV_BROKEN__


#include <vector>
#include <iostream>
#include <cmath>

#define LSB(x) ((x) & -(x))

namespace broken{

// DECLARATIONS

template <typename NodeVal>
class FenwickNode {
	public:
		// Constructor
		FenwickNode() : branch(0.0), node(0.0), parent(nullptr), firstchild(nullptr), sibling(nullptr) {}

		// Functions
		void add(const float value);

		void update(const float value);

		void zero();

		void print();

		// getters
		NodeVal & operator * ();

		float get_branch() const;

		float get_p() const;

		template <typename X>
		friend class reBrokenStick;
	private:
		NodeVal content; // stored value
		double branch; // the weight of the branch
		double node; // the weight of the node
		FenwickNode<NodeVal> * parent;
		FenwickNode<NodeVal> * firstchild;
		FenwickNode<NodeVal> * sibling;

		void _add(const float value);

		// to find largest node with with prefix_sum(i) <= value.
		FenwickNode<NodeVal> * find(double rn);

		double dummyCumsum();

};

template <typename NodeVal>
class reBrokenStick {
	public:
		//reBrokenStick(unsigned int _size);
		reBrokenStick(unsigned int _size): used(0), tree(0){
			if(_size) init(_size);
		}

		void reserve(unsigned int _size){init(_size);}

		void push_back(float p, NodeVal val); 

		inline FenwickNode<NodeVal>& back();
		// Find the largest i with prefix_sum(i) <= value.
		NodeVal draw(double rn);

		double cumsum() const;

		void print();

		//FenwickNode<NodeVal> * operator[] (int i) ;
		FenwickNode<NodeVal> & operator[] (int i) ;

	private:
		unsigned int used; // number of nodes used
		std::vector<FenwickNode<NodeVal>> tree;

		void init(unsigned int newsize);

		void updateRelation(unsigned int index);

		void updateRelations() ;
};

// DEFINITIONS

template <typename NodeVal>
void reBrokenStick<NodeVal>::init(unsigned int newsize){
	if( newsize > tree.size() ){ // if new number of values stored does not exceeds size() no need to increase space
		// calculate neccessary nodes - newsize will be the actual needed size
		if(newsize > 3) newsize= std::pow(2, static_cast<unsigned int>(std::log2(newsize-2)+1))+1;

		if( (newsize > tree.capacity()) && (!tree.empty()) ) throw std::runtime_error("Undefined\n");

		// add empty nodes
		unsigned int from = tree.size(), plus = newsize - tree.size();
		// initialise new ones
		while( plus--){
			tree.emplace_back();
		}

		// set new ones pointers
		for(;from < tree.size(); from++){
			updateRelation(from); // used is size, so used - 1 is last value, new value is used
		}
	}
}

template <typename NodeVal>
void reBrokenStick<NodeVal>::push_back(float p, NodeVal val){ 
	init(used);
	// now node exists for sure, lets assign value to it!
	tree[used].add(p);
	tree[used++].content = val;
}

template <typename NodeVal>
inline FenwickNode<NodeVal>& reBrokenStick<NodeVal>::back(){
	return tree[used-1];
}

// Find the largest i with prefix_sum(i) <= value.
template <typename NodeVal>
NodeVal reBrokenStick<NodeVal>::draw(double rn){
	if(cumsum() == 0.0) return tree[0].content;
	// init correct random value
	if(rn<0.0 || rn >1.0) throw std::invalid_argument("Invalid random number provided to reBrokenStick. It should be between 0 and 1!\n");
	rn *= cumsum();

	// check 0 position 
	if(rn < tree[0].branch) return tree[0].content; // Sometimes this is the biggest
	
	//auto out = tree.back().find(rn);
	return tree.back().find(rn)->content;
    
}

template <typename NodeVal>
double reBrokenStick<NodeVal>::cumsum() const {return tree.back().branch;}

template <typename NodeVal>
void reBrokenStick<NodeVal>::print(){
	std::cout << "ind\tnod\tbra\tcum\tval\n";
	for (int nodei = 0; nodei < used; nodei++) {
		std::cout << nodei << '\t';
		tree[nodei].print(); 
	}
}

//template <typename NodeVal>
//FenwickNode<NodeVal>* reBrokenStick<NodeVal>::operator[] (int i) {
//	return &tree[i];
//}

template <typename NodeVal>
FenwickNode<NodeVal>& reBrokenStick<NodeVal>::operator[] (int i) {
	return tree[i];
}

template <typename NodeVal>
void reBrokenStick<NodeVal>::updateRelation(unsigned int index){
	if(index == 0){
		tree[index].firstchild = nullptr;
		tree[index].sibling= nullptr;
		tree[index].parent = (tree.size()>1)?&tree[1]:nullptr;
		return;
	}

	const unsigned int my_LSB = LSB(index);

	// parent
	unsigned int target = index + my_LSB;
	tree[index].parent = (target < tree.size())?&tree[target]:nullptr;


	if(index & 1){ // there is a bit on the first position -> no child, no sibling 
		       
		// firstchild
		tree[index].firstchild = (index == 1)?&tree[0]:nullptr; // it is the first node -> its chid is the dummy 0

		// sibling
		tree[index].sibling= nullptr;
	} else { // it has at least one child and one sibling
		 
		// firstchild
		tree[index].firstchild = &tree[(index ^ my_LSB) | (my_LSB >> 1) ];

		// sibling
		target = index + (my_LSB >> 1);
		tree[index].sibling = (target < tree.size())?&tree[target]:nullptr;
	}
	//std::cout << "updated" << tree[index].branch << std::endl;
}

template <typename NodeVal>
float FenwickNode<NodeVal>::get_branch() const {return static_cast<float>(branch);}

template <typename NodeVal>
float FenwickNode<NodeVal>::get_p() const {return static_cast<float>(node);}

template <typename NodeVal>
void reBrokenStick<NodeVal>::updateRelations() {
	for(int i = 0; i < tree.size(); ++i) updateRelation(i);
}

template <typename NodeVal>
void FenwickNode<NodeVal>::add(const float value){
	_add(value);
	node += value;
	if(node < 0.0) node = 0.0; 
}

template <typename NodeVal>
void FenwickNode<NodeVal>::update(const float value){
	if(value < 0.0) throw std::runtime_error("ERROR: new value can not be negative!\n");
	add(value - node);
}

template <typename NodeVal>
void FenwickNode<NodeVal>::zero(){
	update(0.0);
}

template <typename NodeVal>
void FenwickNode<NodeVal>::print(){
	std::cout 
		<< node << '\t'; 
	std::cout 
		<< branch << '\t'; 
	std::cout 
		<< dummyCumsum() << '\t';
	std::cout 
		<< content << std::endl;
} 

template <typename NodeVal>
NodeVal & FenwickNode<NodeVal>::operator * (){
	return(content);
}

template <typename NodeVal>
void FenwickNode<NodeVal>::_add(const float value){
	if(value == 0.0) return;

	branch += static_cast<double>(value);
	if(branch < 0) branch = 0.0;

	if(parent != nullptr) parent->_add(value);
}

// to find largest node with with prefix_sum(i) <= value.
template <typename NodeVal>
FenwickNode<NodeVal> * FenwickNode<NodeVal>::find(double rn){
	if(branch > rn){ // it should not be me, someone lower  
		if(firstchild != nullptr){
			return(firstchild->find(rn));
		} else { // it is a ME!
			return this;
		}
	} else { // it is me or my parent or my siblings or their children
		if(sibling != nullptr){ // it is my siblings` families (or my parent)
			return(sibling->find(rn-branch)); // but is is definitely not me nor my branch so dont search in here!
		} else if(parent != nullptr){ // it is my parent!
			return parent;
		} else { // it is a ME!
			return this;
		}
	}
}

template <typename NodeVal>
double FenwickNode<NodeVal>::dummyCumsum(){
	if(parent == nullptr) return branch;

	double sum = branch;
	auto ptr = parent->firstchild; 
	while(ptr != this) {
		sum += ptr->branch;
		ptr = ptr->sibling;
	}
	return(sum);
}

}
#endif
