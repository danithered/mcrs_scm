#include <broken.hpp>


template <typename NodeVal>
void FenwickNode<NodeVal>::add(const double value){
	_add(value);
	node += value;
	if(node < 0.0) throw std::runtime_error("Invald modification of node in reBrokenStick: propensity cant be negative!\n"); 
}

template <typename NodeVal>
void FenwickNode<NodeVal>::update(const double value){
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
void FenwickNode<NodeVal>::_add(const double value){
	if(value == 0.0) return;

	branch += value;

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


template <typename NodeVal>
reBrokenStick<NodeVal>::reBrokenStick(unsigned int _size): used(0), tree(0){
	init(_size);
	//updateRelations();
}

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
void reBrokenStick<NodeVal>::push_back(double p, int val){ 
	//print();
	init(used);
	// now node exists for sure, lets assign value to it!
	tree[used].add(p);
	tree[used++].content = val;
}

// Find the largest i with prefix_sum(i) <= value.
template <typename NodeVal>
int reBrokenStick<NodeVal>::draw(double rn){
	// init correct random value
	if(rn<0.0 || rn >1.0) throw std::invalid_argument("Invalid random number provided to reBrokenStick. It should be between 0 and 1!\n");
	rn *= cumsum();

	// check 0 position 
	if(rn < tree[0].branch) return tree[0].content; // Sometimes this is the biggest
	
	auto out = tree.back().find(rn);
	return (out)->content;
    
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
	int target = index + my_LSB;
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
double FenwickNode<NodeVal>::get_branch() const {return branch;}

template <typename NodeVal>
double FenwickNode<NodeVal>::get_p() const {return node;}

template <typename NodeVal>
void reBrokenStick<NodeVal>::updateRelations() {
	for(int i = 0; i < tree.size(); ++i) updateRelation(i);
}

/*int main(){

	std::vector<std::pair<double, int>> vals{{000, 100},{0,101},{1,102},{0,103},{0,104}, {2,105}, {0,106}, {1,107}, {1,108}, {000, 109},{0,110},{0,111},{0,112},{0,113}, {0,114}, {0,115}, {0,116}}; 
	reBrokenStick<int> x(16);
	for(auto & [p, val] : vals ){x.push_back(p, val);}
	
	x.print();
	std::cout << "draws:\t" << x.draw(0.1) << '\t' << x.draw(0.4) << '\t' << x.draw(0.9) << std::endl;

	auto& p = x[2];
	
	//std::cout << **p << '\t' << p->get_branch() << '\t' << p->get_p() << std::endl;
	std::cout << *p << '\t' << p.get_branch() << '\t' << p.get_p() << std::endl;
	
	//p->update(5);
	p.update(5);
	
	x.print();
	

    return 0;
}
*/

