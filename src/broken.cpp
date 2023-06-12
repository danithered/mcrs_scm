#include <list>
#include <vector>
#include <iostream>

#define LSB(x) ((x) & -(x))
//inline int LSB(int x) {return(x & -x);}


class FenwickNode {
	public:

		int * content; // stored value
		// Constructor
		FenwickNode() : branch(0.0), node(0.0), parent(nullptr), firstchild(nullptr), sibling(nullptr) {std::cout << node << "\n";}

		// Functions
		void add(const double value){
			_add(value);
			node += value;
		}

		void update(const double value){
			add(value - node);
		}

		void zero(){
			update(0.0);
		}

		void print(){
			std::cout 
				<< node << '\t'; 
			std::cout 
				<< branch << '\t'; 
			std::cout 
				<< dummyCumsum() << '\t';
			std::cout 
				<< *content << std::endl;
		}

		friend class reBrokenStick;
	private:
		double branch; // the weight of the branch
		double node; // the weight of the node
		FenwickNode * parent;
		FenwickNode * firstchild;
		FenwickNode * sibling;

		void _add(const double value){
			if(value == 0.0) return;

			branch += value;

			if(node < 0.0) throw std::runtime_error("Invald modification of node in reBrokenStick: propensity cant be negative!\n"); 

			if(parent != nullptr) parent->add(value);
		}

		// to find node, where 
		FenwickNode * find(double rn){
			if(branch > rn){ // it should not be me, someone lower  
				return(firstchild->find(rn));
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

		double dummyCumsum(){
			if(parent == nullptr) return branch;

			double sum = branch;
			auto ptr = parent->firstchild; 
			while(ptr != this) {
				sum += ptr->branch;
				ptr = ptr->sibling;
			}
			return(sum);
		}

};

class reBrokenStick {
	public:
		reBrokenStick(unsigned int _size): used(0), tree(_size){
			init(_size);
			updateRelations();
		}

		void init(const unsigned int no_nodes){
			if( no_nodes > tree.size() ){ // if new number of values stored does not exceeds size() no need to increase space
				// calculate neccessary nodes
				unsigned int newsize = tree.size() + 1;
				if (newsize % 2) {
					newsize = (1 << static_cast<unsigned int>(newsize / 2) + 1);
				}

				if(newsize > tree.capacity()) throw std::runtime_error("Undefined\n");

				// add empty nodes
				if(newsize > tree.size()){
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
		}

		void push_back(double p, int * val){ 
			print();
			init(used+1);
			// now node exists for sure, lets assign value to it!
			tree[used++].add(p);
			tree[used].content = val;
		}

		int * draw(double rn){
			// init correct random value
			if(rn<0.0 || rn >1.0) throw std::invalid_argument("Invalid random number provided to reBrokenStick. It should be between 0 and 1!\n");
			rn *= cumsum();

			// check 0 position 
			if(rn < tree[0].branch) return tree[0].content; // Sometimes this is the biggest
			
			auto out = tree.back().find(rn);
			return (++out)->content;
		    
		}

		double cumsum() const {return tree.back().branch;}

		void print(){
			std::cout << "header\n";
			for (auto & node : tree) node.print(); 
		}
	private:
		unsigned int used; // number of nodes used
		std::vector<FenwickNode> tree;

		void updateRelation(unsigned int index){
			if(index == 0){
				tree[index].firstchild = nullptr;
				tree[index].sibling= nullptr;
				tree[index].parent = (tree.size()>1)?&tree[1]:nullptr;
			}

			const unsigned int my_LSB = LSB(index);

			// parent
			int target = index + my_LSB;
			tree[index].parent = (target < tree.size())?&tree[target]:nullptr;


			if(index & 1){ // there is a bit on the first position -> no child, no sibling 
				       
				// firstchild
				if(index == 1) { // it is the first node -> its chid is the dummy 0
					tree[index].firstchild = &tree[0];
				}
				tree[index].firstchild = nullptr;

				// sibling
				tree[index].sibling= nullptr;
			} else { // it has at least one child and one sibling
				 
				// firstchild
				tree[index].firstchild = &tree[(index ^ my_LSB) & (my_LSB >> 1) ];

				// sibling
				target = index + (1 >> my_LSB);
				tree[index].sibling = (target < tree.size())?&tree[target]:nullptr;
			}
			std::cout << "updated" << tree[index].branch << std::endl;
		}

		void updateRelations() {
			for(int i = 0; i < tree.size(); ++i) updateRelation(i);
		}
};


class FenwickTree {
public:
    std::vector<int> A;
    int SIZE;
    
    FenwickTree(int size) : A(size + 1, 0), SIZE(size+1) {}

    
    int range_sum(int i, int j) {
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
    }

// Add delta to element with index i (zero-based)
void add(int i, int delta) {
	if (i == 0) {
		A[0] += delta;
		return;
	}
	for (; i < SIZE; i+= LSB(i))
		A[i] += delta;
}
    
    // Convert A[] in place to Fenwick tree form
    void init(void) {
    	for (int i = 1; i < SIZE; ++i) {
    		int j = i + LSB(i);
    		if (j < SIZE)
    			A[j] += A[i];
    	}
    }
    
    // Convert back to array of per-element counts
    void fini(void) {
    	for (int i = SIZE - 1; i > 0; --i) {
    		int j = i + LSB(i);
    		if (j < SIZE)
    			A[j] -= A[i];
    	}
    }
    
    // Return a single element's value
    int get(int i) {
        if(i == -1) return(A[0]);
    	return range_sum(i, i + 1);
    }
    
    // Set (as opposed to adjust) a single element's value
    void set(int i, int value) {
    	add(i, value - get(i));
    }
    
    // Find the largest i with prefix_sum(i) <= value.
    // NOTE: Requires that all values are non-negative!
    unsigned int rank_query(int value) {
	if(value < A[0]) return 0; // Claim empty - sometimes this is the biggest
    	int i = 0, j = SIZE - 1;
    	// j is a power of 2.
    
        value -= A[0];
    	for (; j > 0;  j >>= 1) {
		std::cout << "I am on " << (i+j) << std::endl;
    	    //std::cout << "for: j=" << j << ", i=" << i << ", value: " << value << ", fval: " << A[i+j] << std::endl;
    		if (i + j < SIZE && A[i + j] <= value) {
    			value -= A[i + j];
    			i += j;
    		    //std::cout << "\tafter: j=" << j << ", i=" << i << ", value: " << value << std::endl;
    		    std::cout << "found smaller value" << std::endl;
    		}
    	}
    	return i+1;
    }
    
    void print(){
        std::cout << "ind";
        for(int i=0; i<SIZE; ++i) std::cout << '\t' << i;
        std::cout << std::endl;
        
        std::cout << "fva";
        for(auto &e : A) std::cout << '\t' << e;
        std::cout << std::endl;
        
        std::cout << "val"; 
        for(int i = -1; i < (SIZE-1); i++) {
            std::cout << "\t" << get(i); 
        }
        std::cout << std::endl;
    
        std::cout << "cum"; 
        for(int i = -1, cs = 0; i < (SIZE-1); i++) {
            cs += get(i);
            std::cout << "\t" << cs; 
        }
        std::cout << std::endl;
    
    }
};

int main(){
    /*FenwickTree tree(16);
    tree.set(1, 3);
    tree.set(3, 2);
    tree.set(5, 1);
    tree.set(6, 1);
    tree.set(10, 4);
    tree.set(12, 1);
    tree.set(15, 1);

    tree.print();
    
    std::cout << tree.rank_query(10) << std::endl;
    std::cout << tree.rank_query(11) << std::endl;
    std::cout << tree.rank_query(3) << std::endl;
    std::cout << tree.rank_query(12) << std::endl;
    */

	std::vector<int> vals{0,1,2,3,4}; 
	reBrokenStick x(16);
//	x.push_back(0,&vals[0]);
//	x.push_back(0,&vals[1]);
//	x.push_back(1,&vals[2]);
//	x.push_back(0,&vals[3]);
	
	x.print();

    return 0;
}
