#ifndef HEAP_H
#define HEAP_H

#include <cmath>
#include <vector>
#include <iostream>
#include <exception>

using namespace std;
/**
 * @Mainpage
 * HEAP Template Class Implementation
 * Note that position 0 of the heap is not used, a sentinel is placed there
 */

class noSentinelException: public exception {
	virtual const char* what() const throw(){
		return "Heap was not initialized with a sentinel";
	}
};

class emptyHeapException: public exception {
	virtual const char* what() const throw(){
		return "Heap does not contain any elements";
	}
};

template <typename DataType, typename Compare, typename Swap> 
class Heap {
	
public:
	
	//------------------------
	// PUBLIC INTERFACE
	//------------------------
	
	typedef std::vector<DataType> Sequence;
	
   	/**
	 * @Constructor
	 * Default Constructor
	 * Empty Heap
	 */
	Heap(Compare* compare, Swap* swap):compare_(compare), swap_(swap), sentinel_(false){};
	
	
	~Heap(){};
	
	int size() const {
		return (int) contents_.size();
	};
	
	void clear()  {
		contents_.clear();
	};
	
	/**
	 * @return DataType element w. highest priority
	 */
	const DataType& top() const {
		if(sentinel_)
			return contents_[1];
		throw noSentinelException();
	};
	
	void print() const  {
		if(contents_.size() > 1){
			for (unsigned i = 1; i < contents_.size(); i++)
			{
				contents_[i]->print();
			}
		}
		else throw emptyHeapException();
	}
	
	/**
	 * Store the Top Element. 
	 * Put the last (temporaraly) on front as root.
	 * Shrink the Container by Popping it once.
	 * Restore the Heap property by applying heapifyDown() from the root.  
	 * @return DataType element w. highest priority
	 */
	const DataType pop() {
		if(sentinel_){
			if(contents_.size() <= 1) throw emptyHeapException();
			DataType max = contents_[1];
			contents_[1] = contents_.back();
			contents_.pop_back();
			heapifyDown(1);
			return max;
		}
		throw noSentinelException();
	};
	
	/**
	 * Before the heap can be used the sentinel needs to be added. 
	 * @param element of DataType
	 * Store the new element at the bottom of the heap
	 * Restore the Heap property by applying heapifyUp from the new element's index up
	 * NOTE : HeapifyUp implemented as a while loop, returns location of element in heap 
	 * Worst Case : Pushed element becomes the new root => Performance O (log n )
	 */
	void setSentinel(const DataType& element){
		if (contents_.size() == 0)		
			contents_.push_back(element);	
		else contents_[0] = element;
		sentinel_ = true;
	};
	
	/**
	 * Store the Top Element. 
	 * @param element of DataType
	 * Store the new element at the bottom of the heap
	 * Restore the Heap property by applying heapifyUp from the new element's index up
	 * NOTE : HeapifyUp implemented as a while loop, returns location of element in heap 
	 * Worst Case : Pushed element becomes the new root => Performance O (log n )
	 */
	void push(const DataType& element){
		if(sentinel_){
			int i = (int) contents_.size();
			contents_.push_back(element);	
//			cout << i << "\t" << *contents_[i] << endl;

			while((i > 1) &! compare_->smaller(contents_[parent(i)],contents_[i])){
//				cout <<" -->parent " <<  parent(i) << "\t" << *contents_[parent(i)] << endl;
				swap_->exchange(contents_[parent(i)] , contents_[i]);
//				cout << "Shifted up " << endl;
				i = parent(i);
//				cout <<" -->new " <<  i << "\t" << *contents_[i] << endl;
			}
		}
		else throw noSentinelException();
	};
	
	bool empty() const {
		return (contents_.size() <= 1);
	}
	
	
	/**
	 * Compares elements, using the Compare Template param
	 * Default comparison is GreaterThan
	 * Worst Case : element being heapified, descends to the bottom of the Heap
	 * Performance => O ( log n )
	 */
	void heapifyDown(int index){
		if(sentinel_){
//			cout << index << "\t" << *contents_[index] << endl;
			int l = left(index);
			int r = right(index);
			int selected;
			
			// compare left child with the current node-index
			if(l < size() && compare_->smaller(contents_[l] , contents_[index])){
//				cout <<" -->left " <<  l << "\t" << *contents_[l] << endl;
				selected = l;
			}else{
				selected = index;
			}
			// compare right child with the retained node-index of above comparison
			if(r < size() && compare_->smaller(contents_[r] , contents_[selected])){
//				cout <<" -->right " <<  r << "\t" << *contents_[r] << endl;
				selected = r;	
			}
			// if there is an unabalance [violated heap property] we need to swap 
			// elements at selected indexes
			if(index != selected){
				swap_->exchange(contents_[index] , contents_[selected]);
//				cout << "Shifted down " << endl;
				heapifyDown(selected);
			}
		}
		else throw noSentinelException();
	};
	
	/**
	 * Compares elements, using the Compare Template param
	 * Default comparison is GreaterThan
	 * Worst Case : element being heapified, descends from the top to the bottom of the Heap
	 * Performance => O ( log n )
	 */
	void update(int index){
		if(sentinel_){
//			cout << index << "\t" << *contents_[index] << endl;
			int p = parent(index);
//			cout <<" -->parent " <<  p << "\t" << *contents_[p] << endl;
			if(compare_->smaller(contents_[index] , contents_[p])){
				swap_->exchange(contents_[index] , contents_[p]);
//				cout << "Shifted up " << endl;
//				cout <<" -->new " <<  p << "\t" << *contents_[p] << endl;
				update(p);
			}
			else heapifyDown(index);
		}
		else throw noSentinelException();
	}
	
	
	//------------------------
	// PRIVATE METHODS
	//------------------------
	
private:
	/** finds parent of node index.
	 * New Depths are being filled up sequentially, a Heap is a complete Binary Tree
	 * Each new depth has Max.number of Nodes equal to Sum of all Nodes higher in the
         * Heap minus 1. Therefore parentnode for index or index+1 is always the rounded 
         * down value of the index divided by 2.
	 */
	int parent(int index){
		int tmp = (int)floor(index / 2);
		return (tmp >= 1?tmp:1); 
	};
	
	int left(int index){
		return 2 * index;	
	};
	
	int right(int index){
		return (2 * index)+1;	
	};
	
	//-------------
	// DATA MEMBERS
	//-------------
	
	Sequence 	contents_;
	Compare*	compare_;
	Swap*		swap_;
	bool		sentinel_;
	
};
	
#endif
