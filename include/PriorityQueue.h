#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <stdexcept>
#include <vector>
#include "heap.h"

/**
 * @param DataType : type of elements conatined in the p-queue
 * @param Compare Class : sort convention for the elements in the p-queu
 * Default Compare = GreaterThan -> Elements with larger priority on top in the p-queue
 */

template <typename DataType , typename Compare, typename Swap >
class PriorityQueue{
	
	public :
	
	typedef std::vector<DataType> Sequence; 
	
	/**
	 * @Constructor
	 * Default : Empty Priority Queue
	 */
	PriorityQueue(Compare* compare, Swap* swap):data_(compare,swap){};
	
	/*
	 * @Constructor
	 * Optional : PriorityQueue built out of a vector<DataType> of elements 
	 */
	PriorityQueue(Sequence& elements, Compare* compare, Swap* swap):data_(elements, compare, swap){
		
	};
	
	/*
	 * @Destructor
	 */
	~PriorityQueue(){};
	
	int  size() const { return data_.size();}
	
	/*
	 * @error : Priority Queue might still be empty when requesting top
	 */
	const DataType& top() const throw (std::runtime_error){ 
		if(empty()){
	 		throw std::runtime_error("PriorityQueue::Pop -> Empty Priority Queue");
	 	}else{
	 		return data_.top();
	 	}
	};
	
	/*
	 * @error : Priority Queue might still be empty when requesting pop
	 */
	const DataType pop() throw (std::runtime_error) {
		if(empty()){
			throw std::runtime_error("PriorityQueue::Pop -> Empty Priority Queue");
		}else{
			return data_.pop();
		}
	};
	
	void setSentinel(const DataType& element) { data_.setSentinel(element);};

	void push(const DataType& element) { data_.push(element);};
	
	void print() const { return data_.print();}
	
	bool empty() const { return data_.empty();}
	
	void clear() {data_.clear();}
	
	void update(int index) {data_.update(index);}
	
	private : 
	
	// ---------------------------
	// DATA MEMBERS
	// ---------------------------
	// PriorityQueue Uses a Heap internally
	
	Heap<DataType , Compare, Swap> data_;
	
};

#endif
