/*
 * Skeleton_blocker_complex_iterators.h
 *
 *  Created on: Feb 4, 2014
 *      Author: dsalinas
 */

#ifndef SKELETON_BLOCKER_COMPLEX_ITERATORS_H_
#define SKELETON_BLOCKER_COMPLEX_ITERATORS_H_

/**
 *@class Complex_edge_iterator
 *@brief Iterator on the vertices of a simplicial complex
 *
 */
template<typename Traits>
class Skeleton_blocker_complex<Traits>::Complex_vertex_iterator
{
	friend class Skeleton_blocker_complex<Traits> ;
private:
	const Skeleton_blocker_complex<Traits>* complex;
	std::pair<boost_vertex_iterator, boost_vertex_iterator> vertexIterator;


public:
	Complex_vertex_iterator():complex(NULL){
	}

private:
	Complex_vertex_iterator(const Skeleton_blocker_complex<Traits>* complex_):
		complex(complex_),
		vertexIterator(vertices(complex_->skeleton))
{
}

public:
	Complex_vertex_iterator(std::pair<boost_vertex_iterator, boost_vertex_iterator>& pair):
		vertexIterator(pair),complex(NULL){
	}

	bool operator==(const Complex_vertex_iterator& other){
		return (vertexIterator == other.vertexIterator &&  complex == other.complex);
	}

	bool operator!=(const Complex_vertex_iterator& other){
		return(! (*this == other));
	}
private:
	void gotoNextVertex(){
		if(!finished()){
			++vertexIterator.first;
			bool is_active = (complex->skeleton[*vertexIterator.first]).is_active();
			if(!is_active) gotoNextVertex();
		}
	}

public:
	Complex_vertex_iterator& operator++(){
		gotoNextVertex();
		return(*this);
	}

	Complex_vertex_iterator operator++(int){
		Complex_vertex_iterator tmp(*this);
		gotoNextVertex();
		return(tmp);
	}

	Vertex_handle  operator*()	{
		// these two lines are just in case the first vertex isnt active
		bool is_active = (complex->skeleton[*vertexIterator.first]).is_active();
		if(!is_active) gotoNextVertex();
		return(Vertex_handle(*(vertexIterator.first)));
	}

	bool finished() const{
		return vertexIterator.first == vertexIterator.second;
	}
};
//////////////////////////////////////////////////////////////////////



/**
 *@class Complex_edge_iterator
 *@brief Iterator on edges of a simplicial complex
 *
 */
template<typename Traits>
class Skeleton_blocker_complex<Traits>::Complex_edge_iterator
{
	friend class Skeleton_blocker_complex<Traits> ;

	template<typename L>
	friend class Skeleton_blocker_complex<Traits>::Triangle_around_vertex_iterator;

	typedef Skeleton_blocker_complex<Traits> Complex;
	typedef typename Complex::boost_edge_iterator boost_edge_iterator;
	typedef typename Complex::Edge_handle Edge_handle;

private:

	const Complex* complex;
	std::pair<boost_edge_iterator,boost_edge_iterator> edge_iterator ;

public:

	Complex_edge_iterator():complex(NULL){
	}
	Complex_edge_iterator(const Skeleton_blocker_complex<Traits>* complex_):
		complex(complex_),
		edge_iterator(boost::edges(complex_->skeleton))
	{
	}

//	Complex_edge_iterator& operator=(const Complex_edge_iterator& other)
//	{
//		complex = other.complex;
//		edge_iterator = other.edge_iterator;
//		return(*this);
//	}

	bool operator==(const Complex_edge_iterator& other) const{
		return (complex== other.complex) &&(edge_iterator == other.edge_iterator);
	}

	bool operator!=(const Complex_edge_iterator& other)const{
		return(! (*this == other));
	}

	Complex_edge_iterator& operator++(){
		if(edge_iterator.first != edge_iterator.second){
			++(edge_iterator.first);
		}
		return(*this);
	}

	Complex_edge_iterator operator++(int){
		Complex_edge_iterator tmp(*this);
		if(edge_iterator.first != edge_iterator.second){
			++(edge_iterator.first);
		}
		return(tmp);
		// does 'return edge_iterator++' works?
	}

	Edge_handle operator*()	{
		return(*(edge_iterator.first));
	}
};




/**
 * \brief Iterator over the triangles that are
 * adjacent to a vertex of the simplicial complex.
 */
template<typename Traits>
template<typename LinkType>
class Skeleton_blocker_complex<Traits>::Triangle_around_vertex_iterator{

private:
	typedef typename LinkType::Vertex_handle Vertex_handle;
	typedef typename LinkType::Root_vertex_handle Root_vertex_handle;
	typedef typename LinkType::Simplex_handle Simplex_handle;
	typedef Skeleton_blocker_complex Complex;
	typedef typename Skeleton_blocker_complex<Traits>::Complex_edge_iterator Complex_edge_iterator;

	Complex* complex_;
	Vertex_handle v_;
	LinkType link_;
	Complex_edge_iterator current_edge_;
	bool is_end_;
public:
	Triangle_around_vertex_iterator(Skeleton_blocker_complex<Traits>* complex,Vertex_handle v):
		complex_(complex),v_(v),link_(*complex,v_),current_edge_(link_.edge_range().begin()),is_end_(current_edge_ == link_.edge_range().end()){
}
	/**
	 * @brief ugly hack to get an iterator to the end
	 */
	Triangle_around_vertex_iterator(Skeleton_blocker_complex<Traits>* complex,Vertex_handle v,bool is_end):
		complex_(complex),v_(v),link_(),is_end_(true){
	}

	// For now, the operator points to the begin of triangles
	Triangle_around_vertex_iterator& operator=(const Triangle_around_vertex_iterator& other){
		v_ = other.v_;
		complex_ = other.complex_;
		link_ = other.link_;

		current_edge_= link_.edge_range().begin();
		//while (*current_edge_ != *(other.current_edge_))
			//++current_edge_;
		is_end_ = other.is_end_;
		return *this;
	}

	bool operator ==(const Triangle_around_vertex_iterator& other) const{
		return (complex_==other.complex_) && ((finished() &&other.finished()) || current_edge_ == other.current_edge_);
	}

	bool operator !=(const Triangle_around_vertex_iterator& other) const{

		return !(*this == other);
	}

	Simplex_handle operator*(){
		Root_vertex_handle v1 = link_[*current_edge_].first();
		Root_vertex_handle v2 = link_[*current_edge_].second();
		return Simplex_handle(v_,*(complex_->get_address(v1)),*(complex_->get_address(v2)));
	}

	bool finished() const{
		return is_end_;
	}

	Triangle_around_vertex_iterator& operator++(){
		++current_edge_;
		is_end_ = (current_edge_ == link_.edge_range().end());
		return(*this);
	}
};



/**
 * \brief Iterator over the triangles of the
 * simplicial complex.
 *
 */
template<typename Traits>
class Skeleton_blocker_complex<Traits>::Triangle_iterator{

private:
	typedef typename Skeleton_blocker_complex<Traits>::Vertex_handle Vertex_handle;
	typedef typename Skeleton_blocker_complex<Traits>::Root_vertex_handle Root_vertex_handle;
	typedef typename Skeleton_blocker_complex<Traits>::Simplex_handle Simplex_handle;
	typedef typename Skeleton_blocker_complex<Traits>::Complex_vertex_iterator Complex_vertex_iterator;
	typedef Skeleton_blocker_complex<Traits> Complex;
	typedef Complex::Superior_triangle_around_vertex_iterator STAVI;

	Complex* complex_;
	Complex_vertex_iterator current_vertex_;
	STAVI current_triangle_;
	bool is_end_;
public:

	/*
	 * @remark  assume that the complex is non-empty
	 */
	Triangle_iterator(Skeleton_blocker_complex<Traits>* complex):
		complex_(complex),
		current_vertex_(complex->vertex_range().begin()),
		current_triangle_(complex,*current_vertex_), // xxx this line is problematic is the complex is empty
		is_end_(false){

		assert(!complex->empty());

	}

	/**
	 * @brief ugly hack to get an iterator to the end
	 * @remark  assume that the complex is non-empty
	 */
	Triangle_iterator(Skeleton_blocker_complex<Traits>* complex,bool is_end):
		complex_(complex),
		current_vertex_(complex->vertex_range().begin()),
		current_triangle_(complex->superior_triangle_range(*current_vertex_).end()), // xxx this line is problematic is the complex is empty
		is_end_(true){
		assert(!complex->empty());
	}


	Triangle_iterator& operator=(const Triangle_iterator & other){
		complex_ = other.complex_;
		Complex_vertex_iterator current_vertex_;
		STAVI current_triangle_;
		bool is_end_;
		return *this;
	}


	bool operator ==(const Triangle_iterator& other) {

		bool both_are_finished = is_finished() && other.is_finished();
		bool both_arent_finished = !is_finished() && !other.is_finished();
		// if the two iterators are not finished, they must have the same state
		return (complex_==other.complex_) && (both_are_finished	||
				( (both_arent_finished) && current_vertex_ == other.current_vertex_ && current_triangle_ == other.current_triangle_));

	}

	bool operator !=(const Triangle_iterator& other){
		return !(*this == other);
	}

	Simplex_handle operator*(){
		return *current_triangle_;
	}

	void goto_next_vertex(){
		++current_vertex_;
		is_end_ = (current_vertex_ == complex_->vertex_range().end());
		if (!is_end_){
			STAVI other(complex_, *current_vertex_);
			current_triangle_ = other;
		}
	}

	Triangle_iterator& operator++(){
		if(!current_triangle_.finished()){
			++current_triangle_;
			while(current_triangle_.finished()&& !is_finished()){
				goto_next_vertex();
			}
		}
		else{
			goto_next_vertex();
		}
		return(*this);
	}

	bool is_finished() const{
		return is_end_;
	}


};





/**
 * @brief Iterator through the blockers of a vertex
 */
template<typename Traits>
class Skeleton_blocker_complex<Traits>::Complex_const_blocker_iterator{
	friend class Skeleton_blocker_complex<Traits> ;
	typedef Skeleton_blocker_complex<Traits> Complex;
	typedef typename Complex::BlockerMapConstIterator BlockerMapConstIterator;

private:
	const Complex * complex;
	BlockerMapConstIterator currentPosition;
public:

	Complex_const_blocker_iterator():
		complex(0),
		currentPosition()

{}

	Complex_const_blocker_iterator(const Skeleton_blocker_complex * complex_,BlockerMapConstIterator position):
		complex(complex_),
		currentPosition(position)
	{}

	bool operator==(const Complex_const_blocker_iterator& other){
		return
				currentPosition == other.currentPosition
				&& complex == (other.complex);
	}

	bool operator!=(const Complex_const_blocker_iterator& other){
		return(! (*this == other));
	}

	Complex_const_blocker_iterator& operator++(){
		currentPosition++;
		return(*this);
	}

	const Simplex_handle* operator*()	{
		return(currentPosition->second);
	}
};



/**
 * @brief Iterator through the blockers of a vertex
 */
// ReturnType = const Simplex_handle* or Simplex_handle*
// MapIteratorType = BlockerMapConstIterator or BlockerMapIterator
template<typename Traits>
template<typename MapIteratorType, typename ReturnType>
class Skeleton_blocker_complex<Traits>::Complex_blocker_iterator{
	friend class Skeleton_blocker_complex<Traits> ;
private:
	MapIteratorType currentPosition;
public:

	Complex_blocker_iterator():currentPosition(){}

	Complex_blocker_iterator(MapIteratorType position):
		currentPosition(position)
	{}

	bool operator==(const Complex_blocker_iterator& other) const{
		return currentPosition == other.currentPosition;
	}

	bool operator!=(const Complex_blocker_iterator& other){
		return(! (*this == other));
	}

	Complex_blocker_iterator& operator++(){
		currentPosition++;
		return(*this);
	}

	ReturnType operator*()	{
		return(currentPosition->second);
	}
};






template<typename Traits>
class Skeleton_blocker_complex<Traits>::Complex_vertex_range {
private:
	const Skeleton_blocker_complex<Traits>* complex_;
public:
	Complex_vertex_range(const Skeleton_blocker_complex<Traits>* complex):complex_(complex){
	}

	Complex_vertex_iterator begin(){
		return Complex_vertex_iterator(complex_);
	}

	Complex_vertex_iterator end(){
		Complex_vertex_iterator vIt(complex_);
		vIt.vertexIterator.first = vIt.vertexIterator.second ;
		return vIt;
	}
};

template<typename Traits>
class Skeleton_blocker_complex<Traits>::Complex_edge_range {
private:
	const Skeleton_blocker_complex* complex_;
public:
	Complex_edge_range(const Skeleton_blocker_complex* complex):complex_(complex){
	}

	Complex_edge_iterator begin(){
		return (Complex_edge_iterator(complex_));
	}

	Complex_edge_iterator end()
	{
		Complex_edge_iterator edge_it(complex_);
		edge_it.edge_iterator.first = edge_it.edge_iterator.second;
		return edge_it;
	}
};


template<typename Traits>
template<typename LinkType>
class Skeleton_blocker_complex<Traits>::Triangle_around_vertex_range {
private:
	typedef Skeleton_blocker_complex<Traits> Complex;
	typedef Complex::Triangle_around_vertex_iterator<LinkType> Tavi;
	typedef typename Skeleton_blocker_complex<Traits>::Vertex_handle Vertex_handle;
	Complex* complex_;
	Vertex_handle v_;
public:
	Triangle_around_vertex_range(Skeleton_blocker_complex<Traits>* complex,Vertex_handle v):
		complex_(complex),v_(v){
	}

	Tavi begin(){
		return Tavi(complex_,v_);
	}

	Tavi end()
	{
		return Tavi(complex_,v_,true);

	}
};


template<typename Traits>
class Skeleton_blocker_complex<Traits>::Triangle_range {
private:
	typedef Skeleton_blocker_complex<Traits> Complex;
	typedef Complex::Triangle_iterator Triangle_iterator;
	typedef typename Skeleton_blocker_complex<Traits>::Vertex_handle Vertex_handle;
	Complex* complex_;
public:
	Triangle_range(Skeleton_blocker_complex<Traits>* complex):
		complex_(complex){
	}

	Triangle_iterator begin(){
		return Triangle_iterator(complex_);
	}

	Triangle_iterator end()
	{
		return Triangle_iterator(complex_,true);
	}
};





template<typename Traits>
class Skeleton_blocker_complex<Traits>::Complex_const_blocker_range {
private:
	const Skeleton_blocker_complex* complex_;
	Vertex_handle v_;
public:
	Complex_const_blocker_range(const Skeleton_blocker_complex* complex,Vertex_handle v):complex_(complex),v_(v){
	}

	Complex_const_blocker_iterator begin(){
		return Complex_const_blocker_iterator(complex_,complex_->blocker_map.lower_bound(v_));

	}

	Complex_const_blocker_iterator end()
	{
		return Complex_const_blocker_iterator(complex_,complex_->blocker_map.upper_bound(v_));
	}
};









#endif /* SKELETON_BLOCKER_COMPLEX_ITERATORS_H_ */
