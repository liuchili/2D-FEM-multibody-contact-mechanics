// -*- C++ -*-
#ifndef ESSENTIALBOUNDARYCONDITIONS_H
#define ESSENTIALBOUNDARYCONDITIONS_H

#include "Utilities_Liuchi.h"

class EssentialBC {
public:
	EssentialBC() {}
	EssentialBC(const size_t & grainId, const size_t & nodeIdGrainBased, 
		const size_t & direction, const double & value) {
		_grainId = grainId;
		_nodeIdGrainBased = nodeIdGrainBased;
		_direction = direction;
		_value = value;
	}

	size_t
	getGrainId() const {
		return _grainId;
	}

	size_t
	getNodeIdGrainBased() const {
		return _nodeIdGrainBased;
	}

	size_t
	getDirection() const {
		return _direction;
	}

	double
	getValue() const {
		return _value;
	}

	void
	changeBoundaryValue (double & value) {
		_value = value;
	}
private:
	size_t _grainId;
	size_t _nodeIdGrainBased;
	size_t _direction;
	double _value;
};

#endif  // ESSENTIALBOUNDARYCONDITIONS_H
