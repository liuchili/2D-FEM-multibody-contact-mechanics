// -*- C++ -*-
#ifndef RigidStraightWall_H
#define RigidStraightWall_H

#include "Definitions_Liuchi.h"

class RigidStraightWall {
public:
	RigidStraightWall() {}
	RigidStraightWall(const vector<Vector2d> & boundaryPoints, 
		const Vector2d & normal, const size_t wallId, const double & mu) {
		_boundaryPoints = boundaryPoints;
		_normal = normal;
		_wallId = wallId;
		_mu = mu;
	}

	void
	updateWallPositions(const vector<Vector2d> & displacements) {
		_boundaryPoints[0]+= displacements[0];
		_boundaryPoints[1]+= displacements[1];
		Vector2d tangent = _boundaryPoints[1]-_boundaryPoints[0];
		_normal << -tangent(1), tangent(0);
		_normal /= _normal.norm();
	}

	void
	changeMu(const double & mu) {
		_mu = mu;
	}

	Vector2d
	getWallNormal() const {
		return _normal;
	}

	vector<Vector2d>
	getBoundaryPoints() const {
		return _boundaryPoints;
	}

	size_t
	getWallId() const {
		return _wallId;
	}

	double
	getMu() const {
		return _mu;
	}
private:
	vector<Vector2d>    _boundaryPoints;
	Vector2d            _normal;
	size_t              _wallId;
	double              _mu;
};


#endif /* RIGIDSTRAIGHTWALL_H_ */
