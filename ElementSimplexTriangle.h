// -*- C++ -*-
#ifndef ELEMENT_TRIANGLE
#define ELEMENT_TRIANGLE

#include "Definitions.h"
#include "Utilities.h"
#include "MaterialModelLinearElastic.h"




class SimplexTriangle {
public:
  SimplexTriangle () {}
  SimplexTriangle(const vector<Vector2d> & nodePositions, const vector<size_t> & nodeIds,
                  MaterialModelLinearElastic * materialmodel,
                  const size_t & qpnumber) :
                  _nodeIds(nodeIds), _qpnumber(qpnumber), _materialmodel(materialmodel) {

    _quadpoints.resize(qpnumber);
    _quadweights.resize(qpnumber);
    if (qpnumber == 1) {
      _quadweights[0] = 0.5;
      _quadpoints[0] = Vector2d(1./3.,1./3.);
    }
    if (qpnumber == 3) {
      _quadweights[0] = 1./6.;
      _quadweights[1] = 1./6.;
      _quadweights[2] = 1./6.;

      _quadpoints[0] = Vector2d(0.5,0.);
      _quadpoints[1] = Vector2d(0.,0.5);
      _quadpoints[2] = Vector2d(0.5,0.5);
    }

    Matrix<double, 3, 4> derivativeToStrains;
    derivativeToStrains.fill(0.);
    derivativeToStrains(0,0) = 1;
    derivativeToStrains(1,3) = 1;
    derivativeToStrains(2,1) = 1;
    derivativeToStrains(2,2) = 1;

    _bMatrices.resize(qpnumber);
    _weightedJacobians.resize(qpnumber);
    _nodalWeights.resize(nodeIds.size());

    Matrix<double, 4, 6> shapeFunctionDerivatives; shapeFunctionDerivatives.fill(0.);
    Matrix<double, 2, 2> jacobian; jacobian.fill(0.);
    Matrix<double, 4, 4> blockGamma; blockGamma.fill(0.);
    double volume = 0.;

    for (unsigned int qpindex = 0; qpindex < qpnumber; qpindex++) {
      shapeFunctionDerivatives.fill(0.);
      shapeFunctionDerivatives(0,0) = 1; shapeFunctionDerivatives(0,4) = -1;
      shapeFunctionDerivatives(1,2) = 1; shapeFunctionDerivatives(1,4) = -1;
      shapeFunctionDerivatives(2,1) = 1; shapeFunctionDerivatives(2,5) = -1;
      shapeFunctionDerivatives(3,3) = 1; shapeFunctionDerivatives(3,5) = -1;

      jacobian(0,0) = nodePositions[0](0) - nodePositions[2](0);
      jacobian(0,1) = nodePositions[0](1) - nodePositions[2](1);
      jacobian(1,0) = nodePositions[1](0) - nodePositions[2](0);
      jacobian(1,1) = nodePositions[1](1) - nodePositions[2](1);


      // blockGamma, accordingly, is also constant
      blockGamma.fill(0.);
      blockGamma.block<2,2>(0,0) = jacobian.inverse();
      blockGamma.block<2,2>(2,2) = jacobian.inverse();

      _bMatrices[qpindex] = derivativeToStrains*blockGamma*shapeFunctionDerivatives;
      _weightedJacobians[qpindex] = _quadweights[qpindex]*jacobian.determinant()*1.; //1. is the thickness
      volume += _weightedJacobians[qpindex];
    }

    for (unsigned int index = 0; index < _nodeIds.size(); index++) {
      _nodalWeights[index] = volume/double(_nodeIds.size());
    }

  }

  double
  computeEnergy(const vector<Vector2d> & nodalDisplacements) const {
    double energy = 0.0;
    vector<Vector3d> strainsAtGuassianPoints = computeStrainsAtGuassianPoints(nodalDisplacements);
    for (unsigned int qpindex = 0; qpindex < _qpnumber; qpindex++) {
      const Vector3d strain = strainsAtGuassianPoints[qpindex];
      energy += _weightedJacobians[qpindex] * _materialmodel->computeEnergy(strain);
    }
    return energy;
  }

  vector<Vector2d>
  computeForces(const vector<Vector2d> & nodalDisplacements) const {
    Matrix<double, 6, 1> nodalForcesVector;
    nodalForcesVector.fill(0.);
    Vector3d stress; stress.fill(0.);
    vector<Vector3d> strainsAtGuassianPoints = computeStrainsAtGuassianPoints(nodalDisplacements);
    for (unsigned int qpindex = 0; qpindex < _qpnumber; qpindex++) {
      stress = _materialmodel->computeStressPlaneStrain(strainsAtGuassianPoints[qpindex]);
      nodalForcesVector += _weightedJacobians[qpindex]*_bMatrices[qpindex].transpose()*stress;
    }
    vector<Vector2d> nodalForces(3);
    nodalForces[0] = nodalForcesVector.block<2,1>(0,0);
    nodalForces[1] = nodalForcesVector.block<2,1>(2,0);
    nodalForces[2] = nodalForcesVector.block<2,1>(4,0);
    return nodalForces;
  }

  vector<Vector2d>
  computeBodyForces(const Vector2d & acceleration) const {
    vector<Vector2d> nodeForces(3);
    nodeForces[0].fill(0.);
    nodeForces[1].fill(0.);
    nodeForces[2].fill(0.);
    for (unsigned int qpindex = 0; qpindex < _qpnumber; qpindex++) {
      nodeForces[0](0) += _materialmodel->getDensity()*acceleration(0)*_weightedJacobians[qpindex]*_quadpoints[qpindex](0);
      nodeForces[0](1) += _materialmodel->getDensity()*acceleration(1)*_weightedJacobians[qpindex]*_quadpoints[qpindex](0);

      nodeForces[1](0) += _materialmodel->getDensity()*acceleration(0)*_weightedJacobians[qpindex]*_quadpoints[qpindex](1);
      nodeForces[1](1) += _materialmodel->getDensity()*acceleration(1)*_weightedJacobians[qpindex]*_quadpoints[qpindex](1);

      nodeForces[2](0) += _materialmodel->getDensity()*acceleration(0)*_weightedJacobians[qpindex]*(1-_quadpoints[qpindex](0)-_quadpoints[qpindex](1));
      nodeForces[2](1) += _materialmodel->getDensity()*acceleration(1)*_weightedJacobians[qpindex]*(1-_quadpoints[qpindex](0)-_quadpoints[qpindex](1));
    }
    return nodeForces;
  }

  Matrix<double, 6, 6>
  computeStiffnessMatrix() const {
    Matrix<double, 6, 6> stiffnessMatrix; stiffnessMatrix.fill(0);
    Matrix<double, 3, 3> TangentMatrix = _materialmodel->getTangentMatrixPlaneStrain();
    for (unsigned int qpindex = 0; qpindex < _qpnumber; qpindex++) {
      stiffnessMatrix += _weightedJacobians[qpindex]*_bMatrices[qpindex].transpose()*TangentMatrix*_bMatrices[qpindex];
    }
    return stiffnessMatrix;
  }

  vector<Vector3d>
  computeStrainsAtGuassianPoints (const vector<Vector2d> & nodalDisplacements) const {
    vector<Vector3d> strainsAtGuassianPoints(_qpnumber);
    Matrix<double, 6, 1> displacementVector; displacementVector.fill(0.); 
    displacementVector.block<2,1>(0,0) = nodalDisplacements[0];
    displacementVector.block<2,1>(2,0) = nodalDisplacements[1];
    displacementVector.block<2,1>(4,0) = nodalDisplacements[2];

    for (unsigned int qpindex = 0; qpindex < _qpnumber; qpindex++) {
      strainsAtGuassianPoints[qpindex] = _bMatrices[qpindex]*displacementVector;
    }
    return strainsAtGuassianPoints;
  }

  vector<Vector3d>
  computeStressesAtGaussianPoints(const vector<Vector2d> & nodalDisplacements) const {
    vector<Vector3d> stresses(_qpnumber);
    vector<Vector3d> strainsAtGuassianPoints = computeStrainsAtGuassianPoints(nodalDisplacements);
    for (unsigned int qpindex = 0; qpindex < _qpnumber; qpindex++) {
      stresses[qpindex] = _materialmodel->computeStressPlaneStrain(strainsAtGuassianPoints[qpindex]);
    }
    return stresses;
  }

  vector<size_t>
  getNodeIds() const {
    return _nodeIds;
  }

  vector<double>
  getNodalWeights() const {
    return _nodalWeights;
  }

  void
  updateNodalPositions(const vector<Vector2d> & Positions) {
   
    Matrix<double, 2, 2> jacobian; jacobian.fill(0.);
    Matrix<double, 4, 4> blockGamma; blockGamma.fill(0.);
    double volume = 0.;
    Matrix<double, 4, 6> shapeFunctionDerivatives; shapeFunctionDerivatives.fill(0.);

    Matrix<double, 3, 4> derivativeToStrains;
    derivativeToStrains.fill(0.);
    derivativeToStrains(0,0) = 1;
    derivativeToStrains(1,3) = 1;
    derivativeToStrains(2,1) = 1;
    derivativeToStrains(2,2) = 1;

    for (unsigned int qpindex = 0; qpindex < _qpnumber; qpindex++) {
      shapeFunctionDerivatives.fill(0.);
      shapeFunctionDerivatives(0,0) = 1; shapeFunctionDerivatives(0,4) = -1;
      shapeFunctionDerivatives(1,2) = 1; shapeFunctionDerivatives(1,4) = -1;
      shapeFunctionDerivatives(2,1) = 1; shapeFunctionDerivatives(2,5) = -1;
      shapeFunctionDerivatives(3,3) = 1; shapeFunctionDerivatives(3,5) = -1;

      jacobian(0,0) = Positions[0](0) - Positions[2](0);
      jacobian(0,1) = Positions[0](1) - Positions[2](1);
      jacobian(1,0) = Positions[1](0) - Positions[2](0);
      jacobian(1,1) = Positions[1](1) - Positions[2](1);


      // blockGamma, accordingly, is also constant
      blockGamma.fill(0.);
      blockGamma.block<2,2>(0,0) = jacobian.inverse();
      blockGamma.block<2,2>(2,2) = jacobian.inverse();

      _bMatrices[qpindex] = derivativeToStrains*blockGamma*shapeFunctionDerivatives;
      _weightedJacobians[qpindex] = _quadweights[qpindex]*jacobian.determinant()*1.;
      volume += _weightedJacobians[qpindex];
    }

    for (unsigned int index = 0; index < _nodeIds.size(); index++) {
      _nodalWeights[index] = volume/double(_nodeIds.size());
    }
    
  }


private:
  vector<size_t>                               _nodeIds;
  size_t                                       _qpnumber;
  MaterialModelLinearElastic                   * _materialmodel;
  vector<double>                               _weightedJacobians;
  vector<double>                               _nodalWeights;
  vector<Matrix<double, 3, 6> >                _bMatrices;
  vector<Vector2d>                             _quadpoints;
  vector<double>                               _quadweights;

};


#endif //ELEMENT_TRIANGLE

