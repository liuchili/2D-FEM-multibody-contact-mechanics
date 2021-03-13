#ifndef MATERIAL_MODEL_LINEAR_ELASTIC
#define MATERIAL_MODEL_LINEAR_ELASTIC

#include "Definitions.h"
#include "Utilities.h"


class MaterialModelLinearElastic {

public:

  MaterialModelLinearElastic(const double & youngsModulus, const double & poissonsRatio, const double & density) {
    _tangentMatrixPlaneStrain.fill(0.);
    _tangentMatrixPlaneStrain(0,0) = 1-poissonsRatio;
    _tangentMatrixPlaneStrain(0,1) = poissonsRatio;
    _tangentMatrixPlaneStrain(1,0) = poissonsRatio;
    _tangentMatrixPlaneStrain(1,1) = 1-poissonsRatio;
    _tangentMatrixPlaneStrain(2,2) = (1-2*poissonsRatio)/2.;
    _tangentMatrixPlaneStrain *= youngsModulus/(1+poissonsRatio)/(1-2*poissonsRatio);

    _tangentMatrixPlaneStress.fill(0.);
    _tangentMatrixPlaneStress(0,0) = 1.;
    _tangentMatrixPlaneStress(0,1) = poissonsRatio;
    _tangentMatrixPlaneStress(1,0) = poissonsRatio;
    _tangentMatrixPlaneStress(1,1) = 1.;
    _tangentMatrixPlaneStress(2,2) = (1-poissonsRatio)/2;
    _tangentMatrixPlaneStress *= youngsModulus/(1-poissonsRatio*poissonsRatio);

    _density = density;

  }
  
  Vector3d
  computeStressPlaneStrain (const Vector3d & strain) {
    return _tangentMatrixPlaneStrain*strain;
    //return _tangentMatrixPlaneStress*strain;
  }

  double
  computeEnergy(const Vector3d & strain) {
    return 0.5*strain.dot(computeStressPlaneStrain(strain));
  }

  Matrix<double, 3, 3>
  const getTangentMatrixPlaneStrain() const {
    return _tangentMatrixPlaneStrain;
    //return _tangentMatrixPlaneStress;
  }

  double
  getDensity() const {
    return _density;
  }


private:
  Matrix<double, 3, 3>   _tangentMatrixPlaneStrain;
  Matrix<double, 3, 3>   _tangentMatrixPlaneStress;
  double                 _density;
};

#endif 