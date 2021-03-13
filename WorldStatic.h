// -*- C++ -*-
#ifndef WORLDSTATIC_H
#define WORLDSTATIC_H

#include "Grain.h"

class WorldStatic {
public:
  WorldStatic () {}
  WorldStatic(const vector<StaticGrain> & grains, const vector<RigidStraightWall> & walls) {
  	_grains = grains;
  	_ngrains = grains.size();
  	_totalNodeNumber = 0;
  	_walls = walls;
  	_blockSizes.resize(_ngrains);
  	for (unsigned int index = 0; index < grains.size(); index++) {
  		_totalNodeNumber += grains[index].getNumberOfNodes();
  		_blockSizes[index] = grains[index].getNumberOfNodes();
  	}
  }

  VectorXd
  assembleWorldBodyForceVector(const Vector2d & acceleration) {
    VectorXd worldBodyForceVector;
    worldBodyForceVector.resize(_totalNodeNumber*2);
    worldBodyForceVector.fill(0.);
    size_t startingPoint = 0;
    for (unsigned int gindex = 0; gindex < _ngrains; gindex++) {
      startingPoint = accumulate(_blockSizes.begin(),_blockSizes.begin()+gindex,0)*2;
      worldBodyForceVector.block(startingPoint, 0, _blockSizes[gindex]*2,
        1) = _grains[gindex].assembleBodyForceVector(acceleration);
    }
    return worldBodyForceVector;
  }

  
  vector<Triplet<double,size_t> >
  assembleWorldStiffnessMatrix(const vector<EssentialBC> & essentialBoundaryConditions) { 
    SparseMatrix<double> worldStiffnessMatrix(_totalNodeNumber*2, _totalNodeNumber*2);
    vector<Triplet<double,size_t> > worldTripletList;
    worldTripletList.reserve(_totalNodeNumber);
    
    vector<Triplet<double,size_t> > singleGrainTripletList;
    size_t startingPoint = 0;
    size_t grainId = 0;
    for (size_t gindex = 0; gindex < _ngrains; gindex++) {
      startingPoint = accumulate(_blockSizes.begin(),_blockSizes.begin()+gindex,0)*2;
      singleGrainTripletList = _grains[gindex].assembleStiffnessMatrix(essentialBoundaryConditions);
      for (size_t tindex = 0; tindex < singleGrainTripletList.size(); tindex++) {
        worldTripletList.push_back(Triplet<double,size_t>(startingPoint+singleGrainTripletList[tindex].row(),
          startingPoint+singleGrainTripletList[tindex].col(),singleGrainTripletList[tindex].value()));
      }
    }
    return worldTripletList;
  }
  
  VectorXd
  assembleWorldContactForceVectorAgainstRigidStraightBoundaries(const vector< vector<Vector2d> > & worldDisplacements, 
    const vector< vector<Vector2d> > & worldDisplacementsAccumulate,
    const vector< vector<Vector2d> > & worldWallDisplacements, const vector< vector<Vector2d> > & worldWallDisplacementsAccumulate,
    const size_t & updateFlag, vector<Triplet<double,size_t> > & worldTripletList, const double & penaltyratio) {
    VectorXd worldContactForceVector;
    worldContactForceVector.resize(_totalNodeNumber*2);
    worldContactForceVector.fill(0.);

    worldTripletList.clear();
    worldTripletList.reserve(_totalNodeNumber);    

    size_t startingPoint = 0;
    vector<Vector2d> singleGrainDisplacements;
    vector<Vector2d> singleWallDisplacements;
    VectorXd singleContactForce;
    VectorXd singleContactForcePerturbed;
    double h;
    size_t boundaryNodeId;
    size_t updatecheck = 0;
    vector<size_t> boundaryNodeIdCheck;
    vector<size_t> boundaryNodeIdCheckTemp;
    for (unsigned int gindex = 0; gindex < _ngrains; gindex++) {
      singleGrainDisplacements = worldDisplacements[gindex];
      startingPoint = accumulate(_blockSizes.begin(),_blockSizes.begin()+gindex,0)*2;
      for (unsigned int windex = 0; windex < _walls.size(); windex++) {
        singleWallDisplacements = worldWallDisplacements[windex];
        if (_grains[gindex].bCircleCheckWithRigidStraightWall(_walls[windex])) {
          singleContactForce = _grains[gindex].assembleContactForceVectorRigidStraightBoundary(singleGrainDisplacements, worldDisplacementsAccumulate[gindex],
            updateFlag, _walls[windex],singleWallDisplacements, worldWallDisplacementsAccumulate[windex], boundaryNodeIdCheck,penaltyratio);
          worldContactForceVector.block(startingPoint,0,_blockSizes[gindex]*2,1) += singleContactForce;
          for (unsigned int bindex = 0; bindex < boundaryNodeIdCheck.size(); bindex++) {
              boundaryNodeId = boundaryNodeIdCheck[bindex];
              for (unsigned int dofIndex = 0; dofIndex < 2; dofIndex++) {
                 h = max(cbrt(1.0e-33), cbrt(1.0e-33)*fabs(singleGrainDisplacements[boundaryNodeId](dofIndex)));
                 singleGrainDisplacements[boundaryNodeId](dofIndex) += h;
                 singleContactForcePerturbed = _grains[gindex].assembleContactForceVectorRigidStraightBoundary(singleGrainDisplacements, worldDisplacementsAccumulate[gindex],
                   updatecheck, _walls[windex], singleWallDisplacements, 
                   worldWallDisplacementsAccumulate[windex],boundaryNodeIdCheckTemp,penaltyratio);
                 singleContactForcePerturbed = (singleContactForcePerturbed-singleContactForce)/h;
                 for (size_t cindex = 0; cindex < singleContactForcePerturbed.size(); cindex++) {
                    if (singleContactForcePerturbed(cindex)!=0) {
                      worldTripletList.push_back(Triplet<double,size_t>(startingPoint+cindex,
                        startingPoint+boundaryNodeId*2+dofIndex,-1.*singleContactForcePerturbed(cindex)));                    
                    }
                 } 
                 singleGrainDisplacements[boundaryNodeId](dofIndex) -= h;
              }
          }
        }
      }
    }   
    return worldContactForceVector;
  }
    
  
  VectorXd
  assembleWorldContactForceBetweenGrains(const vector< vector<Vector2d> > & worldDisplacements, const vector< vector<Vector2d> > & worldDisplacementsAccumulate, 
    const size_t & updateFlag,  vector<Triplet<double,size_t> > & worldTripletList, const double & penaltyratio) {
    VectorXd worldContactForceVector;
    worldContactForceVector.resize(_totalNodeNumber*2);
    worldContactForceVector.fill(0.);

    worldTripletList.clear();
    worldTripletList.reserve(_totalNodeNumber);


    size_t startingPointThis = 0;
    size_t startingPointOther = _grains[0].getNumberOfNodes()*2;
    size_t blockSizeThis = 0;
    size_t blockSizeOther = 0;

    vector<Vector2d> displacementsThis;
    vector<Vector2d> displacementsOther;

    size_t updatecheck = 0;
    size_t boundaryNodeIdThis = 0;
    size_t boundaryNodeIdOther = 0;
    double h = 0.;


    vector<VectorXd> contactForceUsingThisAsMaster;
    vector<VectorXd> contactForceUsingOtherAsMaster;
    vector<VectorXd> contactForceUsingThisAsMasterPerturbed;
    vector<VectorXd> contactForceUsingOtherAsMasterPerturbed;

    vector<size_t>   boundaryNodeIdsCheckThis;
    vector<size_t>   boundaryNodeIdsCheckOther;
    vector<size_t>   boundaryCheckThisTemp;
    vector<size_t>   boundaryCheckOtherTemp;

    VectorXd contactForceThis, contactForceOther;
    VectorXd contactForceThisPerturbed, contactForceOtherPerturbed;

    vector<size_t>   boundaryNodeIdsThis;
    vector<size_t>   boundaryNodeIdsOther;

    for (unsigned int i = 0; i < _ngrains; i++) {
      displacementsThis = worldDisplacements[i];
      blockSizeThis = _blockSizes[i]*2;
      contactForceThis.resize(blockSizeThis);
      contactForceThis.fill(0); 
      contactForceThisPerturbed.resize(blockSizeThis);
      contactForceThisPerturbed.fill(0);
      startingPointThis = accumulate(_blockSizes.begin(),_blockSizes.begin()+i,0)*2;
      for (unsigned int j = i+1; j < _ngrains; j++) {
        displacementsOther = worldDisplacements[j];
        blockSizeOther = _blockSizes[j]*2;
        contactForceOther.resize(blockSizeOther);
        contactForceOther.fill(0);
        contactForceOtherPerturbed.resize(blockSizeOther);
        contactForceOtherPerturbed.fill(0);
        startingPointOther = accumulate(_blockSizes.begin(),_blockSizes.begin()+j,0)*2;

        if (_grains[i].bCircleCheckWithGrain(_grains[j])) {


           contactForceUsingThisAsMaster = _grains[i].assembleContactForceVectorWithAnotherGrain(displacementsThis,
             _grains[j],displacementsOther,worldDisplacementsAccumulate[i], worldDisplacementsAccumulate[j], updateFlag,boundaryNodeIdsCheckThis,penaltyratio);

           contactForceUsingOtherAsMaster = _grains[j].assembleContactForceVectorWithAnotherGrain(displacementsOther,
             _grains[i],displacementsThis,worldDisplacementsAccumulate[j], worldDisplacementsAccumulate[i],updateFlag,boundaryNodeIdsCheckOther,penaltyratio);

           contactForceThis = (contactForceUsingThisAsMaster[0]+contactForceUsingOtherAsMaster[1])/2.;
           contactForceOther = (contactForceUsingThisAsMaster[1]+contactForceUsingOtherAsMaster[0])/2.;

           worldContactForceVector.block(startingPointThis,0,blockSizeThis,1) += contactForceThis;
           worldContactForceVector.block(startingPointOther,0,blockSizeOther,1) += contactForceOther;

           for (unsigned int nindex = 0; nindex < boundaryNodeIdsCheckThis.size(); nindex++) {
               boundaryNodeIdThis = boundaryNodeIdsCheckThis[nindex];
               for (unsigned int dofIndex = 0; dofIndex < 2; dofIndex++) {
                   h = max(cbrt(1.0e-33), cbrt(1.0e-33)*fabs(displacementsThis[boundaryNodeIdThis](dofIndex)));
                   displacementsThis[boundaryNodeIdThis](dofIndex) += h;
                   contactForceUsingThisAsMasterPerturbed = _grains[i].assembleContactForceVectorWithAnotherGrain(displacementsThis,
                       _grains[j],displacementsOther,worldDisplacementsAccumulate[i], worldDisplacementsAccumulate[j],updatecheck,boundaryCheckThisTemp,penaltyratio);
                   
                   contactForceUsingOtherAsMasterPerturbed = _grains[j].assembleContactForceVectorWithAnotherGrain(displacementsOther,
                       _grains[i],displacementsThis,worldDisplacementsAccumulate[j], worldDisplacementsAccumulate[i], updatecheck,boundaryCheckOtherTemp,penaltyratio);


                   contactForceThisPerturbed = (contactForceUsingThisAsMasterPerturbed[0]+contactForceUsingOtherAsMasterPerturbed[1])/2.;
                   contactForceOtherPerturbed = (contactForceUsingThisAsMasterPerturbed[1]+contactForceUsingOtherAsMasterPerturbed[0])/2.;

                   contactForceThisPerturbed = (contactForceThisPerturbed-contactForceThis)/h;
                   contactForceOtherPerturbed = (contactForceOtherPerturbed-contactForceOther)/h;
                   for (size_t cindex = 0; cindex < blockSizeThis; cindex++) {
                      if (contactForceThisPerturbed(cindex)!=0) {
                        worldTripletList.push_back(Triplet<double,size_t>(startingPointThis+cindex,
                          startingPointThis+boundaryNodeIdThis*2+dofIndex,-1.*contactForceThisPerturbed(cindex)));

                      }
                   }
                   for (size_t cindex = 0; cindex < blockSizeOther; cindex++) {
                      if (contactForceOtherPerturbed(cindex)!=0) {
                        worldTripletList.push_back(Triplet<double,size_t>(startingPointOther+cindex,
                          startingPointThis+boundaryNodeIdThis*2+dofIndex,-1.*contactForceOtherPerturbed(cindex)));      
                      }                      
                   }
                   displacementsThis[boundaryNodeIdThis](dofIndex) -= h;
                }
           }

           for (unsigned int nindex = 0; nindex < boundaryNodeIdsCheckOther.size(); nindex++) {
               boundaryNodeIdOther = boundaryNodeIdsCheckOther[nindex];
               for (unsigned int dofIndex = 0; dofIndex < 2; dofIndex++) {
                   h = max(cbrt(1.0e-33), cbrt(1.0e-33)*fabs(displacementsOther[boundaryNodeIdOther](dofIndex)));
                   displacementsOther[boundaryNodeIdOther](dofIndex) += h;
                   contactForceUsingThisAsMasterPerturbed = _grains[i].assembleContactForceVectorWithAnotherGrain(displacementsThis,
                       _grains[j],displacementsOther,worldDisplacementsAccumulate[i], worldDisplacementsAccumulate[j],updatecheck,boundaryCheckThisTemp,penaltyratio);
                   contactForceUsingOtherAsMasterPerturbed = _grains[j].assembleContactForceVectorWithAnotherGrain(displacementsOther,
                       _grains[i],displacementsThis,worldDisplacementsAccumulate[j], worldDisplacementsAccumulate[i],updatecheck,boundaryCheckOtherTemp,penaltyratio);

                   contactForceThisPerturbed = (contactForceUsingThisAsMasterPerturbed[0]+contactForceUsingOtherAsMasterPerturbed[1])/2.;
                   contactForceOtherPerturbed = (contactForceUsingThisAsMasterPerturbed[1]+contactForceUsingOtherAsMasterPerturbed[0])/2.;

                   contactForceThisPerturbed = (contactForceThisPerturbed-contactForceThis)/h;
                   contactForceOtherPerturbed = (contactForceOtherPerturbed-contactForceOther)/h;
                   for (size_t cindex = 0; cindex < blockSizeThis; cindex++) {
                      if (contactForceThisPerturbed(cindex)!=0) {
                        worldTripletList.push_back(Triplet<double,size_t>(startingPointThis+cindex,
                          startingPointOther+boundaryNodeIdOther*2+dofIndex,-1.*contactForceThisPerturbed(cindex)));                     
                      }
                   }
                   for (size_t cindex = 0; cindex < blockSizeOther; cindex++) {
                      if (contactForceOtherPerturbed(cindex)!=0) {
                        worldTripletList.push_back(Triplet<double,size_t>(startingPointOther+cindex,
                          startingPointOther+boundaryNodeIdOther*2+dofIndex,-1.*contactForceOtherPerturbed(cindex)));                      
                      }                      
                   }                    
                   displacementsOther[boundaryNodeIdOther](dofIndex) -= h;
                }                           
           }
        }
      }
    }
    return worldContactForceVector;
  }
  

  vector< vector<Vector3d> >
  computeWorldElementStresses(const vector< vector<Vector2d> > & worldDisplacements) {
    vector< vector<Vector3d> > worldElementStresses(_ngrains);
    for (unsigned int gindex = 0; gindex < _ngrains; gindex++) {
      worldElementStresses[gindex] = _grains[gindex].computeElementStresses(worldDisplacements[gindex]);
    }
    return worldElementStresses;
  }

  vector< vector<Vector3d> >
  computeWorldNodalStresses(const vector< vector<Vector2d> > & worldDisplacements) {
    vector< vector<Vector3d> > worldNodalStresses(_ngrains);
    for (unsigned int gindex = 0; gindex < _ngrains; gindex++) {
      worldNodalStresses[gindex] = _grains[gindex].computeNodalStresses(worldDisplacements[gindex]);
    }
    return worldNodalStresses;
  }

  vector< vector<Vector3d> >
  computeWorldNodalStrains(const vector< vector<Vector2d> > & worldDisplacements) {
    vector< vector<Vector3d> > worldNodalStrains(_ngrains);
    for (unsigned int gindex = 0; gindex < _ngrains; gindex++) {
      worldNodalStrains[gindex] = _grains[gindex].computeNodalStrains(worldDisplacements[gindex]);
    }
    return worldNodalStrains;
  }

  vector< vector<Vector2d> >
  computeWorldNodalForces(const vector<vector<Vector2d> > & worldDisplacements) {
    vector <vector<Vector2d> > worldNodalForces(_ngrains);
    vector<Vector2d> singleGrainNodalForces;
    for (unsigned int gindex = 0; gindex < _ngrains; gindex++) {
       singleGrainNodalForces = Utilities::distributeGlobalVectorsToLocalVectors(_grains[gindex].assembleForceVector(worldDisplacements[gindex]));
       worldNodalForces[gindex] = singleGrainNodalForces;
    }
    return worldNodalForces;
  }

  void
  updateWorldGrainAndWallPositions(const vector< vector<Vector2d> > & worldDisplacements, const vector< vector<Vector2d> > & worldWallDisplacements) {
  	for (unsigned int gindex = 0; gindex < _ngrains; gindex++) {
  		_grains[gindex].updateNodePositions(worldDisplacements[gindex]);
  	}
    for (unsigned int windex = 0; windex < _walls.size(); windex++) {
      _walls[windex].updateWallPositions(worldWallDisplacements[windex]);
    }
  }


  void
  moveOneRigidStraightWall(const size_t & id, const vector<Vector2d> & displacements) {
  	_walls[id].updateWallPositions(displacements);
  }

  void 
  eraseBoundaryStressRecord() {
    for (unsigned int gindex = 0; gindex < _ngrains; gindex++) {
      _grains[gindex].eraseNodalNormalAndShear();
    }
  }

  void
  changeWallFriction (const vector<double> & wallmu) {
    for (unsigned int index = 0; index < _walls.size(); index++) {
      _walls[index].changeMu(wallmu[index]);
    }
  }

  void
  changeAllGrainFriction(const double & mu) {
    for (unsigned int index = 0; index < _ngrains; index++) {
      _grains[index].changeFrictionCoefficient(mu);
    }
  }

  void
  changeSingleGrainFriction(const double & mu, const size_t & id) {
      _grains[id].changeFrictionCoefficient(mu);  
  }


  size_t
  getNumberOfGrains() const {
  	return _ngrains;
  }

  vector<RigidStraightWall> 
  getWalls() const {
  	return _walls;
  }

  vector<StaticGrain>
  getGrains() const {
  	return _grains;
  }

  vector<size_t>
  getBlockSizes() const {
  	return _blockSizes;
  }

  size_t
  getTotalNodeNumber() const {
    return _totalNodeNumber;
  }


private:
  vector<StaticGrain>           _grains;
  size_t                        _ngrains;
  size_t                        _totalNodeNumber;
  vector<RigidStraightWall>     _walls;
  vector<size_t>                _blockSizes;

};
#endif 
