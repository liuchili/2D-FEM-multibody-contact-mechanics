// -*- C++ -*-
#ifndef GRAIN_H
#define GRAIN_H

#include "Definitions.h"
#include "Utilities.h"
#include "MaterialModelLinearElastic.h"
#include "ElementSimplexTriangle.h"
#include "RigidStraightWall.h"
#include "EssentialBoundaryConditions.h"

class StaticGrain {
public:
  StaticGrain () {}
  StaticGrain (const size_t & grainId, const vector<SimplexTriangle> & elements, const size_t numberOfNodes, const vector<Vector2d> & nodePositions,
    const vector<size_t> & boundaryNodeIds, const vector< vector<size_t> > & boundaryNodeConnectivity, 
    const size_t & pnumber, const double & youngsModulus, const double & frictionCoefficient) {
    _grainId = grainId;
    _elements = elements;
    _numberOfNodes = numberOfNodes;
    _nodePositions = nodePositions;
    _boundaryNodeIds = boundaryNodeIds;
    _pnumber = pnumber;
    _boundaryNodeConnectivity = boundaryNodeConnectivity;
    _youngsModulus = youngsModulus;
    _frictionCoefficient = frictionCoefficient;


    Vector2d centricPoint; centricPoint.fill(0.);
    double totalArea = 0.;
    double singleElementArea = 0.;
    Vector2d weightedPoint; weightedPoint.fill(0.);
    vector<Vector2d> elementPositions(3);

    for (unsigned int eindex = 0; eindex < elements.size(); eindex++) {
      centricPoint.fill(0.);
      elementPositions[0] = nodePositions[elements[eindex].getNodeIds()[0]];
      elementPositions[1] = nodePositions[elements[eindex].getNodeIds()[1]];
      elementPositions[2] = nodePositions[elements[eindex].getNodeIds()[2]];
      singleElementArea = (elementPositions[1](0)-elementPositions[0](0))*(elementPositions[2](1)-elementPositions[0](1))/2.
        - (elementPositions[2](0)-elementPositions[0](0))*(elementPositions[1](1)-elementPositions[0](1))/2.;
      totalArea += singleElementArea;
      for (unsigned int pindex = 0; pindex < elementPositions.size(); pindex++) {
        centricPoint += elementPositions[pindex]/double(elementPositions.size());
      }
      weightedPoint += centricPoint*singleElementArea;
    }

    _centroid = weightedPoint/totalArea;

    double largestDistance = 0.;
    for (unsigned int nindex = 0; nindex < numberOfNodes; nindex++) {
      if ((nodePositions[nindex]-_centroid).norm() > largestDistance) {
        largestDistance = (nodePositions[nindex]-_centroid).norm();
      }
    }
    _largestRadius = largestDistance;
    
    _boundaryNodeNormalStress.resize(boundaryNodeIds.size());
    _boundaryNodeContactId.resize(boundaryNodeIds.size());
    _boundaryNodeShearStress.resize(boundaryNodeIds.size());


    _subConnectionNodeNormalStress.resize(boundaryNodeConnectivity.size());
    _subConnectionNodeShearStress.resize(boundaryNodeConnectivity.size());
    _subConnectionContactId.resize(boundaryNodeConnectivity.size());
    _subConnectionNodePositions.resize(boundaryNodeConnectivity.size());

    for (unsigned int bindex = 0; bindex < boundaryNodeIds.size(); bindex++) {
      _boundaryNodeNormalStress[bindex] = Vector2d(0.,0.);
      _boundaryNodeShearStress[bindex] = Vector2d(0.,0.);
      _boundaryNodeContactId[bindex] = 50000;
    }
    double edgeLength = 0;
    Vector2d tangentVector; tangentVector.fill(0.);
    for (unsigned int cindex = 0; cindex < boundaryNodeConnectivity.size(); cindex++) {
      _subConnectionNodeNormalStress[cindex].resize(pnumber);
      _subConnectionNodeShearStress[cindex].resize(pnumber);
      _subConnectionContactId[cindex].resize(pnumber);
      _subConnectionNodePositions[cindex].resize(pnumber);
      tangentVector = nodePositions[boundaryNodeIds[boundaryNodeConnectivity[cindex][1]]] - 
      nodePositions[boundaryNodeIds[boundaryNodeConnectivity[cindex][0]]];
      edgeLength = tangentVector.norm();
      tangentVector /= edgeLength;
      for (unsigned int ii = 0; ii < pnumber; ii++) {
        _subConnectionContactId[cindex][ii] = 50000;
        _subConnectionNodeNormalStress[cindex][ii] = Vector2d(0.,0.);
        _subConnectionNodeShearStress[cindex][ii] = Vector2d(0.,0.);
        _subConnectionNodePositions[cindex][ii] = nodePositions[boundaryNodeIds[boundaryNodeConnectivity[cindex][0]]]+tangentVector*edgeLength*(double(ii)+1.)/(double(pnumber)+1.);
      }
    }
  }

  bool
  bCircleCheckWithGrain(const StaticGrain & other) {
    bool checkFlag = false;
    if ((_centroid-other.getCentroid()).norm() < 1.01*(_largestRadius+other.getLargestRadius())) {
      checkFlag = true;
    }
    return checkFlag;
  }

  bool
  bCircleCheckWithRigidStraightWall(const RigidStraightWall & wall) {
    bool checkFlag = false;
    double distance = fabs((_centroid - wall.getBoundaryPoints()[0]).dot(wall.getWallNormal()));
    if (distance < 1.01*_largestRadius) {
      checkFlag = true;
    }
    return checkFlag;
  }

  double
  assembleEnergy(const vector<Vector2d> & displacements) const {
    double energy = 0.;
    vector<Vector2d> elementdisplacements(3);
    vector<size_t> elementNodeIds;
    for(unsigned int eindex = 0; eindex < _elements.size(); eindex++) {
      elementNodeIds = _elements[eindex].getNodeIds();
      elementdisplacements = Utilities::getElementDisplacementsFromGlobalList(elementNodeIds, displacements);
      energy += _elements[eindex].computeEnergy(elementdisplacements);
    }
    return energy;
  }

  VectorXd
  assembleForceVector(const vector<Vector2d> & displacements) const {
    size_t numberDofs = displacements.size()*2;
    VectorXd forceVector(numberDofs);
    forceVector.fill(0.);

    vector<size_t> elementNodeIds;
    vector<Vector2d> elementdisplacements(3);
    vector<Vector2d> elementForces;
    size_t nodeId;

    for(unsigned int eindex = 0; eindex < _elements.size(); eindex++) {
      elementNodeIds = _elements[eindex].getNodeIds();
      elementdisplacements = Utilities::getElementDisplacementsFromGlobalList(elementNodeIds, displacements);
      elementForces = _elements[eindex].computeForces(elementdisplacements);
      for(unsigned int nindex = 0; nindex < elementForces.size(); nindex++) {
         nodeId = elementNodeIds[nindex];
        for(unsigned int i = 0; i < 2; i++) {
           forceVector(nodeId*2+i) += elementForces[nindex](i);
        }
      }
    }
    return forceVector;
  }

  VectorXd
  assembleBodyForceVector(const Vector2d & acceleration) const {
    size_t numberDofs = _numberOfNodes*2;
    VectorXd forceVector(numberDofs);
    forceVector.fill(0.);

    vector<size_t> elementNodeIds;
    vector<Vector2d> elementForces;
    size_t nodeId;

    for(unsigned int eindex = 0; eindex < _elements.size(); eindex++) {
      elementNodeIds = _elements[eindex].getNodeIds();
      elementForces = _elements[eindex].computeBodyForces(acceleration);
      for(unsigned int nindex = 0; nindex < elementForces.size(); nindex++) {
         nodeId = elementNodeIds[nindex];
        for(unsigned int i = 0; i < 2; i++) {
           forceVector(nodeId*2+i) += elementForces[nindex](i);
        }
      }
    }
    return forceVector;
  }

  VectorXd
  assembleContactForceVectorRigidStraightBoundary( const vector<Vector2d> & displacements, const vector<Vector2d> & displacementsAccumulate, 
    const size_t & updateFlag, const RigidStraightWall & wall, const vector<Vector2d> & wallDisplacements,
    const vector<Vector2d> & wallDisplacementsAccumulate, vector<size_t> & boundaryIdGlobalCheck, const double & penaltyratio) {

    boundaryIdGlobalCheck.resize(0);
    double kn = _youngsModulus/penaltyratio;
    double kt = kn;
    double mu = wall.getMu();

    VectorXd globalContactForce; globalContactForce.resize(_numberOfNodes*2); 
    globalContactForce.fill(0);
    double distance1, disance2;
    Vector2d tangent;
    double relativeTagentialDisp;
    Vector2d nodePosition1, nodePosition2;

    vector<bool>     penetrationFlag(_boundaryNodeIds.size());
    vector<bool>     connectionPenetrationFlag(_boundaryNodeConnectivity.size());


    vector<Vector2d> boundaryNormalIncrementUpdated(_boundaryNodeIds.size()); 
    vector<Vector2d> boundaryShearIncrementUpdated(_boundaryNodeIds.size()); 


    Vector2d shearIncrement; shearIncrement.fill(0.);
    Vector2d normalIncrement; normalIncrement.fill(0.);

    vector<Vector2d> boundaryPoints = wall.getBoundaryPoints();
    size_t boundaryId = wall.getWallId();

    Vector2d normal; normal.fill(0.);
    Vector2d firstPoint = boundaryPoints[0]+wallDisplacements[0];
    Vector2d secondPoint = boundaryPoints[1]+wallDisplacements[1];
    tangent = secondPoint - firstPoint;
    tangent /= tangent.norm();
    normal << -tangent(1), tangent(0);

    double ratio = 0.;

    
    for (unsigned int bNodeIndex = 0; bNodeIndex < _boundaryNodeIds.size(); bNodeIndex++) {
      normalIncrement.fill(0.);
      shearIncrement.fill(0.);
      penetrationFlag[bNodeIndex] = false;
      nodePosition1 = _nodePositions[_boundaryNodeIds[bNodeIndex]]+displacements[_boundaryNodeIds[bNodeIndex]];
      distance1 = (nodePosition1-firstPoint).dot(normal);
      boundaryNormalIncrementUpdated[bNodeIndex].fill(0.);
      boundaryShearIncrementUpdated[bNodeIndex].fill(0.);
      if (distance1 <= max(0.2,_largestRadius/20.)) {
        boundaryIdGlobalCheck.push_back(_boundaryNodeIds[bNodeIndex]);
      }
      if (distance1 < 0) {
        penetrationFlag[bNodeIndex] = true;
        ratio = (nodePosition1-firstPoint).dot(tangent)/(firstPoint-secondPoint).norm();
        normalIncrement = kn*fabs(distance1)*normal;
        relativeTagentialDisp = (displacements[_boundaryNodeIds[bNodeIndex]]+displacementsAccumulate[_boundaryNodeIds[bNodeIndex]]
          -(1.-ratio)*(wallDisplacements[0]+wallDisplacementsAccumulate[0])-ratio*(wallDisplacements[1]+wallDisplacementsAccumulate[0])).dot(tangent);
        shearIncrement = kt*relativeTagentialDisp*(-1.*tangent);
        boundaryNormalIncrementUpdated[bNodeIndex] = normalIncrement;
        boundaryShearIncrementUpdated[bNodeIndex] = min(shearIncrement.norm(),mu*normalIncrement.norm())*shearIncrement/(shearIncrement.norm()+DBL_MIN);
      }
    }


    size_t firstNodeIdLocal, secondNodeIdLocal;
    size_t firstNodeIdGlobal, secondNodeIdGlobal;
    vector<double>     singleConnectionShapeFunctionFirstNode(_pnumber+2);
    vector<Vector2d>   singleConnectionDeformedPositions(_pnumber+2);

    vector<Vector2d>   singleConnectionNormalUpdatedIncrements(_pnumber+2);
    vector<Vector2d>   singleConnectionShearUpdatedIncrements(_pnumber+2);

    vector<bool>       subPointPenetrationFlags(_pnumber);

    Vector2d singleConnectionFirstNodeForce; singleConnectionFirstNodeForce.fill(0.);
    Vector2d singleConnectionSecondNodeForce; singleConnectionSecondNodeForce.fill(0.);
    Vector2d startTotalStress; startTotalStress.fill(0.);
    Vector2d endTotalStress; endTotalStress.fill(0.);
    double subEdgeLength = 0.;
    double wallRatio = 0.;

    for (unsigned int bcIndex = 0; bcIndex < _boundaryNodeConnectivity.size(); bcIndex++) {
      firstNodeIdLocal = _boundaryNodeConnectivity[bcIndex][0];
      secondNodeIdLocal = _boundaryNodeConnectivity[bcIndex][1];
      firstNodeIdGlobal = _boundaryNodeIds[firstNodeIdLocal];
      secondNodeIdGlobal = _boundaryNodeIds[secondNodeIdLocal];
      connectionPenetrationFlag[bcIndex] = false;
      if (penetrationFlag[firstNodeIdLocal] || penetrationFlag[secondNodeIdLocal]) {
        singleConnectionFirstNodeForce.fill(0.);
        singleConnectionSecondNodeForce.fill(0.);

        singleConnectionShapeFunctionFirstNode[0] = 1.;
        singleConnectionShapeFunctionFirstNode[_pnumber+1] = 0.;
        singleConnectionDeformedPositions[0] = _nodePositions[firstNodeIdGlobal]+displacements[firstNodeIdGlobal];
        singleConnectionDeformedPositions[_pnumber+1] = _nodePositions[secondNodeIdGlobal]+displacements[secondNodeIdGlobal];

        singleConnectionNormalUpdatedIncrements[0] = boundaryNormalIncrementUpdated[firstNodeIdLocal];
        singleConnectionShearUpdatedIncrements[0] = boundaryShearIncrementUpdated[firstNodeIdLocal];
        singleConnectionNormalUpdatedIncrements[_pnumber+1] = boundaryNormalIncrementUpdated[secondNodeIdLocal];
        singleConnectionShearUpdatedIncrements[_pnumber+1] = boundaryShearIncrementUpdated[secondNodeIdLocal];

        for (unsigned int subIndex = 0; subIndex < _pnumber; subIndex++) {
          ratio = (_subConnectionNodePositions[bcIndex][subIndex]-_nodePositions[firstNodeIdGlobal]).norm()/
            (_nodePositions[firstNodeIdGlobal]-_nodePositions[secondNodeIdGlobal]).norm();
          singleConnectionDeformedPositions[subIndex+1] = _subConnectionNodePositions[bcIndex][subIndex]
            +(1.-ratio)*displacements[firstNodeIdGlobal]+ratio*displacements[secondNodeIdGlobal];
          distance1 = (singleConnectionDeformedPositions[subIndex+1]-firstPoint).dot(normal);
          subPointPenetrationFlags[subIndex] = false;
          singleConnectionShapeFunctionFirstNode[subIndex+1] = (singleConnectionDeformedPositions[subIndex+1]-singleConnectionDeformedPositions[_pnumber+1]).norm()/
            (singleConnectionDeformedPositions[0]-singleConnectionDeformedPositions[_pnumber+1]).norm();
          singleConnectionNormalUpdatedIncrements[subIndex+1].fill(0.);
          singleConnectionShearUpdatedIncrements[subIndex+1].fill(0.);
          if (distance1 < 0) {
            subPointPenetrationFlags[subIndex] = true;
            connectionPenetrationFlag[bcIndex] = true;
            wallRatio = (singleConnectionDeformedPositions[subIndex+1]-firstPoint).dot(tangent)/(firstPoint-secondPoint).norm();
            normalIncrement = kn*fabs(distance1)*normal;
            relativeTagentialDisp = ((1.-ratio)*(displacements[firstNodeIdGlobal]+displacementsAccumulate[firstNodeIdGlobal])
              +ratio*(displacements[secondNodeIdGlobal]+displacementsAccumulate[secondNodeIdGlobal])
              -(1.-wallRatio)*(wallDisplacements[0]+wallDisplacementsAccumulate[0])-wallRatio*(wallDisplacements[1]+wallDisplacementsAccumulate[1])).dot(tangent);
            shearIncrement = kt*relativeTagentialDisp*(-1.*tangent);
            singleConnectionNormalUpdatedIncrements[subIndex+1] = normalIncrement;
            singleConnectionShearUpdatedIncrements[subIndex+1] = min(shearIncrement.norm(),mu*normalIncrement.norm())*shearIncrement/(shearIncrement.norm()+DBL_MIN);                                  
          }
        }
        for (unsigned int summationIndex = 0; summationIndex < _pnumber+1; summationIndex++) {

          subEdgeLength = (singleConnectionDeformedPositions[summationIndex]-singleConnectionDeformedPositions[summationIndex+1]).norm();

          startTotalStress = (singleConnectionNormalUpdatedIncrements[summationIndex]+singleConnectionShearUpdatedIncrements[summationIndex])
                              *singleConnectionShapeFunctionFirstNode[summationIndex];
          endTotalStress = (singleConnectionNormalUpdatedIncrements[summationIndex+1]+singleConnectionShearUpdatedIncrements[summationIndex+1])
                              *singleConnectionShapeFunctionFirstNode[summationIndex+1];            
          singleConnectionFirstNodeForce += (startTotalStress+endTotalStress)*subEdgeLength*0.5;

          startTotalStress = (singleConnectionNormalUpdatedIncrements[summationIndex]+singleConnectionShearUpdatedIncrements[summationIndex])
                              *(1.-singleConnectionShapeFunctionFirstNode[summationIndex]);
          endTotalStress = (singleConnectionNormalUpdatedIncrements[summationIndex+1]+singleConnectionShearUpdatedIncrements[summationIndex+1])
                              *(1.-singleConnectionShapeFunctionFirstNode[summationIndex+1]);
          singleConnectionSecondNodeForce += (startTotalStress+endTotalStress)*subEdgeLength*0.5;                               
        }
        globalContactForce(firstNodeIdGlobal*2) += singleConnectionFirstNodeForce(0);
        globalContactForce(firstNodeIdGlobal*2+1) += singleConnectionFirstNodeForce(1);
        globalContactForce(secondNodeIdGlobal*2) += singleConnectionSecondNodeForce(0);
        globalContactForce(secondNodeIdGlobal*2+1) += singleConnectionSecondNodeForce(1);
        if (updateFlag == 1) {
          for (unsigned int sindex = 0; sindex < _pnumber; sindex++) {
            if (subPointPenetrationFlags[sindex]) {
              _subConnectionNodeNormalStress[bcIndex][sindex] = singleConnectionNormalUpdatedIncrements[sindex+1];
              _subConnectionNodeShearStress[bcIndex][sindex] = singleConnectionShearUpdatedIncrements[sindex+1];
              _subConnectionContactId[bcIndex][sindex] = boundaryId;
            }
            else if (_subConnectionContactId[bcIndex][sindex] == boundaryId) {
              _subConnectionContactId[bcIndex][sindex] = 50000;
              _subConnectionNodeNormalStress[bcIndex][sindex].fill(0.);
              _subConnectionNodeShearStress[bcIndex][sindex].fill(0.);
            }
          }
        }        
      }
    }
    
    if (updateFlag == 1) {
      size_t bcIndexFirst = 0;
      size_t bcIndexSecond = 0; 
      for (unsigned int bcIndex = 0; bcIndex < _boundaryNodeConnectivity.size(); bcIndex++) {
        firstNodeIdLocal = _boundaryNodeConnectivity[bcIndex][0];
        if (bcIndex == 0) {
           bcIndexSecond = bcIndex;
           bcIndexFirst = _boundaryNodeConnectivity.size()-1;
        }
        else {
           bcIndexSecond = bcIndex;
           bcIndexFirst = bcIndex-1;
        }
        if (penetrationFlag[firstNodeIdLocal]) {
          _boundaryNodeNormalStress[firstNodeIdLocal] = boundaryNormalIncrementUpdated[firstNodeIdLocal];
          _boundaryNodeShearStress[firstNodeIdLocal] = boundaryShearIncrementUpdated[firstNodeIdLocal];
          _boundaryNodeContactId[firstNodeIdLocal] = boundaryId;
        }
        else if (connectionPenetrationFlag[bcIndexFirst] || connectionPenetrationFlag[bcIndexSecond]) {
          _boundaryNodeNormalStress[firstNodeIdLocal].fill(0.);
          _boundaryNodeShearStress[firstNodeIdLocal].fill(0.);
          _boundaryNodeContactId[firstNodeIdLocal] = boundaryId;
        }
        else if (_boundaryNodeContactId[firstNodeIdLocal] == boundaryId) {
          _boundaryNodeNormalStress[firstNodeIdLocal].fill(0.);
          _boundaryNodeShearStress[firstNodeIdLocal].fill(0.);
          _boundaryNodeContactId[firstNodeIdLocal] = 50000;
        }

      }

    }
    return globalContactForce;
  }
 
  
  vector<VectorXd>
  assembleContactForceVectorWithAnotherGrain(const vector<Vector2d> & thisDisplacements,
    StaticGrain & other, const vector<Vector2d> & otherDisplacements,
    const vector<Vector2d> & thisDisplacementsAccumulate, const vector<Vector2d> & otherDisplacementsAccumulate,
    const size_t & updateFlag, vector<size_t> & boundaryIdGlobalCheck, const double & penaltyratio) {

    double kn = min(_youngsModulus, other.getYoungsModulus())/penaltyratio;
    double kt = kn;
    double mu = min(_frictionCoefficient, other.getFrictionCoefficient());

    boundaryIdGlobalCheck.resize(0);

    vector<double>     projectionRatioFromOther(_boundaryNodeIds.size()); 
    vector<size_t>     projectionIdsFromOther(_boundaryNodeIds.size()); 
    vector<size_t>     projectionConnectionIdsFromOther(_boundaryNodeIds.size());
    vector<size_t>     projectionFlagsFromOther(_boundaryNodeIds.size()); 
    vector<bool>       penetrationFlagsFromOther(_boundaryNodeIds.size());

    vector<bool>       connectionPenetrationFlag(_boundaryNodeConnectivity.size());


    vector<size_t>     projectionIdsAlternative(_boundaryNodeIds.size());
    vector<Vector2d>   normalVectorsAlternative(_boundaryNodeIds.size());
    vector<Vector2d>   tangentialVectorsAlternative(_boundaryNodeIds.size());
    vector<double>     ratiosAlternative(_boundaryNodeIds.size());


    vector<Vector2d>   boundaryNormalIncrementUpdated(_boundaryNodeIds.size());
    vector<Vector2d>   boundaryShearIncrementUpdated(_boundaryNodeIds.size());

    const vector<size_t> otherBoundaryNodeIds = other.getBoundaryNodeId();
    const vector< vector<size_t> > otherboundaryConnections = other.getBoundaryNodeConnectivity();

    const double thresDistance = 0.1;
    const double thresPortion = 0.001;

    vector<Vector2d> otherBoundaryPositions(otherBoundaryNodeIds.size());
    vector<Vector2d> otherBoundaryDisplacements(otherBoundaryNodeIds.size());
    for (unsigned int index = 0; index < otherBoundaryNodeIds.size(); index++) {
      otherBoundaryPositions[index] = other.getNodePositions()[otherBoundaryNodeIds[index]];
      otherBoundaryDisplacements[index] = otherDisplacements[otherBoundaryNodeIds[index]];
    }
    vector <size_t> boundaryConnectionIndex; boundaryConnectionIndex.resize(otherboundaryConnections.size());
    for (unsigned int index = 0; index < boundaryConnectionIndex.size(); index++) {
      boundaryConnectionIndex[index] = index;
    }

    size_t otherConnectionIndexFirst = 0;
    size_t otherConnectionIndexSecond = 0;


    Vector2d singleNodePosition; singleNodePosition.fill(0.);
    Vector2d singleNodeDisplacement; singleNodeDisplacement.fill(0.);
    Vector2d singleNodeDisplacementAccumulate; singleNodeDisplacementAccumulate.fill(0.);
    size_t singleNodeIdGlobal = 0;


    vector<Vector2d> threeOtherNodePositions(3);
    vector<Vector2d> threeOtherNodeDisplacements(3);
    vector<Vector2d> threeOhterNodeDisplacementsAccumulate(3);

    vector<size_t>   threeOtherBoundaryNodeIds(3);
    vector<size_t>   twoOtherBoundaryConnectionIds(2); 

    size_t           projectionFlag = 0;
    size_t           projectionId = 0;
    size_t           connectionId = 0;
    double           ratio = 0.;
    bool             penetrationFlag = false;
    double           normalDistanceTemp = 0.;
    Vector2d         normal; normal.fill(0.);
    Vector2d         tangentialDistance; tangentialDistance.fill(0.);


    size_t           projectionIdAlternative = 0;
    double           ratioAlternative = 0.;
    Vector2d         normalVectorAlternative;
    Vector2d         tangentialVectorAlternative;


    Vector2d shearIncrement; shearIncrement.fill(0.);
    Vector2d normalIncrement; normalIncrement.fill(0.);


    size_t otherNodeIdFirstLocal = 0;
    size_t otherNodeIdSecondLocal = 0;
    size_t otherNodeIdFirstGlobal = 0;
    size_t otherNodeIdSecondGlobal = 0;

    for (unsigned int bNodeIndex = 0; bNodeIndex < _boundaryNodeIds.size(); bNodeIndex++) {
      singleNodeIdGlobal = _boundaryNodeIds[bNodeIndex];
      singleNodePosition = _nodePositions[singleNodeIdGlobal];
      singleNodeDisplacement = thisDisplacements[singleNodeIdGlobal];
      singleNodeDisplacementAccumulate = thisDisplacementsAccumulate[singleNodeIdGlobal];

      shearIncrement.fill(0.);
      normalIncrement.fill(0.);
      boundaryNormalIncrementUpdated[bNodeIndex].fill(0.);
      boundaryShearIncrementUpdated[bNodeIndex].fill(0.);
      projectionRatioFromOther[bNodeIndex] = 0.;
      projectionIdsFromOther[bNodeIndex] = _numberOfNodes+1;
      projectionConnectionIdsFromOther[bNodeIndex] = 0;
      projectionFlagsFromOther[bNodeIndex] = 0;
      penetrationFlagsFromOther[bNodeIndex] = false;

      projectionIdsAlternative[bNodeIndex] = otherboundaryConnections.size()+1;
      normalVectorsAlternative[bNodeIndex].fill(0.);
      tangentialVectorsAlternative[bNodeIndex].fill(0.);
      ratiosAlternative[bNodeIndex] = 0.;

      if (Utilities::identifyOnePointProjectionToLineSegments(thresDistance,singleNodePosition, singleNodeDisplacement,otherboundaryConnections,
           boundaryConnectionIndex, otherBoundaryPositions, otherBoundaryDisplacements, threeOtherBoundaryNodeIds, twoOtherBoundaryConnectionIds)) {
         threeOtherNodePositions[0] = otherBoundaryPositions[threeOtherBoundaryNodeIds[0]];
         threeOtherNodePositions[1] = otherBoundaryPositions[threeOtherBoundaryNodeIds[1]];
         threeOtherNodePositions[2] = otherBoundaryPositions[threeOtherBoundaryNodeIds[2]];         
         threeOtherBoundaryNodeIds[0] = otherBoundaryNodeIds[threeOtherBoundaryNodeIds[0]];
         threeOtherBoundaryNodeIds[1] = otherBoundaryNodeIds[threeOtherBoundaryNodeIds[1]];
         threeOtherBoundaryNodeIds[2] = otherBoundaryNodeIds[threeOtherBoundaryNodeIds[2]];
         threeOtherNodeDisplacements[0] = otherDisplacements[threeOtherBoundaryNodeIds[0]];
         threeOtherNodeDisplacements[1] = otherDisplacements[threeOtherBoundaryNodeIds[1]];
         threeOtherNodeDisplacements[2] = otherDisplacements[threeOtherBoundaryNodeIds[2]];
         threeOhterNodeDisplacementsAccumulate[0] = otherDisplacementsAccumulate[threeOtherBoundaryNodeIds[0]];
         threeOhterNodeDisplacementsAccumulate[1] = otherDisplacementsAccumulate[threeOtherBoundaryNodeIds[1]];
         threeOhterNodeDisplacementsAccumulate[2] = otherDisplacementsAccumulate[threeOtherBoundaryNodeIds[2]];

         if (Utilities::computePointToSegments(singleNodePosition, singleNodeDisplacement, threeOtherNodePositions,
              threeOtherNodeDisplacements, thresDistance, thresPortion, singleNodeDisplacementAccumulate,threeOhterNodeDisplacementsAccumulate,
              projectionFlag, projectionId, ratio, penetrationFlag, connectionId,normalDistanceTemp, normal, tangentialDistance, 
              projectionIdAlternative, ratioAlternative, normalVectorAlternative, tangentialVectorAlternative)) {
            if (projectionFlag > 0) {
                projectionRatioFromOther[bNodeIndex] = ratio;
                if (projectionFlag == 1) {projectionIdsFromOther[bNodeIndex] = twoOtherBoundaryConnectionIds[projectionId];
                                          projectionFlagsFromOther[bNodeIndex] = 1;}
                                          
                if (projectionFlag == 2) {projectionIdsFromOther[bNodeIndex] = threeOtherBoundaryNodeIds[projectionId];
                                          projectionConnectionIdsFromOther[bNodeIndex] = twoOtherBoundaryConnectionIds[connectionId];
                                          projectionFlagsFromOther[bNodeIndex] = 2;}
                if (penetrationFlag) {
                    normalIncrement = fabs(normalDistanceTemp)*normal*kn;
                    shearIncrement = kt*tangentialDistance;
                    penetrationFlagsFromOther[bNodeIndex] = true;

                    if (normalVectorAlternative.norm() > 0) {
                       ratiosAlternative[bNodeIndex] = ratioAlternative;
                       normalVectorsAlternative[bNodeIndex] = normalVectorAlternative;
                       tangentialVectorsAlternative[bNodeIndex] = tangentialVectorAlternative;
                       projectionIdsAlternative[bNodeIndex] = twoOtherBoundaryConnectionIds[projectionIdAlternative];
                    }                   
                }               
            }
         }
      }
      if (penetrationFlagsFromOther[bNodeIndex]) {
          shearIncrement = min(shearIncrement.norm(),mu*normalIncrement.norm())*shearIncrement/(shearIncrement.norm()+DBL_MIN);
          boundaryShearIncrementUpdated[bNodeIndex] = shearIncrement;
          boundaryNormalIncrementUpdated[bNodeIndex] = normalIncrement;
      }
      if (projectionFlagsFromOther[bNodeIndex] > 0) {boundaryIdGlobalCheck.push_back(_boundaryNodeIds[bNodeIndex]);}     
    }
    
    size_t nodeIdCurrentLocal, nodeIdForwardLocal, nodeIdBackwardLocal;
    size_t nodeIdFurtherForwardLocal, nodeIdFurtherBackwardLocal;
    Vector2d normalCandidate1, normalCandidate2;
    size_t nodeIdCurrentGlobal;
    for (unsigned int bcIndex = 0; bcIndex < _boundaryNodeConnectivity.size(); bcIndex++) {
      nodeIdCurrentLocal = _boundaryNodeConnectivity[bcIndex][0];
      nodeIdForwardLocal = _boundaryNodeConnectivity[bcIndex][1];
      nodeIdCurrentGlobal = _boundaryNodeIds[nodeIdCurrentLocal];
      if (bcIndex == 0) {
        nodeIdBackwardLocal = _boundaryNodeConnectivity[_boundaryNodeConnectivity.size()-1][0];
      }
      else {
        nodeIdBackwardLocal = _boundaryNodeConnectivity[bcIndex-1][0];
      }
      if (bcIndex < 2) {
        nodeIdFurtherBackwardLocal = _boundaryNodeConnectivity[_boundaryNodeConnectivity.size()-2+bcIndex][0];
      }
      else {
        nodeIdFurtherBackwardLocal = _boundaryNodeConnectivity[bcIndex-2][0];
      }
      if (bcIndex > _boundaryNodeConnectivity.size()-3) {
        nodeIdFurtherForwardLocal = _boundaryNodeConnectivity[bcIndex+2-_boundaryNodeConnectivity.size()][0];
      }
      else {
        nodeIdFurtherForwardLocal = _boundaryNodeConnectivity[bcIndex+2][0];
      }

      if (normalVectorsAlternative[nodeIdCurrentLocal].norm()>0) {
        if (penetrationFlagsFromOther[nodeIdForwardLocal] && (!penetrationFlagsFromOther[nodeIdBackwardLocal])) {
           if (penetrationFlagsFromOther[nodeIdFurtherForwardLocal]) { 
             normal = boundaryNormalIncrementUpdated[nodeIdForwardLocal]/boundaryNormalIncrementUpdated[nodeIdForwardLocal].norm();
             normalCandidate1 = boundaryNormalIncrementUpdated[nodeIdCurrentLocal]/(boundaryNormalIncrementUpdated[nodeIdCurrentLocal].norm()+DBL_MIN);
             normalCandidate2 = normalVectorsAlternative[nodeIdCurrentLocal]/(normalVectorsAlternative[nodeIdCurrentLocal].norm()+DBL_MIN);
             if (normalCandidate2.dot(normal) > normalCandidate1.dot(normal)) {
  	           projectionIdsFromOther[nodeIdCurrentLocal] = projectionIdsAlternative[nodeIdCurrentLocal];
  	           projectionRatioFromOther[nodeIdCurrentLocal] = ratiosAlternative[nodeIdCurrentLocal];
  	           penetrationFlagsFromOther[nodeIdCurrentLocal] = true;
  	           normalIncrement = normalVectorsAlternative[nodeIdCurrentLocal]*kn;
  	           shearIncrement = tangentialVectorsAlternative[nodeIdCurrentLocal]*kt;

    		       if (penetrationFlagsFromOther[nodeIdCurrentLocal]) {
                 boundaryShearIncrementUpdated[nodeIdCurrentLocal] = min(shearIncrement.norm(),normalIncrement.norm()*mu)*shearIncrement/(shearIncrement.norm()+DBL_MIN);
    		         boundaryNormalIncrementUpdated[nodeIdCurrentLocal] = normalIncrement;
    		       }
             }
           }
        }
        else if ((!penetrationFlagsFromOther[nodeIdForwardLocal]) && penetrationFlagsFromOther[nodeIdBackwardLocal]) {
          if (penetrationFlagsFromOther[nodeIdFurtherBackwardLocal]) {
          	normal = boundaryNormalIncrementUpdated[nodeIdBackwardLocal]/boundaryNormalIncrementUpdated[nodeIdBackwardLocal].norm();
          	normalCandidate1 = boundaryNormalIncrementUpdated[nodeIdCurrentLocal]/(boundaryNormalIncrementUpdated[nodeIdCurrentLocal].norm()+DBL_MIN);
            normalCandidate2 = normalVectorsAlternative[nodeIdCurrentLocal]/(normalVectorsAlternative[nodeIdCurrentLocal].norm()+DBL_MIN);
            if (normalCandidate2.dot(normal) > normalCandidate1.dot(normal)) {
  	          projectionIdsFromOther[nodeIdCurrentLocal] = projectionIdsAlternative[nodeIdCurrentLocal];
  	          projectionRatioFromOther[nodeIdCurrentLocal] = ratiosAlternative[nodeIdCurrentLocal];
  	          penetrationFlagsFromOther[nodeIdCurrentLocal] = true;
  	          normalIncrement = normalVectorsAlternative[nodeIdCurrentLocal]*kn;
  	          shearIncrement = tangentialVectorsAlternative[nodeIdCurrentLocal]*kt;

              if (penetrationFlagsFromOther[nodeIdCurrentLocal]) {
                boundaryShearIncrementUpdated[nodeIdCurrentLocal] = min(shearIncrement.norm(),normalIncrement.norm()*mu)*shearIncrement/(shearIncrement.norm()+DBL_MIN);
                boundaryNormalIncrementUpdated[nodeIdCurrentLocal] = normalIncrement;
              }
            }            
          }
        }
      }
    }
    
    
    size_t thisNodeIdFirstLocal = 0;
    size_t thisNodeIdSecondLocal = 0;
    size_t thisNodeIdFirstGlobal = 0;
    size_t thisNodeIdSecondGlobal = 0;

    vector<double>     singleConnectionShapeFunctionFirstNode(_pnumber+2);
    vector<Vector2d>   singleConnectionDeformedPositions(_pnumber+2);

    vector<Vector2d>   singleConnectionNormalUpdatedIncrements(_pnumber+2);
    vector<Vector2d>   singleConnectionShearUpdatedIncrements(_pnumber+2);

    vector<bool>       subPointPenetrationFlags(_pnumber);


    Vector2d singleConnectionFirstNodeForce; singleConnectionFirstNodeForce.fill(0);
    Vector2d singleConnectionSecondNodeForce; singleConnectionSecondNodeForce.fill(0);


    size_t otherConnectionIdStart = 0;
    size_t otherConnectionIdEnd = 0;
    vector<size_t> otherConnectionIdsForCheck;
    vector< vector<size_t> > otherConnectionsForCheck;
    vector<size_t> connectionPathOne;
    vector<size_t> connectionPathTwo;
    size_t numberConnection;


    Vector2d startTotalStress;
    Vector2d endTotalStress;
    double   subEdgeLength;

    VectorXd contactForceThis(_numberOfNodes*2); contactForceThis.fill(0);
    VectorXd contactForceOther(other.getNumberOfNodes()*2); contactForceOther.fill(0);

    size_t nodeIdFirstOther, nodeIdSecondOther;

    size_t outsideFlagDeformed;

    Vector2d firstNodeNormal;
    Vector2d secondNodeNormal;


    for (unsigned int connectionIndex = 0; connectionIndex < _boundaryNodeConnectivity.size(); connectionIndex++) {

        thisNodeIdFirstLocal = _boundaryNodeConnectivity[connectionIndex][0];
        thisNodeIdSecondLocal = _boundaryNodeConnectivity[connectionIndex][1];
        thisNodeIdFirstGlobal = _boundaryNodeIds[thisNodeIdFirstLocal];
        thisNodeIdSecondGlobal = _boundaryNodeIds[thisNodeIdSecondLocal];

        connectionPenetrationFlag[connectionIndex] = false;

        if (projectionFlagsFromOther[thisNodeIdFirstLocal] > 0 && projectionFlagsFromOther[thisNodeIdSecondLocal] > 0) {
           singleConnectionFirstNodeForce.fill(0.);
           singleConnectionSecondNodeForce.fill(0.);

           singleConnectionShapeFunctionFirstNode[0] = 1.;
           singleConnectionShapeFunctionFirstNode[_pnumber+1] = 0.;
           singleConnectionDeformedPositions[0] = _nodePositions[thisNodeIdFirstGlobal]+thisDisplacements[thisNodeIdFirstGlobal];
           singleConnectionDeformedPositions[_pnumber+1] = _nodePositions[thisNodeIdSecondLocal]+thisDisplacements[thisNodeIdSecondGlobal];

           singleConnectionNormalUpdatedIncrements[0] = boundaryNormalIncrementUpdated[thisNodeIdFirstLocal];
           singleConnectionShearUpdatedIncrements[0] = boundaryShearIncrementUpdated[thisNodeIdFirstLocal];
           singleConnectionNormalUpdatedIncrements[_pnumber+1] = boundaryNormalIncrementUpdated[thisNodeIdSecondLocal];
           singleConnectionShearUpdatedIncrements[_pnumber+1] = boundaryShearIncrementUpdated[thisNodeIdSecondLocal];


           if (projectionFlagsFromOther[thisNodeIdFirstLocal] == 1) {otherConnectionIdStart = projectionIdsFromOther[thisNodeIdFirstLocal];}
           else if (projectionFlagsFromOther[thisNodeIdFirstLocal] == 2) {otherConnectionIdStart = projectionConnectionIdsFromOther[thisNodeIdFirstLocal];}
            
           if (projectionFlagsFromOther[thisNodeIdSecondLocal] == 1) {otherConnectionIdEnd = projectionIdsFromOther[thisNodeIdSecondLocal];}
           else if (projectionFlagsFromOther[thisNodeIdSecondLocal] == 2) {otherConnectionIdEnd = projectionConnectionIdsFromOther[thisNodeIdSecondLocal];}

           if (otherConnectionIdStart == otherConnectionIdEnd) {

              otherConnectionIdsForCheck.resize(1);
              otherConnectionIdsForCheck[0] = otherConnectionIdStart;
           	  
            }
            else {
                numberConnection = max(otherConnectionIdStart, otherConnectionIdEnd)-min(otherConnectionIdStart, otherConnectionIdEnd)+1;
                connectionPathOne.resize(numberConnection);
                for (unsigned int ii = 0; ii < numberConnection; ii++) {
                    connectionPathOne[ii] = min(otherConnectionIdStart, otherConnectionIdEnd)+ii;
                }
               
                numberConnection = otherboundaryConnections.size()-max(otherConnectionIdStart, otherConnectionIdEnd)+min(otherConnectionIdStart, otherConnectionIdEnd)+1;
                connectionPathTwo.resize(numberConnection);
                for (unsigned int ii = 0; ii < numberConnection; ii++) {
                    connectionPathTwo[ii] = max(otherConnectionIdStart, otherConnectionIdEnd)+ii;
                    if (connectionPathTwo[ii] > otherboundaryConnections.size()-1) {
                       connectionPathTwo[ii] -= otherboundaryConnections.size();
                    }
                }
                if (connectionPathOne.size() > connectionPathTwo.size()) {otherConnectionIdsForCheck = connectionPathTwo;}
                else {otherConnectionIdsForCheck = connectionPathOne;}
            }

            otherConnectionsForCheck.resize(otherConnectionIdsForCheck.size());
            for (unsigned int index = 0; index < otherConnectionIdsForCheck.size(); index++) {
              otherConnectionsForCheck[index].resize(2);
              otherConnectionsForCheck[index][0] = otherboundaryConnections[otherConnectionIdsForCheck[index]][0];
              otherConnectionsForCheck[index][1] = otherboundaryConnections[otherConnectionIdsForCheck[index]][1];
            }

            for (unsigned int subIndex = 0; subIndex < _pnumber; subIndex++) {
                normalIncrement.fill(0.);
                shearIncrement.fill(0.);
                subPointPenetrationFlags[subIndex] = false;
                singleConnectionShearUpdatedIncrements[subIndex+1].fill(0.);
                singleConnectionNormalUpdatedIncrements[subIndex+1].fill(0.);
                singleNodePosition = _subConnectionNodePositions[connectionIndex][subIndex];
                ratio = (singleNodePosition-_nodePositions[thisNodeIdFirstGlobal]).norm()/(_nodePositions[thisNodeIdSecondGlobal]-_nodePositions[thisNodeIdFirstGlobal]).norm();
                singleNodeDisplacement = thisDisplacements[thisNodeIdFirstGlobal]*(1.-ratio)+thisDisplacements[thisNodeIdSecondGlobal]*ratio;
                singleNodeDisplacementAccumulate = thisDisplacementsAccumulate[thisNodeIdFirstGlobal]*(1.-ratio)+thisDisplacementsAccumulate[thisNodeIdSecondGlobal]*ratio;
                singleConnectionDeformedPositions[subIndex+1] = singleNodePosition+singleNodeDisplacement;
                singleConnectionShapeFunctionFirstNode[subIndex+1] = (singleNodePosition+singleNodeDisplacement-
                                                                       _nodePositions[thisNodeIdSecondGlobal]-thisDisplacements[thisNodeIdSecondGlobal]).norm();

                singleConnectionShapeFunctionFirstNode[subIndex+1] /= (_nodePositions[thisNodeIdSecondGlobal]+thisDisplacements[thisNodeIdSecondGlobal]-
                                                                        _nodePositions[thisNodeIdFirstGlobal]-thisDisplacements[thisNodeIdFirstGlobal]).norm();
                if (Utilities::identifyOnePointProjectionToLineSegments(thresDistance,singleNodePosition, singleNodeDisplacement,otherConnectionsForCheck,
                                  otherConnectionIdsForCheck, otherBoundaryPositions, otherBoundaryDisplacements, threeOtherBoundaryNodeIds, twoOtherBoundaryConnectionIds)) {
                   if (threeOtherBoundaryNodeIds.size()==3) {
                      threeOtherNodePositions[0] = otherBoundaryPositions[threeOtherBoundaryNodeIds[0]];
                      threeOtherNodePositions[1] = otherBoundaryPositions[threeOtherBoundaryNodeIds[1]];
                      threeOtherNodePositions[2] = otherBoundaryPositions[threeOtherBoundaryNodeIds[2]];                    
                      threeOtherBoundaryNodeIds[0] = otherBoundaryNodeIds[threeOtherBoundaryNodeIds[0]];
                      threeOtherBoundaryNodeIds[1] = otherBoundaryNodeIds[threeOtherBoundaryNodeIds[1]];
                      threeOtherBoundaryNodeIds[2] = otherBoundaryNodeIds[threeOtherBoundaryNodeIds[2]];
                      threeOtherNodeDisplacements[0] = otherDisplacements[threeOtherBoundaryNodeIds[0]];
                      threeOtherNodeDisplacements[1] = otherDisplacements[threeOtherBoundaryNodeIds[1]];
                      threeOtherNodeDisplacements[2] = otherDisplacements[threeOtherBoundaryNodeIds[2]];
                      threeOhterNodeDisplacementsAccumulate[0] = otherDisplacementsAccumulate[threeOtherBoundaryNodeIds[0]];
                      threeOhterNodeDisplacementsAccumulate[1] = otherDisplacementsAccumulate[threeOtherBoundaryNodeIds[1]];
                      threeOhterNodeDisplacementsAccumulate[2] = otherDisplacementsAccumulate[threeOtherBoundaryNodeIds[2]];

                   
                      if (Utilities::computePointToSegments(singleNodePosition, singleNodeDisplacement, threeOtherNodePositions,
                          threeOtherNodeDisplacements, thresDistance, thresPortion, singleNodeDisplacementAccumulate,threeOhterNodeDisplacementsAccumulate,
                          projectionFlag, projectionId, ratio, penetrationFlag, connectionId,normalDistanceTemp, normal, tangentialDistance, projectionIdAlternative, 
                          ratioAlternative, normalVectorAlternative, tangentialVectorAlternative)) {

                        if (penetrationFlag) {
                            subPointPenetrationFlags[subIndex] = true;
                            connectionPenetrationFlag[connectionIndex] = true;
                            normalIncrement = kn*fabs(normalDistanceTemp)*normal;
                            shearIncrement = kt*tangentialDistance;
                            
                            
                            if (penetrationFlagsFromOther[thisNodeIdFirstLocal] && (!penetrationFlagsFromOther[thisNodeIdSecondLocal])) {
                              if (normalVectorAlternative.norm()> 0) {
                                firstNodeNormal = boundaryNormalIncrementUpdated[thisNodeIdFirstLocal]/(boundaryNormalIncrementUpdated[thisNodeIdFirstLocal].norm()+DBL_MIN);
                                secondNodeNormal = normalVectorAlternative/normalVectorAlternative.norm();                                
                                if (firstNodeNormal.dot(normal) < firstNodeNormal.dot(secondNodeNormal)) {
                                  normalIncrement = kn*normalVectorAlternative;
                                  shearIncrement = kt*tangentialVectorAlternative;
                                }
                              }
                            }
                            else if ((!penetrationFlagsFromOther[thisNodeIdFirstLocal]) && penetrationFlagsFromOther[thisNodeIdSecondLocal]) {
                              if (normalVectorAlternative.norm()> 0) {
                                secondNodeNormal = boundaryNormalIncrementUpdated[thisNodeIdSecondLocal]/(boundaryNormalIncrementUpdated[thisNodeIdSecondLocal].norm()+DBL_MIN);
                                firstNodeNormal = normalVectorAlternative/normalVectorAlternative.norm();
                                if (secondNodeNormal.dot(normal) < secondNodeNormal.dot(firstNodeNormal)) {
                                  normalIncrement = kn*normalVectorAlternative;
                                  shearIncrement = kt*tangentialVectorAlternative;
                                }
                              }
                            }
                                                                                                        
                        }
                      }
                   }
                   else if (threeOtherBoundaryNodeIds.size()==2) {
                      threeOtherNodePositions[0] = otherBoundaryPositions[threeOtherBoundaryNodeIds[0]];
                      threeOtherNodePositions[1] = otherBoundaryPositions[threeOtherBoundaryNodeIds[1]];
                      threeOtherNodeDisplacements[0] = otherDisplacements[otherBoundaryNodeIds[threeOtherBoundaryNodeIds[0]]];
                      threeOtherNodeDisplacements[1] = otherDisplacements[otherBoundaryNodeIds[threeOtherBoundaryNodeIds[1]]];
                      threeOhterNodeDisplacementsAccumulate[0] = otherDisplacementsAccumulate[otherBoundaryNodeIds[threeOtherBoundaryNodeIds[0]]];
                      threeOhterNodeDisplacementsAccumulate[1] = otherDisplacementsAccumulate[otherBoundaryNodeIds[threeOtherBoundaryNodeIds[1]]];

                      Utilities::computePointToSingleSegment(singleNodePosition, singleNodeDisplacement,singleNodeDisplacementAccumulate, threeOtherNodePositions[0],threeOtherNodePositions[1],
                         threeOtherNodeDisplacements[0], threeOtherNodeDisplacements[1], threeOhterNodeDisplacementsAccumulate[0],threeOhterNodeDisplacementsAccumulate[1],
                         penetrationFlag, normal, tangentialDistance, tangentialVectorAlternative, 
                         outsideFlagDeformed, ratio);
                      if (penetrationFlag) {
                          subPointPenetrationFlags[subIndex] = true;
                          connectionPenetrationFlag[connectionIndex] = true;
                          normalIncrement = kn*normal;
                          shearIncrement = kt*(-1.*tangentialDistance)*tangentialVectorAlternative.dot(tangentialDistance);
                      }                      
                   }
                }
                if (subPointPenetrationFlags[subIndex]) {          
                   shearIncrement = min(shearIncrement.norm(),mu*normalIncrement.norm())*shearIncrement/(shearIncrement.norm()+DBL_MIN);
                   singleConnectionShearUpdatedIncrements[subIndex+1] = shearIncrement;
                   singleConnectionNormalUpdatedIncrements[subIndex+1] = normalIncrement;
                }
            }

            for (unsigned int summationIndex = 0; summationIndex < _pnumber+1; summationIndex++) {

                subEdgeLength = (singleConnectionDeformedPositions[summationIndex]-singleConnectionDeformedPositions[summationIndex+1]).norm();

                startTotalStress = (singleConnectionNormalUpdatedIncrements[summationIndex]+singleConnectionShearUpdatedIncrements[summationIndex])
                                    *singleConnectionShapeFunctionFirstNode[summationIndex];
                endTotalStress = (singleConnectionNormalUpdatedIncrements[summationIndex+1]+singleConnectionShearUpdatedIncrements[summationIndex+1])
                                    *singleConnectionShapeFunctionFirstNode[summationIndex+1];            
                singleConnectionFirstNodeForce += (startTotalStress+endTotalStress)*subEdgeLength*0.5;

                startTotalStress = (singleConnectionNormalUpdatedIncrements[summationIndex]+singleConnectionShearUpdatedIncrements[summationIndex])
                                    *(1.-singleConnectionShapeFunctionFirstNode[summationIndex]);
                endTotalStress = (singleConnectionNormalUpdatedIncrements[summationIndex+1]+singleConnectionShearUpdatedIncrements[summationIndex+1])
                                    *(1.-singleConnectionShapeFunctionFirstNode[summationIndex+1]);
                singleConnectionSecondNodeForce += (startTotalStress+endTotalStress)*subEdgeLength*0.5;                               
            }


            if (projectionFlagsFromOther[thisNodeIdFirstLocal]==1) {

               contactForceThis(thisNodeIdFirstGlobal*2) += singleConnectionFirstNodeForce(0);
               contactForceThis(thisNodeIdFirstGlobal*2+1) += singleConnectionFirstNodeForce(1);

               nodeIdFirstOther = otherBoundaryNodeIds[otherboundaryConnections[projectionIdsFromOther[thisNodeIdFirstLocal]][0]];
               nodeIdSecondOther = otherBoundaryNodeIds[otherboundaryConnections[projectionIdsFromOther[thisNodeIdFirstLocal]][1]];
               contactForceOther(nodeIdFirstOther*2) += singleConnectionFirstNodeForce(0)*(-1.)*(1.-projectionRatioFromOther[thisNodeIdFirstLocal]);
               contactForceOther(nodeIdFirstOther*2+1) += singleConnectionFirstNodeForce(1)*(-1.)*(1.-projectionRatioFromOther[thisNodeIdFirstLocal]);
               contactForceOther(nodeIdSecondOther*2) += singleConnectionFirstNodeForce(0)*(-1.)*projectionRatioFromOther[thisNodeIdFirstLocal];
               contactForceOther(nodeIdSecondOther*2+1) += singleConnectionFirstNodeForce(1)*(-1.)*projectionRatioFromOther[thisNodeIdFirstLocal];                   
            }
            else if (projectionFlagsFromOther[thisNodeIdFirstLocal]==2) {

               contactForceThis(thisNodeIdFirstGlobal*2) += singleConnectionFirstNodeForce(0);
               contactForceThis(thisNodeIdFirstGlobal*2+1) += singleConnectionFirstNodeForce(1);

               nodeIdFirstOther = projectionIdsFromOther[thisNodeIdFirstLocal];
               contactForceOther(nodeIdFirstOther*2) += singleConnectionFirstNodeForce(0)*(-1.);
               contactForceOther(nodeIdFirstOther*2+1) += singleConnectionFirstNodeForce(1)*(-1.);         
            }

            if (projectionFlagsFromOther[thisNodeIdSecondLocal]==1) {

               contactForceThis(thisNodeIdSecondGlobal*2) += singleConnectionSecondNodeForce(0);
               contactForceThis(thisNodeIdSecondGlobal*2+1) += singleConnectionSecondNodeForce(1);

               nodeIdFirstOther = otherBoundaryNodeIds[otherboundaryConnections[projectionIdsFromOther[thisNodeIdSecondLocal]][0]];
               nodeIdSecondOther = otherBoundaryNodeIds[otherboundaryConnections[projectionIdsFromOther[thisNodeIdSecondLocal]][1]];
               contactForceOther(nodeIdFirstOther*2) += singleConnectionSecondNodeForce(0)*(-1.)*(1.-projectionRatioFromOther[thisNodeIdSecondLocal]);
               contactForceOther(nodeIdFirstOther*2+1) += singleConnectionSecondNodeForce(1)*(-1.)*(1.-projectionRatioFromOther[thisNodeIdSecondLocal]);
               contactForceOther(nodeIdSecondOther*2) += singleConnectionSecondNodeForce(0)*(-1.)*projectionRatioFromOther[thisNodeIdSecondLocal];
               contactForceOther(nodeIdSecondOther*2+1) += singleConnectionSecondNodeForce(1)*(-1.)*projectionRatioFromOther[thisNodeIdSecondLocal];                   
            }
            else if (projectionFlagsFromOther[thisNodeIdSecondLocal]==2) {
          
               contactForceThis(thisNodeIdSecondGlobal*2) += singleConnectionSecondNodeForce(0);
               contactForceThis(thisNodeIdSecondGlobal*2+1) += singleConnectionSecondNodeForce(1);

               nodeIdFirstOther = projectionIdsFromOther[thisNodeIdSecondLocal];
               contactForceOther(nodeIdFirstOther*2) += singleConnectionSecondNodeForce(0)*(-1.);
               contactForceOther(nodeIdFirstOther*2+1) += singleConnectionSecondNodeForce(1)*(-1.);          
            }
            if (updateFlag == 1) {
               for (unsigned int sindex = 0; sindex < _pnumber; sindex++) {
                   if (subPointPenetrationFlags[sindex]) {
                      _subConnectionNodeNormalStress[connectionIndex][sindex] = singleConnectionNormalUpdatedIncrements[sindex+1];
                      _subConnectionNodeShearStress[connectionIndex][sindex] = singleConnectionShearUpdatedIncrements[sindex+1];
                      _subConnectionContactId[connectionIndex][sindex] = other.getGrainId();
                   }
                   else if (_subConnectionContactId[connectionIndex][sindex] == other.getGrainId()) {
                      _subConnectionContactId[connectionIndex][sindex] = 50000;
                      _subConnectionNodeNormalStress[connectionIndex][sindex].fill(0.);
                      _subConnectionNodeShearStress[connectionIndex][sindex].fill(0.);

                   }
               }
            }
        } 
    } 

    if (updateFlag == 1) {
       size_t bcIndexFirst = 0;
       size_t bcIndexSecond = 0;      

       for (unsigned int bcIndex = 0; bcIndex < _boundaryNodeConnectivity.size(); bcIndex++) {
         thisNodeIdFirstLocal = _boundaryNodeConnectivity[bcIndex][0];
         if (bcIndex == 0) {
            bcIndexSecond = bcIndex;
            bcIndexFirst = _boundaryNodeConnectivity.size() - 1;
         }
         else {
            bcIndexSecond = bcIndex;
            bcIndexFirst = bcIndex - 1;
         }
         if (penetrationFlagsFromOther[thisNodeIdFirstLocal]) {
            _boundaryNodeNormalStress[thisNodeIdFirstLocal] = boundaryNormalIncrementUpdated[thisNodeIdFirstLocal];
            _boundaryNodeShearStress[thisNodeIdFirstLocal] = boundaryShearIncrementUpdated[thisNodeIdFirstLocal];
            _boundaryNodeContactId[thisNodeIdFirstLocal] = other.getGrainId();
         }
         else if (connectionPenetrationFlag[bcIndexFirst] || connectionPenetrationFlag[bcIndexSecond]) {
            _boundaryNodeNormalStress[thisNodeIdFirstLocal].fill(0.);
            _boundaryNodeShearStress[thisNodeIdFirstLocal].fill(0.);
            _boundaryNodeContactId[thisNodeIdFirstLocal] = other.getGrainId();
         }
         else if (_boundaryNodeContactId[thisNodeIdFirstLocal] == other.getGrainId()) {
            _boundaryNodeContactId[thisNodeIdFirstLocal] = 50000;
            _boundaryNodeNormalStress[thisNodeIdFirstLocal].fill(0.);
            _boundaryNodeShearStress[thisNodeIdFirstLocal].fill(0.);
         }
       }
    }
    vector<VectorXd> globalContactForces(2);
    globalContactForces[0] = contactForceThis;
    globalContactForces[1] = contactForceOther;

    return globalContactForces;
  }

  double
  calculateContactRadiusRigidStraightBoundary(const RigidStraightWall & wall) const {
    double results = 0;

    double kn = _youngsModulus/0.1;

    vector<Vector2d> intersectionPoints;
    Vector2d oneIntersectionPoint;oneIntersectionPoint.fill(0.);
    double distance1; double distance2;
    vector<size_t> elementNodeIds;
    double integrationLength;
    Vector2d edgeNormal; edgeNormal.fill(0.);
    Vector2d nodePosition1, nodePosition2;

    size_t nodeId1, nodeId2;

    Vector2d normal = wall.getWallNormal();
    vector<Vector2d> boundaryPoints = wall.getBoundaryPoints();




    for (unsigned int bindex = 0; bindex < _boundaryNodeConnectivity.size(); bindex++) {
        
        nodeId1 = _boundaryNodeIds[_boundaryNodeConnectivity[bindex][0]];
        nodePosition1 = _nodePositions[nodeId1];
        distance1 = (nodePosition1-boundaryPoints[0]).dot(normal);
        if (distance1 < 0) {
          intersectionPoints.push_back(nodePosition1);
        }

        for (unsigned int cindex = 0; cindex < _pnumber; cindex++) {
          nodePosition1 = _subConnectionNodePositions[bindex][cindex];
          distance1 = (nodePosition1-boundaryPoints[0]).dot(normal);
          if (distance1 < 0) {
            intersectionPoints.push_back(nodePosition1);
          }
        }
    }

    double largestDistance = 0;
    double currentDistance = 0;
    for (unsigned int ii = 0; ii < intersectionPoints.size(); ii++) {
      for (unsigned int jj = ii+1; jj < intersectionPoints.size(); jj++) {
        currentDistance = (intersectionPoints[ii]-intersectionPoints[jj]).norm();
        if (currentDistance > largestDistance) {
          largestDistance = currentDistance;
        }
      }
    }
    results = largestDistance/2.;
    return results;
  }

  double
  computeContactRadiusWithOtherGrain (const StaticGrain & other) {
    vector<Vector2d> nodesInContact;
    for (unsigned int nindex = 0; nindex < _boundaryNodeIds.size(); nindex++) {
       if (_boundaryNodeContactId[nindex] == other.getGrainId()) {
          nodesInContact.push_back(_nodePositions[_boundaryNodeIds[nindex]]);
       }
       for (unsigned int sindex = 0; sindex < _pnumber; sindex++) {
          if (_subConnectionContactId[nindex][sindex] == other.getGrainId()) {
            nodesInContact.push_back(_subConnectionNodePositions[nindex][sindex]);
          }
       }
    }
    double largestDistance = 0.;
    double distanceTemp = 0;
    for (unsigned int ii = 0; ii < nodesInContact.size(); ii++) {
      for (unsigned int jj = ii+1; jj < nodesInContact.size(); jj++) {
          distanceTemp = (nodesInContact[ii]-nodesInContact[jj]).norm();
          if (largestDistance < distanceTemp) {
              largestDistance = distanceTemp;
          }
      }
    }
    return largestDistance/2.;
  }
  

  vector<Triplet<double,size_t> >
  assembleStiffnessMatrix(const vector<EssentialBC> & essentialBoundaryConditions) const {
    size_t numberDofs = _numberOfNodes*2;

    MatrixXd stiffnessMatrix(numberDofs, numberDofs);
    stiffnessMatrix.fill(0);
    vector<Triplet<double,size_t> > tripletList;
    tripletList.reserve(_numberOfNodes);

    vector<size_t> elementNodeIds;
    Matrix<double, 6, 6> elementStiffnessMatrix; elementStiffnessMatrix.fill(0.);
    size_t nodeId1;
    size_t nodeId2;
    for (size_t eindex = 0; eindex < _elements.size(); eindex++) {
      elementStiffnessMatrix = _elements[eindex].computeStiffnessMatrix();
      elementNodeIds = _elements[eindex].getNodeIds();
      for (size_t nodeindex1 = 0; nodeindex1 < elementNodeIds.size(); nodeindex1++) {
        nodeId1 = elementNodeIds[nodeindex1];
        for (size_t nodeindex2 = 0; nodeindex2 < elementNodeIds.size(); nodeindex2++) {
          nodeId2 = elementNodeIds[nodeindex2];
          for (size_t i = 0; i < 2; i++) {
            for (size_t j = 0; j < 2; j++) {
              stiffnessMatrix(nodeId1*2+i,nodeId2*2+j) += elementStiffnessMatrix(nodeindex1*2+i,nodeindex2*2+j);
            }
          }
        }
      }
    }
    for (unsigned int bcIndex = 0; bcIndex < essentialBoundaryConditions.size(); bcIndex++) {
      if (essentialBoundaryConditions[bcIndex].getGrainId()==_grainId) {
        nodeId1 = essentialBoundaryConditions[bcIndex].getNodeIdGrainBased()*2;
        nodeId2 = essentialBoundaryConditions[bcIndex].getDirection();
        stiffnessMatrix.row(nodeId1+nodeId2).fill(0);
        stiffnessMatrix(nodeId1+nodeId2,nodeId1+nodeId2) = 1.;
      }
    }
    for (size_t rowIndex = 0; rowIndex < numberDofs; rowIndex++) {
      for (size_t colIndex = 0; colIndex < numberDofs; colIndex++) {
        if (stiffnessMatrix(rowIndex,colIndex)!=0) {
          tripletList.push_back(Triplet<double,size_t>(rowIndex,colIndex,stiffnessMatrix(rowIndex,colIndex)));
        }
      }
    }
    return tripletList;
  }


  vector<Vector3d>
  computeElementStresses(const vector<Vector2d> & displacements) const {
    vector<Vector3d> allElementAvgStresses; allElementAvgStresses.resize(_elements.size());
    vector<size_t> elementNodeIds;
    vector<Vector2d> elementdisplacements;
    vector<Vector3d> elementStresses;
    Vector3d avgStress; avgStress.fill(0.);    
    for (unsigned int elementIndex = 0; elementIndex < _elements.size(); elementIndex++) {
      elementNodeIds = _elements[elementIndex].getNodeIds();
      elementdisplacements = Utilities::getElementDisplacementsFromGlobalList(elementNodeIds, displacements);
      elementStresses = _elements[elementIndex].computeStressesAtGaussianPoints(elementdisplacements);
      avgStress.fill(0.);
      for (unsigned int sindex = 0; sindex < elementStresses.size(); sindex++) {
        avgStress += elementStresses[sindex];
      }
      avgStress /= double(elementStresses.size());
      allElementAvgStresses[elementIndex] = avgStress;
    }
    return allElementAvgStresses;    
  }


  vector<Vector3d>
  computeElementStrains(const vector<Vector2d> & displacements) const {
    vector<Vector3d> allElementAvgStrains(_elements.size());
    vector<size_t> elementNodeIds;
    vector<Vector2d> elementdisplacements;
    vector<Vector3d> elementStrains;
    Vector3d avgStrain; avgStrain.fill(0.);
    for (unsigned int elementIndex = 0; elementIndex < _elements.size(); elementIndex++) {
      elementNodeIds = _elements[elementIndex].getNodeIds();
      elementdisplacements = Utilities::getElementDisplacementsFromGlobalList(elementNodeIds, displacements);
      elementStrains = _elements[elementIndex].computeStrainsAtGuassianPoints(elementdisplacements);
      avgStrain.fill(0.);
      for (unsigned int sindex = 0; sindex < elementStrains.size(); sindex++) {
        avgStrain += elementStrains[sindex];
      }
      avgStrain /= double(elementStrains.size());
      allElementAvgStrains[elementIndex] = avgStrain;
    }
    return allElementAvgStrains;
  }

  vector<Vector3d>
  computeNodalStresses(const vector<Vector2d> & displacements) const {
    vector<Vector3d> nodalStresses(displacements.size());
    vector<double> volumeSums(displacements.size());

    vector<size_t> elementNodeIds;
    vector<double> elementWeights;
    vector<Vector3d> elementStresses = computeElementStresses(displacements);

    for (unsigned int index = 0; index < displacements.size(); index++) {
      nodalStresses[index].fill(0.);
      volumeSums[index] = 0.;
    }

    for (unsigned int eindex = 0; eindex < _elements.size(); eindex++) {
      elementNodeIds = _elements[eindex].getNodeIds();
      elementWeights = _elements[eindex].getNodalWeights();
      for (unsigned int nindex = 0; nindex < elementNodeIds.size(); nindex++) {
        nodalStresses[elementNodeIds[nindex]] += elementStresses[eindex]*elementWeights[nindex];
        volumeSums[elementNodeIds[nindex]] += elementWeights[nindex];
      }
    }

    for (unsigned int nindex = 0; nindex < displacements.size(); nindex++) {
      nodalStresses[nindex] /= volumeSums[nindex];
    }
    return nodalStresses;
  }

  vector<Vector3d>
  computeNodalStrains(const vector<Vector2d> & displacements) const {
    vector<Vector3d> nodalStrains(displacements.size());
    vector<double> volumeSums(displacements.size());

    vector<size_t> elementNodeIds;
    vector<double> elementWeights;
    vector<Vector3d> elementStrains = computeElementStrains(displacements);

    for (unsigned int index = 0; index < displacements.size(); index++) {
      nodalStrains[index].fill(0.);
      volumeSums[index] = 0.;
    }

    for (unsigned int eindex = 0; eindex < _elements.size(); eindex++) {
      elementNodeIds = _elements[eindex].getNodeIds();
      elementWeights = _elements[eindex].getNodalWeights();
      for (unsigned int nindex = 0; nindex < elementNodeIds.size(); nindex++) {
        nodalStrains[elementNodeIds[nindex]] += elementStrains[eindex]*elementWeights[nindex];
        volumeSums[elementNodeIds[nindex]] += elementWeights[nindex];
      }
    }

    for (unsigned int nindex = 0; nindex < displacements.size(); nindex++) {
      nodalStrains[nindex] /= volumeSums[nindex];
    }
    return nodalStrains;
  }



  void
  updateNodePositions (const vector<Vector2d> & displacements) {
    vector<size_t> elementNodeIds;
    vector<Vector2d> nodalPositions;
    Vector2d centricPoint; centricPoint.fill(0.);
    double singleElementArea = 0.;
    Vector2d weightedPoint; weightedPoint.fill(0.);
    double totalArea = 0;
    for (unsigned int eindex = 0; eindex < _elements.size(); eindex++) {
      elementNodeIds = _elements[eindex].getNodeIds();
      nodalPositions.resize(elementNodeIds.size());
      centricPoint.fill(0.);
      for (unsigned int nindex = 0; nindex < elementNodeIds.size(); nindex++) {
        nodalPositions[nindex] = _nodePositions[elementNodeIds[nindex]]+displacements[elementNodeIds[nindex]];
        centricPoint += nodalPositions[nindex]/double(elementNodeIds.size());
      }
      _elements[eindex].updateNodalPositions(nodalPositions);
      singleElementArea = (nodalPositions[1](0)-nodalPositions[0](0))*(nodalPositions[2](1)-nodalPositions[0](1))/2.
        - (nodalPositions[2](0)-nodalPositions[0](0))*(nodalPositions[1](1)-nodalPositions[0](1))/2.;
      totalArea += singleElementArea;
      weightedPoint += centricPoint*singleElementArea;  
    }

    _centroid = weightedPoint/totalArea;
    double largestDistance = 0.;
    Vector2d singleNodePosition; singleNodePosition.fill(0.);
    for (unsigned int nindex = 0; nindex < _numberOfNodes; nindex++) {
      singleNodePosition = _nodePositions[nindex]+displacements[nindex];
      largestDistance = max(largestDistance, (singleNodePosition-_centroid).norm());     
    }
    _largestRadius = largestDistance;

    size_t firstNodeIdLocal;
    size_t secondNodeIdLocal;
    size_t firstNodeIdGlobal;
    size_t secondNodeIdGlobal;
    double ratio;
    for (unsigned int bIndex = 0; bIndex < _boundaryNodeIds.size(); bIndex++) {
      firstNodeIdLocal = _boundaryNodeConnectivity[bIndex][0];
      secondNodeIdLocal = _boundaryNodeConnectivity[bIndex][1];
      firstNodeIdGlobal = _boundaryNodeIds[firstNodeIdLocal];
      secondNodeIdGlobal = _boundaryNodeIds[secondNodeIdLocal];
      for (unsigned int sIndex = 0; sIndex < _pnumber; sIndex++) {
        ratio = (_subConnectionNodePositions[bIndex][sIndex]-_nodePositions[firstNodeIdGlobal]).norm()/
                  (_nodePositions[secondNodeIdGlobal]-_nodePositions[firstNodeIdGlobal]).norm();
        _subConnectionNodePositions[bIndex][sIndex] += displacements[firstNodeIdGlobal]*(1.-ratio)+ratio*displacements[secondNodeIdGlobal];
      }
    }
    for (unsigned int nindex = 0; nindex < _numberOfNodes; nindex++) {
      _nodePositions[nindex] += displacements[nindex];
    }
  }

  vector<Vector2d>
  getNodePositions() const {
    return _nodePositions;
  }

  double
  getYoungsModulus() const {
    return _youngsModulus;
  }

  size_t
  getGrainId() const {
    return _grainId;
  }

  size_t
  getPNumber() const {
    return _pnumber;
  }

  vector<size_t>
  getBoundaryNodeId() const {
    return _boundaryNodeIds;
  }

  vector< vector<size_t> >
  getBoundaryNodeConnectivity() const {
    return _boundaryNodeConnectivity;
  }


  size_t
  getNumberOfNodes() const {
    return _numberOfNodes;
  }
  
  Vector2d
  getCentroid() const {
    return _centroid;
  }

  double
  getLargestRadius() const {
    return _largestRadius;
  }

  double
  getFrictionCoefficient() const {
    return _frictionCoefficient;
  }

  void
  changeFrictionCoefficient(const double & mu) {
    _frictionCoefficient = mu;
  }

  size_t
  getNumberOfElements() const {
    return _elements.size();
  }

  vector<size_t>
  getBoundaryContactId() const {
    return _boundaryNodeContactId;
  }

  vector<Vector2d>
  getBoundaryNodeNormalStress() const {
    return _boundaryNodeNormalStress;
  }

  vector<Vector2d>
  getBoundaryNodeShearStress() const {
    return _boundaryNodeShearStress;
  }

  vector< vector<Vector2d> > 
  getSubConnectionNodeNormalStress() const {
    return _subConnectionNodeNormalStress;
  }

  vector< vector<Vector2d> >
  getSubConnectionNodeShearStress() const {
    return _subConnectionNodeShearStress;
  }

  vector< vector<size_t> >
  getSubConnectionNodeContactId() const {
    return _subConnectionContactId;
  }

  vector< vector<Vector2d> >
  getSubConnectionNodePosition() const {
    return _subConnectionNodePositions;
  }

  vector< vector<size_t> >
  getBoundaryNodeIdsInContact() const {
    vector<size_t> singleContact(2);
    vector< vector<size_t> > nodeIdsWithContact;
    for (unsigned int nindex = 0; nindex < _boundaryNodeIds.size(); nindex++) {
       if (_boundaryNodeContactId[nindex] != 50000) {
          singleContact[0] = _boundaryNodeIds[nindex];
          singleContact[1] = _boundaryNodeContactId[nindex];
          nodeIdsWithContact.push_back(singleContact);       
       }
    }
    return nodeIdsWithContact;
  }

  void
  eraseNodalNormalAndShear() {
    for (unsigned int nindex = 0; nindex < _boundaryNodeIds.size(); nindex++) {
      _boundaryNodeNormalStress[nindex].fill(0.);
      _boundaryNodeShearStress[nindex].fill(0.);
      _boundaryNodeContactId[nindex] = 50000;
    }
    for (unsigned int cindex = 0; cindex < _boundaryNodeConnectivity.size(); cindex++) {
      for (unsigned int subIndex = 0; subIndex < _pnumber; subIndex++) {
        _subConnectionNodeNormalStress[cindex][subIndex].fill(0.);
        _subConnectionNodeShearStress[cindex][subIndex].fill(0.);
        _subConnectionContactId[cindex][subIndex] = 50000;
      }
    }
  }



private:
  size_t                      _grainId;
  vector<SimplexTriangle>     _elements;
  size_t                      _numberOfNodes;
  vector<Vector2d>            _nodePositions;                              // all node positions
  vector<size_t>              _boundaryNodeIds;
  vector<Vector2d>            _boundaryNodeShearStress;                    // this got updated after each newton rahpson iteration converged          
  vector<Vector2d>            _boundaryNodeNormalStress;
  vector<size_t>              _boundaryNodeContactId;                      // Id of the other object *this object is in contact with  
  vector< vector<size_t> >    _boundaryNodeConnectivity;                   // boundary edge mesh using boundary node local Id
  double                      _youngsModulus;
  double                      _frictionCoefficient;
                                       

  Vector2d                    _centroid;
  double                      _largestRadius;

  size_t                      _pnumber;                                    // number of points along each connectivity
  vector< vector<Vector2d> >  _subConnectionNodeShearStress;               // points along each connectivity store the shear stress (stores history)
  vector< vector<Vector2d> >  _subConnectionNodeNormalStress;
  vector< vector<size_t> >    _subConnectionContactId;                     // points along each connectivity store the id of the other grain in contact
  vector< vector<Vector2d> >  _subConnectionNodePositions;                 // positions of these points
  

};
#endif
