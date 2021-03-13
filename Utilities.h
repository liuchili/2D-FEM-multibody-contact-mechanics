// -*- C++ -*-
#ifndef UTILITIES_H
#define UTILITIES_H

#include "Definitions.h"

namespace Utilities {

	vector<Vector2d>
	getElementDisplacementsFromGlobalList(const vector<size_t> & elementNodeIds,
		                                  const vector<Vector2d> & displacements) {
		vector<Vector2d> elementDisplacements(elementNodeIds.size());
		for (unsigned int index = 0; index < elementNodeIds.size(); index++) {
			elementDisplacements[index] = displacements[elementNodeIds[index]];
		}
	return elementDisplacements;
	}

	vector<Vector2d>
	distributeGlobalVectorsToLocalVectors(const VectorXd & globalVector) {
		const unsigned int Dof = 2;
		unsigned int numberOfNodes = globalVector.size()/Dof;
		vector<Vector2d> localVectors(numberOfNodes);
		for (unsigned int nindex = 0; nindex < numberOfNodes; nindex++) {
			for (unsigned int coordinate = 0; coordinate < Dof; coordinate++) {
				localVectors[nindex](coordinate) = globalVector(nindex*Dof+coordinate);
			}
		}
		return localVectors; 
	}

	VectorXd
	convertLocalNodalDisplacementToGlobalList(const vector<Vector2d> & localDisplacements) {
		VectorXd globalDisplacements; globalDisplacements.resize(localDisplacements.size()*2);
		for (unsigned int nindex = 0; nindex < localDisplacements.size(); nindex++) {
			globalDisplacements[nindex*2] = localDisplacements[nindex](0);
			globalDisplacements[nindex*2+1] = localDisplacements[nindex](1);
		}
		return globalDisplacements;
	}

	vector< vector<Vector2d> > 
	convertWorldDisplacementVectorToGrainForm(const VectorXd & worldDisplacement, const vector<size_t> & blockSizes) {
		vector< vector<Vector2d> > worldDisplacementsGrainForm(blockSizes.size());
		vector<Vector2d> singleGrainDisplacements;
		VectorXd singleGrainDisplacementsGlobalForm;
		size_t startingPoint = 0;
		for (unsigned int gindex = 0; gindex < blockSizes.size(); gindex++) {
			//singleGrainDisplacementsGlobalForm = worldDisplacement.block<blockSizes[gindex],1>(startingPoint,0);
			singleGrainDisplacementsGlobalForm = worldDisplacement.block(startingPoint, 0, blockSizes[gindex]*2, 1);
			singleGrainDisplacements = distributeGlobalVectorsToLocalVectors(singleGrainDisplacementsGlobalForm);
			worldDisplacementsGrainForm[gindex] = singleGrainDisplacements;
			startingPoint += blockSizes[gindex]*2;
		}
		return worldDisplacementsGrainForm;
	}

	VectorXd
	convertWorldDisplacementGrainFormToVectorForm(const vector< vector<Vector2d> > & worldDisplacement, const vector<size_t> & blockSizes) {
		size_t nDof = 0;
		for (unsigned int gindex = 0; gindex < blockSizes.size(); gindex++) {
			nDof += blockSizes[gindex]*2;
		}
		VectorXd worldDisplacementVector(nDof);
		worldDisplacementVector.fill(0.);
		size_t startingPoint = 0;
		for (unsigned int gindex = 0; gindex < blockSizes.size(); gindex++) {
			worldDisplacementVector.block(startingPoint,0,blockSizes[gindex]*2,1) 
			  = convertLocalNodalDisplacementToGlobalList(worldDisplacement[gindex]);
			startingPoint += blockSizes[gindex]*2;
		}
		return worldDisplacementVector;
	}
	
	Matrix<double, 2, 2>
	convertVoigtStressToMatrixForm (const Vector3d & voigtStress) {
		Matrix<double, 2, 2> stress; stress.fill(0.);
		stress(0,0) = voigtStress(0);
		stress(0,1) = voigtStress(2);
		stress(1,0) = voigtStress(2);
		stress(1,1) = voigtStress(1);
		return stress;
	}

	
    bool
    computePointToSegments(const Vector2d & point, const Vector2d & pointDisplacement, 
        const vector<Vector2d> & otherNodePositions, const vector<Vector2d> & otherNodeDisplacements, const double & thresholdDistance, const double & thresPortion,
        const Vector2d & pointDisplacementAccumulate, const vector<Vector2d> & otherNodeDisplacementsAccumulate,
        size_t & projectionFlag, size_t & projectionId, double & ratio, bool & penetrationFlag, size_t & connectionId,
        double & normalDistance, Vector2d & normal, Vector2d & tangentialDistance, size_t & projectionIdAlternative,
        double & ratioAlternative, Vector2d & normalDistanceVectorAlternative, Vector2d & tangentialDistanceVectorAlternative) {


    	bool checkFlag = false;

        projectionFlag = 0;
        penetrationFlag = false;
        projectionId = 0;
        connectionId = 0;
        ratio = 0.;
        normalDistance = 10000.;
        tangentialDistance.fill(0.);
        normal.fill(0.);


        projectionIdAlternative = 0;
        ratioAlternative = 0.;
        normalDistanceVectorAlternative.fill(0.);
        tangentialDistanceVectorAlternative.fill(0.);

        Vector2d projectedPoint; projectedPoint.fill(0.);

        vector<Vector2d> distanceVectors(3);
        vector<double>   distances(3);
        for (unsigned int index = 0; index < 3; index++) {
            distanceVectors[index] = point+pointDisplacement-otherNodePositions[index]-otherNodeDisplacements[index];
            distances[index] = distanceVectors[index].norm();
        }
        double edgeLength0, edgeLength1;
        double normalDistance0, normalDistance1;    
        Vector2d tangent0 = otherNodePositions[1]+otherNodeDisplacements[1]-otherNodePositions[0]-otherNodeDisplacements[0];
        edgeLength0 = tangent0.norm();
        tangent0 /= tangent0.norm();
        Vector2d normal0; normal0 << -tangent0(1), tangent0(0);

        Vector2d tangent1 = otherNodePositions[2]+otherNodeDisplacements[2]-otherNodePositions[1]-otherNodeDisplacements[1];
        edgeLength1 = tangent1.norm();
        tangent1 /= tangent1.norm();
        Vector2d normal1; normal1 << -tangent1(1), tangent1(0);

        bool projectableFirst, projectableSecond;
        projectableFirst = false;
        projectableSecond = false;


        
        if (!(distances[0]>thresholdDistance && distances[1]>thresholdDistance && distances[2]>thresholdDistance)
              && distanceVectors[0].dot(tangent0)*distanceVectors[2].dot(tangent1) <= 0) { 
            checkFlag = true;
            normalDistance0 = distanceVectors[0].dot(normal0);
            normalDistance1 = distanceVectors[1].dot(normal1);
            if (distances[0]*distances[1]*distances[2]==0) { 
                projectionFlag = 2;
                ratio = 1.;
                normalDistance = 0.;
                penetrationFlag = false;
                normal.fill(0.);
                tangentialDistance.fill(0.);
                if (distances[0]==0) {
                    projectionId = 0;
                }
                else if (distances[1]==0) {
                    projectionId = 1;
                }
                else {
                    projectionId = 2;
                }
            }
            else if (fabs(normalDistance0)==0 && distances[0]+distances[1] <= edgeLength0) {
                projectionFlag = 1;
                penetrationFlag = false;
                ratio = distances[0]/edgeLength0;
                projectionId = 0;
                normalDistance = 0.;
                normal.fill(0.);
                tangentialDistance.fill(0.);               
                if (distanceVectors[1].dot(tangent1)*distanceVectors[2].dot(tangent1) <= 0 && normalDistance1 > 0) {
                	penetrationFlag = true;
                	projectionIdAlternative = 1;
                	ratioAlternative = distanceVectors[1].dot(tangent1)/edgeLength1;
                	normalDistanceVectorAlternative = normalDistance1*(-1.)*normal1;
                	tangentialDistanceVectorAlternative << normal1(1), -normal1(0);
                }
                else if (distanceVectors[1].dot(tangent1)*distanceVectors[2].dot(tangent1) > 0 && normalDistance1 > 0) {
                	if (fabs(distanceVectors[1].dot(tangent1))<thresPortion*edgeLength1 && tangent0.dot(tangent1) < sqrt(3)/2) {	
                		penetrationFlag = true;
                		projectionIdAlternative = 1;
                		ratioAlternative = 0;
                		normalDistanceVectorAlternative = normalDistance1*(-1.)*normal1;
                		tangentialDistanceVectorAlternative << normal1(1), -normal1(0);
                	}
                }                
            }
            else if (fabs(normalDistance1)==0 && distances[1]+distances[2] <= edgeLength1) {
                projectionFlag = 1;
                penetrationFlag = false;
                ratio = distances[1]/edgeLength1;
                projectionId = 1;
                normalDistance = 0.;
                normal.fill(0.);
                tangentialDistance.fill(0.);                
                if (distanceVectors[0].dot(tangent0)*distanceVectors[1].dot(tangent0) <= 0 && normalDistance0 > 0) {
                	penetrationFlag = true;
                	projectionIdAlternative = 0;
                	ratioAlternative = distanceVectors[0].dot(tangent0)/edgeLength0;
                	normalDistanceVectorAlternative = normalDistance0*(-1.)*normal0;
                	tangentialDistanceVectorAlternative << normal0(1), -normal0(0);               	
                }
                else if (distanceVectors[0].dot(tangent0)*distanceVectors[1].dot(tangent0) > 0 && normalDistance0 > 0) {
                	if (fabs(distanceVectors[1].dot(tangent0))<thresPortion*edgeLength0 && tangent0.dot(tangent1) < sqrt(3)/2) {
                		penetrationFlag = true;
                		projectionIdAlternative = 0;
                		ratioAlternative = 1;
                		normalDistanceVectorAlternative = normalDistance0*(-1.)*normal0;
                		tangentialDistanceVectorAlternative << normal0(1), -normal0(0);
                	}
                }
                                
            }
            else {

                if (distanceVectors[0].dot(tangent0)*distanceVectors[1].dot(tangent0) <= 0) {
                    projectableFirst = true;
                }
                else {
                    projectableFirst = false;
                }
                if (distanceVectors[1].dot(tangent1)*distanceVectors[2].dot(tangent1) <= 0) {
                    projectableSecond = true;
                }
                else {
                    projectableSecond = false;
                }

                if (projectableFirst || projectableSecond) {
                    if (projectableFirst && (!projectableSecond)) {
                        projectionFlag = 1;
                        projectionId = 0;
                        ratio = distanceVectors[0].dot(tangent0)/edgeLength0;
                        penetrationFlag = false;
                        normalDistance = -1.*normalDistance0;
                        if (normalDistance < 0) {                          
                            normal = -1.*normal0;
                            tangentialDistance << -normal(1), normal(0);
                            penetrationFlag = true;
                            if (normalDistance1 > 0) {
                            	if (fabs(distanceVectors[1].dot(tangent1))<thresPortion*edgeLength1 && tangent0.dot(tangent1) < sqrt(3)/2) {		
                            		projectionIdAlternative = 1;
                            		ratioAlternative = 0;
                            		normalDistanceVectorAlternative = normalDistance1*(-1.)*normal1;
                            		tangentialDistanceVectorAlternative << normal1(1), -normal1(0);
                            	}
                            }
                        }
                    }
                    else if ((!projectableFirst) && projectableSecond) {
                        projectionFlag = 1;
                        projectionId = 1;
                        ratio = distanceVectors[1].dot(tangent1)/edgeLength1;
                        penetrationFlag = false;
                        normalDistance = -1.*normalDistance1;
                        if (normalDistance < 0) {
                            penetrationFlag = true;
                            normal = -1.*normal1;
                            tangentialDistance << -normal(1), normal(0);
                            if (normalDistance0 > 0) {
                            	if (fabs(distanceVectors[1].dot(tangent0))<thresPortion*edgeLength0 && tangent0.dot(tangent1) < sqrt(3)/2) {	
                            		projectionIdAlternative = 0;
                            		ratioAlternative = 1;
                            		normalDistanceVectorAlternative = normalDistance0*(-1.)*normal0;
                            		tangentialDistanceVectorAlternative << normal0(1), -normal0(0);
                            	}
                            }
                        }
                    }
                    else if (projectableFirst && projectableSecond) {
                        if (normalDistance0 > 0 && normalDistance1 > 0) {
                            if (normalDistance0 > normalDistance1) {
                                projectionFlag = 1;
                                projectionId = 1;
                                penetrationFlag = true;
                                ratio = distanceVectors[1].dot(tangent1)/edgeLength1;
                                normalDistance = -1.*normalDistance1;
                                normal = -1.*normal1;
                                tangentialDistance << -normal(1), normal(0);

                                projectionIdAlternative = 0;
                                ratioAlternative = distanceVectors[0].dot(tangent0)/edgeLength0;
                	            normalDistanceVectorAlternative = normalDistance0*(-1.)*normal0;
                	            tangentialDistanceVectorAlternative << normal0(1),-normal0(0);
                            }
                            else {
                                projectionFlag = 1;
                                projectionId = 0;
                                penetrationFlag = true;
                                ratio = distanceVectors[0].dot(tangent0)/edgeLength0;
                                normalDistance = -1.*normalDistance0;
                                normal = -1.*normal0;
                                tangentialDistance << -normal(1), normal(0);
                                projectionIdAlternative = 1;
                                ratioAlternative = distanceVectors[1].dot(tangent1)/edgeLength1;
                	            normalDistanceVectorAlternative = normalDistance1*(-1.)*normal1;
                	            tangentialDistanceVectorAlternative << normal1(1), -normal1(0); 
                            }
                        }
                        else if (normalDistance0 > 0 && normalDistance1 < 0) {
                            projectionFlag = 1;
                            projectionId = 0;
                            ratio = distanceVectors[0].dot(tangent0)/edgeLength0;
                            penetrationFlag = true;
                            normalDistance = -1.*normalDistance0;
                            normal = -1.*normal0;
                            tangentialDistance << -normal(1), normal(0);
                        }
                        else if (normalDistance0 < 0 && normalDistance1 > 0) {
                            projectionFlag = 1;
                            projectionId = 1;
                            ratio = distanceVectors[1].dot(tangent1)/edgeLength1;
                            penetrationFlag = true;
                            normalDistance = -1.*normalDistance1;
                            normal = -1.*normal1;
                            tangentialDistance << -normal(1), normal(0);
                        } 
                        else {
                            if (normalDistance0 > normalDistance1) {
                                projectionFlag = 1;
                                projectionId = 0;
                                ratio = distanceVectors[0].dot(tangent0)/edgeLength0;
                                penetrationFlag = false;
                                normalDistance = -1.*normalDistance0;
                            }
                            else {
                                projectionFlag = 1;
                                projectionId = 1;
                                ratio = distanceVectors[1].dot(tangent1)/edgeLength1;
                                penetrationFlag = false;
                                normalDistance = -1.*normalDistance1;
                            }                        
                        }
                    }
                }
                else {

                    double smallestDistance = min(min(distances[0],distances[1]),distances[2]);
                    projectionFlag = 2;
                    ratio = 1.;
                    if (distanceVectors[0].norm() == smallestDistance) {
                        projectionId = 0;
                        normal = -1.*normal0;
                    }
                    else if (distanceVectors[1].norm() == smallestDistance) {
                        projectionId = 1;
                        normal = -1.*(edgeLength0*normal0+edgeLength1*normal1)/(edgeLength0+edgeLength1);
                        normal /= normal.norm();                        
                    }
                    else {
                       projectionId = 2;
                       normal = -1.*normal1;
                    }
                    if (normalDistance1 > 0 && normalDistance0 > 0) {
                       penetrationFlag = true;
                       normalDistance = -1.*smallestDistance;
                       tangentialDistance << -normal(1), normal(0);
                    }
                    else if (normalDistance1 < 0 && normalDistance0 < 0) {
                       penetrationFlag = false;
                       normalDistance = smallestDistance;
                       tangentialDistance.fill(0.);
                       normal.fill(0.);
                    }
                }
            }
	        if (projectionFlag == 2) {
	            if (projectionId == 0) {
	                connectionId = 0;
	                ratio = 0.;
	            }
	            else if (projectionId == 1){
	                connectionId = 0;
	                ratio = 1.;
	            }
	            else if (projectionId ==2) {
	            	connectionId = 1;
	            	ratio = 1.;
	            }
	        }
            
	        if (penetrationFlag) {
	        	double relativeDisplacement;
	        	double relativeDisplacementAlternative = 0;
	        	if (projectionFlag == 2) {
	        		relativeDisplacement = (pointDisplacement+pointDisplacementAccumulate-otherNodeDisplacements[projectionId]-
	        			otherNodeDisplacementsAccumulate[projectionId]).dot(tangentialDistance);
	        	}
	        	else if (projectionFlag == 1) {
	        		relativeDisplacement = (pointDisplacement+pointDisplacementAccumulate-((1-ratio)*(otherNodeDisplacements[projectionId]+otherNodeDisplacementsAccumulate[projectionId])
	        			+ratio*(otherNodeDisplacements[projectionId+1]+otherNodeDisplacementsAccumulate[projectionId+1]))).dot(tangentialDistance);

	        		relativeDisplacementAlternative = (pointDisplacement+pointDisplacementAccumulate-((1-ratioAlternative)*(otherNodeDisplacements[projectionIdAlternative]+otherNodeDisplacementsAccumulate[projectionIdAlternative])+
	        			ratioAlternative*(otherNodeDisplacements[projectionIdAlternative+1]+otherNodeDisplacementsAccumulate[projectionIdAlternative+1]))).dot(tangentialDistanceVectorAlternative);
	        	}
	            if (relativeDisplacement >= 0) {tangentialDistance = relativeDisplacement*(-1.*tangentialDistance);}
	            else {tangentialDistance = fabs(relativeDisplacement)*tangentialDistance;}

	            if (relativeDisplacementAlternative >= 0) {tangentialDistanceVectorAlternative = relativeDisplacementAlternative*(-1.*tangentialDistanceVectorAlternative);}
	            else {tangentialDistanceVectorAlternative = fabs(relativeDisplacementAlternative)*tangentialDistanceVectorAlternative;} 	              
	        }
	        
	        	        
	    }
	    return checkFlag;
    }

    void
    computePointToSingleSegment(const Vector2d & point, const Vector2d & pointDisplacement, const Vector2d & pointDisplacementAccumulate, const Vector2d & nodePosition1, const Vector2d & nodePosition2,
    	const Vector2d & nodeDisplacement1, const Vector2d nodeDisplacement2, const Vector2d & nodeDisplacementAccumulate1, const Vector2d & nodeDisplacementAccumulate2, bool & checkFlag,
    	Vector2d & normalDistance, Vector2d & tangentialDirection, Vector2d & relativeDisplacement, size_t & outsideFlagDeformed, double & ratio) {

    	tangentialDirection.fill(0.);
    	Vector2d tangent = nodePosition2+nodeDisplacement2-nodePosition1-nodeDisplacement1;
    	double edgeLength = tangent.norm();
    	tangent /= tangent.norm();
    	tangentialDirection = tangent;
    	Vector2d normal; normal << -tangent(1), tangent(0);


    	normalDistance.fill(0.);
    	relativeDisplacement.fill(0.);

    	Vector2d distanceVector1 = point+pointDisplacement-nodePosition1-nodeDisplacement1;
    	Vector2d distanceVector2 = point+pointDisplacement-nodePosition2-nodeDisplacement2;

        checkFlag = false;
    	outsideFlagDeformed = 1000;
    	ratio = 0.;

    	if (distanceVector1.dot(tangent)*distanceVector2.dot(tangent) <= 0) {
    		outsideFlagDeformed = 1;
    		ratio = distanceVector1.dot(tangent)/edgeLength;
    	}
    	else if (distanceVector1.dot(tangent) < 0 && distanceVector2.dot(tangent) < 0) {
    		outsideFlagDeformed = 0;
    		ratio = 0;
    	}
    	else if (distanceVector1.dot(tangent) > 0 && distanceVector2.dot(tangent) < 0) {
    		outsideFlagDeformed = 2;
    		ratio = 1;
    	}
    	if (distanceVector1.dot(normal) >= 0) {
    		checkFlag = true;
    	}

		normalDistance = -1.*normal*distanceVector1.dot(normal);
		relativeDisplacement = pointDisplacement+pointDisplacementAccumulate-
		  ((1.-ratio)*(nodeDisplacement1+nodeDisplacementAccumulate1)+ratio*(nodeDisplacement2+nodeDisplacementAccumulate2))
					  
    }
    
    bool
    identifyOnePointProjectionToLineSegments(const double & thresDistance, const Vector2d & point, const Vector2d & pointDisplacement, 
                                             const vector< vector<size_t> > & boundaryConnectionsOfOtherToCheck, const vector<size_t> & boundaryConnectionIndex,
                                             const vector<Vector2d> & allBoundaryNodePositionsOfOther, const vector<Vector2d> & allBoundaryNodeDisplacementsOfOther,
                                             vector<size_t> & localSearchNodeIds, vector<size_t> & localSearchConnections) {

        // first do in deformed configuration
        // Loop through all other node in need of checking
        double smallestDistance = 10000.;
        double distanceTemp = 0.;
        size_t otherNodeIdWithSmallestDistance;
        size_t connectionIndexLocal = 0;
        bool   considerable = false;

        Vector2d otherNodePosition;

        localSearchNodeIds.resize(0);
        localSearchConnections.resize(0);

        if (boundaryConnectionsOfOtherToCheck[boundaryConnectionsOfOtherToCheck.size()-1][1]!=
            boundaryConnectionsOfOtherToCheck[0][0]) { // open connections
            for (unsigned int bcIndex = 0; bcIndex < boundaryConnectionsOfOtherToCheck.size()+1; bcIndex++) {
                if (bcIndex < boundaryConnectionsOfOtherToCheck.size()) {
                    otherNodePosition = allBoundaryNodePositionsOfOther[boundaryConnectionsOfOtherToCheck[bcIndex][0]]+
                      allBoundaryNodeDisplacementsOfOther[boundaryConnectionsOfOtherToCheck[bcIndex][0]];                      
                }
                else {
                    otherNodePosition = allBoundaryNodePositionsOfOther[boundaryConnectionsOfOtherToCheck[bcIndex-1][1]]+
                      allBoundaryNodeDisplacementsOfOther[boundaryConnectionsOfOtherToCheck[bcIndex-1][1]];                     
                }
                distanceTemp = (point+pointDisplacement-otherNodePosition).norm();
                if (smallestDistance > distanceTemp) {
                    smallestDistance = distanceTemp;
                    connectionIndexLocal = bcIndex;
                }
            }
        }
        else { // close connections
            for (unsigned int bcIndex = 0; bcIndex < boundaryConnectionsOfOtherToCheck.size(); bcIndex++) {
                otherNodePosition = allBoundaryNodePositionsOfOther[boundaryConnectionsOfOtherToCheck[bcIndex][0]]+
                    allBoundaryNodeDisplacementsOfOther[boundaryConnectionsOfOtherToCheck[bcIndex][0]];
                distanceTemp = (point+pointDisplacement-otherNodePosition).norm();
                if (smallestDistance > distanceTemp) {
                    smallestDistance = distanceTemp;
                    connectionIndexLocal = bcIndex;
                }                
            }
        }

        if (smallestDistance < thresDistance) {
            considerable = true;
            if (boundaryConnectionsOfOtherToCheck[boundaryConnectionsOfOtherToCheck.size()-1][1]!=
                boundaryConnectionsOfOtherToCheck[0][0]) {
                if (connectionIndexLocal == 0) {
                    localSearchNodeIds.resize(2);
                    localSearchNodeIds[0] = boundaryConnectionsOfOtherToCheck[0][0];
                    localSearchNodeIds[1] = boundaryConnectionsOfOtherToCheck[0][1];
                    localSearchConnections.resize(1);
                    localSearchConnections[0] = boundaryConnectionIndex[0];
                }
                else if (connectionIndexLocal < boundaryConnectionsOfOtherToCheck.size()) {
                    localSearchNodeIds.resize(3);
                    localSearchNodeIds[0] = boundaryConnectionsOfOtherToCheck[connectionIndexLocal-1][0];
                    localSearchNodeIds[1] = boundaryConnectionsOfOtherToCheck[connectionIndexLocal][0];
                    localSearchNodeIds[2] = boundaryConnectionsOfOtherToCheck[connectionIndexLocal][1];
                    localSearchConnections.resize(2);
                    localSearchConnections[0] = boundaryConnectionIndex[connectionIndexLocal-1];
                    localSearchConnections[1] = boundaryConnectionIndex[connectionIndexLocal];

                }
                else {
                    localSearchNodeIds.resize(2);
                    localSearchNodeIds[0] = boundaryConnectionsOfOtherToCheck[boundaryConnectionsOfOtherToCheck.size()-1][0];
                    localSearchNodeIds[1] = boundaryConnectionsOfOtherToCheck[boundaryConnectionsOfOtherToCheck.size()-1][1];
                    localSearchConnections.resize(1);
                    localSearchConnections[0] = boundaryConnectionIndex[boundaryConnectionsOfOtherToCheck.size()-1];
                }
            }
            else {
                localSearchNodeIds.resize(3);
                localSearchConnections.resize(2);
                localSearchNodeIds[1] = boundaryConnectionsOfOtherToCheck[connectionIndexLocal][0];
                localSearchNodeIds[2] = boundaryConnectionsOfOtherToCheck[connectionIndexLocal][1];
                localSearchConnections[1] = boundaryConnectionIndex[connectionIndexLocal];
                if (connectionIndexLocal == 0) {
                    localSearchNodeIds[0] = boundaryConnectionsOfOtherToCheck[boundaryConnectionsOfOtherToCheck.size()-1][0];
                    localSearchConnections[0] = boundaryConnectionIndex[boundaryConnectionsOfOtherToCheck.size()-1];
                }
                else {
                    localSearchNodeIds[0] = boundaryConnectionsOfOtherToCheck[connectionIndexLocal-1][0];
                    localSearchConnections[0] = boundaryConnectionIndex[connectionIndexLocal-1];
                }                
            }
        }
        return considerable;
    }    

	Vector2d
	RotateVectorToCurrentFrame (const Vector2d & normalPrevious, const Vector2d & normalCurrent, const Vector2d & quantity) {
		Vector2d rotatedQuantity; rotatedQuantity.fill(0.);

		double k = 0;
		double j = 0;

        k = normalPrevious(0)*normalCurrent(1)-normalPrevious(1)*normalCurrent(0);
        j = normalPrevious.dot(normalCurrent);
        double theta;
        if (k >=0 && j >= 0) {theta = asin(k);}
        if (k >=0 && j <= 0) {theta = acos(j);}
        if (k <=0 && j >= 0) {theta = 2*M_PI-fabs(asin(k));}
        if (k <=0 && j <= 0) {theta = 2*M_PI-fabs(acos(j));}

        double sint = sin(theta);
        double cost = cos(theta);

        rotatedQuantity(0) = cost*quantity(0)-sint*quantity(1);
        rotatedQuantity(1) = sint*quantity(0)+cost*quantity(1);       

        return rotatedQuantity;
	}
}

#endif  // UTILITIES_H
