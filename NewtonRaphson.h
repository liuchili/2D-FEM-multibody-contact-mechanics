// -*- C++ -*-
#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include "Utilities.h"



void
solveDisplacementFieldThroughNewtonRaphson(WorldStatic & world, 
	const size_t & niter, const vector< vector<Vector2d> > & worldDisplacementAccumulate,
	vector< vector<Vector2d> > & worldDisplacementIncrementInitialGuess,
	const vector< vector<Vector2d> > & worldWallDisplacements,
	const vector< vector<Vector2d> > & worldWallDisplacementAccumulate,
	const vector<EssentialBC>  & essentialBoundaryConditions,
	const VectorXd  & externalForceVector,
	const Vector2d & acceleration) {
	const size_t nDof = world.getTotalNodeNumber()*2;
	const vector<size_t> blockSizes = world.getBlockSizes();


	size_t globalLocation = 0;
	size_t grainId = 0;

    vector<Triplet<double,size_t> > stiffnessTriplet;	

    stiffnessTriplet =  world.assembleWorldStiffnessMatrix(essentialBoundaryConditions);
    SparseMatrix<double> worldStiffNessMatrix(nDof,nDof);
    worldStiffNessMatrix.setFromTriplets(stiffnessTriplet.begin(), stiffnessTriplet.end());
    worldStiffNessMatrix.makeCompressed();


	VectorXd worldDisplacementIncrementVector = Utilities::convertWorldDisplacementGrainFormToVectorForm(worldDisplacementIncrementInitialGuess, blockSizes);
	VectorXd worldDisplacementAccumulateVector = Utilities::convertWorldDisplacementGrainFormToVectorForm(worldDisplacementAccumulate, blockSizes);
	VectorXd worldDisplacementVectorTotalVector = worldDisplacementAccumulateVector+worldDisplacementIncrementVector;
	VectorXd worldBodyForceVector = world.assembleWorldBodyForceVector(acceleration);

    
	VectorXd worldForceVector(nDof);
	worldForceVector.fill(0.);

	

	SparseMatrix<double> worldSparseTangent(nDof,nDof);

	vector<Triplet<double,size_t> > wallTriplet;
	vector<Triplet<double,size_t> > grainTriplet;
	vector<Triplet<double,size_t> > totalTriplet;

	

	size_t updateFlag = 0;

	VectorXd wallForceVector(nDof);
	VectorXd internalForceVector(nDof);
	VectorXd grainForceVector(nDof);

	double penaltyratio = 0.1;


	for(unsigned int iter = 0; iter < niter; iter++) {

		if(iter == niter - 1) {
			updateFlag = 1;
		}

		internalForceVector = worldStiffNessMatrix*worldDisplacementVectorTotalVector;
		wallForceVector = world.assembleWorldContactForceVectorAgainstRigidStraightBoundaries(worldDisplacementIncrementInitialGuess,
			worldDisplacementAccumulate,worldWallDisplacements,worldWallDisplacementAccumulate,updateFlag,wallTriplet,penaltyratio);
		grainForceVector = world.assembleWorldContactForceBetweenGrains(worldDisplacementIncrementInitialGuess,worldDisplacementAccumulate,
		  	updateFlag,grainTriplet,penaltyratio);



       totalTriplet.clear();
       totalTriplet.insert(totalTriplet.end(), stiffnessTriplet.begin(), stiffnessTriplet.end());
       totalTriplet.insert(totalTriplet.end(), wallTriplet.begin(), wallTriplet.end());
       totalTriplet.insert(totalTriplet.end(), grainTriplet.begin(), grainTriplet.end());

		
		worldForceVector = internalForceVector-wallForceVector-grainForceVector-worldBodyForceVector-externalForceVector;

        worldSparseTangent.setFromTriplets(totalTriplet.begin(),totalTriplet.end());
        worldSparseTangent.makeCompressed();

		for(unsigned int bcIndex = 0; bcIndex < essentialBoundaryConditions.size(); bcIndex++) {
			grainId = essentialBoundaryConditions[bcIndex].getGrainId();
			globalLocation = accumulate(blockSizes.begin(),
				  blockSizes.begin()+grainId,0)*2+
			      essentialBoundaryConditions[bcIndex].getNodeIdGrainBased()*2+
			      essentialBoundaryConditions[bcIndex].getDirection();			
			worldForceVector(globalLocation) = 0.;
		}

		SparseLU<SparseMatrix<double> > solver;
		solver.analyzePattern(worldSparseTangent);
		solver.factorize(worldSparseTangent); 
		worldDisplacementVectorTotalVector -= solver.solve(worldForceVector);

		worldDisplacementIncrementVector = worldDisplacementVectorTotalVector-worldDisplacementAccumulateVector;
		worldDisplacementIncrementInitialGuess = Utilities::convertWorldDisplacementVectorToGrainForm(worldDisplacementIncrementVector,blockSizes);
		cout << "At " << iter+1 << " iteration the penaltyratio is " << penaltyratio << " with the residue " << worldForceVector.norm() << endl;
	}
}

#endif  // NEWTONRAPHSON_H
