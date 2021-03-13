// -*- C++ -*-
#ifndef READINPUTFILES_H
#define READINPUTFILES_H

#include "Grain.h"
#include "EssentialBoundaryConditions.h"

//namespace Utilities {
	vector<Vector2d>
	readNodePositions(string Positions) {
		string line_position;
		string partial;
		istringstream iss;

		ifstream file_position(Positions.c_str());
		getline(file_position, line_position);
		size_t numberOfNodes = atoi(line_position.c_str());
		//cout << "number of nodes is " << numberOfNodes << endl;

		Vector2d position;
		vector<Vector2d> nodePositions(numberOfNodes);
		for (unsigned int index = 0; index < numberOfNodes; index++) {
    	    getline(file_position, line_position);
    	    iss.str(line_position);
    	    getline(iss, partial, ' ');
    	    position(0) = atof(partial.c_str());
    	    getline(iss, partial, ' ');
    	    position(1) = atof(partial.c_str());
    	    nodePositions[index] = position;
    	    //cout << position(0) << " " << position(1) << endl;
    	    iss.clear();			
		}
		return nodePositions;
	}

	vector< vector<size_t> > 
	readElementMeshes(string Meshes) {
		string line_mesh;
		string partial;
		istringstream iss;

		ifstream file_mesh(Meshes.c_str());
		getline(file_mesh, line_mesh);
		size_t numberOfElements = atoi(line_mesh.c_str());
		//cout << "number of elements is " << numberOfElements << endl; 
		vector<size_t> elementmesh(3);
		vector< vector<size_t> > allElementMeshes(numberOfElements);

		for (unsigned int index = 0; index < numberOfElements; index++) {
			getline(file_mesh, line_mesh);
			iss.str(line_mesh);
			for (unsigned int eindex = 0; eindex < 3; eindex++) {
				getline(iss, partial, ' ');
				elementmesh[eindex] = atoi(partial.c_str());
			}
			iss.clear();
			allElementMeshes[index] = elementmesh;
		}
		return allElementMeshes;
	}

	vector< vector<size_t> >
	readBoundaryConnectivities(string Connectivities) {
		string line_connectivity;
		string partial;
		istringstream iss;

		ifstream file_connectivity(Connectivities.c_str());
		getline(file_connectivity, line_connectivity);
		size_t numberOfElements = atoi(line_connectivity.c_str());
		vector<size_t> singleConnection(2);
		vector< vector<size_t> > Connections(numberOfElements);

		for (unsigned int index = 0; index < numberOfElements; index++) {
			getline(file_connectivity, line_connectivity);
			iss.str(line_connectivity);
			for (unsigned int eindex = 0; eindex < 2; eindex++) {
				getline(iss, partial, ' ');
				singleConnection[eindex] = atoi(partial.c_str());
			}
			iss.clear();
			Connections[index] = singleConnection;
		}
		return Connections;		
	}

	vector<size_t>
	readBoundaryNodeIds(string NodeIds) {
		string line_id;
		ifstream file_id(NodeIds.c_str());
		getline(file_id, line_id);
		size_t numberOfNodes = atoi(line_id.c_str());
		vector<size_t> ids(numberOfNodes);
		for (unsigned int index = 0; index < numberOfNodes; index++) {
			getline(file_id, line_id);
			ids[index] = atoi(line_id.c_str());
		}
		return ids;
	}


	vector<StaticGrain>
	constructGrainsFromFile(string Positions, string Meshes, 
		string Connectivities, string NodeIds, string Pnumbers, 
		string MaterialConstants, const size_t numberOfGrains) {

		string line_position;
		string line_mesh;
		string line_connectivity;
		string line_id;
		string line_number;
		string line_material;

		string partial;
		istringstream iss;

		ifstream file_position(Positions.c_str());
		ifstream file_mesh(Meshes.c_str());
		ifstream file_connectivity(Connectivities.c_str());
		ifstream file_id(NodeIds.c_str());
		ifstream file_number(Pnumbers.c_str());
		ifstream file_material(MaterialConstants.c_str());

		vector<Vector2d>         singleGrainNodePositions;
		vector< vector<size_t> > singleGrainMeshes;
		vector< vector<size_t> > singleGrainBoundaryConnectivities;
		vector<size_t>           singleGrainBoundaryNodeIds;
		vector<Vector3d>         singleGrainProperties;

		vector<size_t>           singleElementMesh(3);
		vector<size_t>           singleBoundaryConnectivity(2);
		Vector2d                 position;

		size_t                   number;
		size_t                   numberOfNodes;
		size_t                   pnumber;
	    size_t                   qpnumber = 3;

		double                   youngsModulus;
		double                   poissonsRatio;
		double                   density;
		double                   frictionCoefficient;

		vector<SimplexTriangle>  singleGrainElements;
		vector<Vector2d>         elementNodePositions(3);

		vector<StaticGrain> grainList(numberOfGrains);                 

		for (unsigned int gindex = 0; gindex < numberOfGrains; gindex++) {
			// do positions for a single grain
			getline(file_position, line_position);
		    number = atoi(line_position.c_str());
		    numberOfNodes = number;
		    singleGrainNodePositions.resize(number);
		    for (unsigned int index = 0; index < number; index++) {
    	        getline(file_position, line_position);
    	        iss.str(line_position);
    	        getline(iss, partial, ' ');
    	        position(0) = atof(partial.c_str());
    	        getline(iss, partial, ' ');
    	        position(1) = atof(partial.c_str());
    	        singleGrainNodePositions[index] = position;
    	        //cout << position(0) << " " << position(1) << endl;
    	        iss.clear();			
		    }

            // do meshes for a single grain
		    getline(file_mesh, line_mesh);
		    number = atoi(line_mesh.c_str());
		    singleGrainElements.resize(number);
		    singleGrainMeshes.resize(number);
		    singleGrainProperties.resize(number);
		    for (unsigned int index = 0; index < number; index++) {
		    	getline(file_mesh, line_mesh);
		    	iss.str(line_mesh);
		    	getline(iss, partial, ' ');
		    	singleElementMesh[0] = atoi(partial.c_str());
		    	getline(iss, partial, ' ');
		    	singleElementMesh[1] = atoi(partial.c_str());
		    	getline(iss, partial, ' ');
		    	singleElementMesh[2] = atoi(partial.c_str());
		    	singleGrainMeshes[index] = singleElementMesh;
		    	getline(iss, partial, ' ');
		    	singleGrainProperties[index](0) = atof(partial.c_str());
		    	getline(iss, partial, ' ');
		    	singleGrainProperties[index](1) = atof(partial.c_str());
		    	getline(iss, partial, ' ');
		    	singleGrainProperties[index](2) = atof(partial.c_str());
		    	iss.clear();
		    }

		    // do connectivities for a single grain
		    getline(file_connectivity, line_connectivity);
		    number = atoi(line_connectivity.c_str());
		    singleGrainBoundaryConnectivities.resize(number);
		    for (unsigned int index = 0; index < number; index++) {
		    	getline(file_connectivity, line_connectivity);
		    	iss.str(line_connectivity);
		    	getline(iss, partial, ' ');
		    	singleBoundaryConnectivity[0] = atoi(partial.c_str());
		    	getline(iss, partial, ' ');
		    	singleBoundaryConnectivity[1] = atoi(partial.c_str());
		    	singleGrainBoundaryConnectivities[index] = singleBoundaryConnectivity;
		    	iss.clear();
		    }

		    // do boundary ids for a single grain
		    getline(file_id, line_id);
		    number = atoi(line_id.c_str());
		    singleGrainBoundaryNodeIds.resize(number);
		    for (unsigned int index = 0; index < number; index++) {
		    	getline(file_id, line_id);
		    	singleGrainBoundaryNodeIds[index] = atoi(line_id.c_str());
		    }


		    // do material constants
		    getline(file_material, line_material);
		    iss.str(line_material);
		    getline(iss, partial, ' ');
		    youngsModulus = atof(partial.c_str());
		    getline(iss, partial, ' ');
		    frictionCoefficient = atof(partial.c_str());
		    //cout << "grain " << gindex << "has friction coe " << frictionCoefficient << endl;
		    iss.clear();

            // now create elements for a single grain
		    for (unsigned int index = 0; index < singleGrainElements.size(); index++) {
		    	elementNodePositions[0] = singleGrainNodePositions[singleGrainMeshes[index][0]];
		    	elementNodePositions[1] = singleGrainNodePositions[singleGrainMeshes[index][1]];
		    	elementNodePositions[2] = singleGrainNodePositions[singleGrainMeshes[index][2]];

		    	MaterialModelLinearElastic * lelastic = new MaterialModelLinearElastic(singleGrainProperties[index](0), singleGrainProperties[index](1), singleGrainProperties[index](2));


		    	SimplexTriangle element = SimplexTriangle(elementNodePositions, singleGrainMeshes[index], lelastic, qpnumber);
		    	singleGrainElements[index] = element;
		    }

		    // do pnumber for a single grain
		    getline(file_number, line_number);
		    pnumber = atoi(line_number.c_str());


            // now construct a grain
		    grainList[gindex] = StaticGrain(gindex, singleGrainElements, numberOfNodes, singleGrainNodePositions,
		    	singleGrainBoundaryNodeIds, singleGrainBoundaryConnectivities, pnumber,youngsModulus,frictionCoefficient);
		    //cout << 2222 << endl;

		}
		return grainList;
	}

	vector<EssentialBC>
	readEssentialBoundaryConditionsFromFile(string essentialBCs) {
		string line_bc;
		string partial;
		istringstream iss;

		ifstream file_bc(essentialBCs.c_str());

		getline(file_bc, line_bc);
		size_t numberOfBCs = atoi(line_bc.c_str());

		size_t grainId;
		size_t nodeIdGrainBased;
		size_t direction;
		double value;
		vector<EssentialBC> EssentialBoundaryConditions(numberOfBCs);

		for (unsigned int index = 0; index < numberOfBCs; index++) {
			getline(file_bc, line_bc);
			iss.str(line_bc);
			getline(iss, partial, ' ');
			grainId = atoi(partial.c_str());
			getline(iss, partial, ' ');
			nodeIdGrainBased = atoi(partial.c_str());
			getline(iss, partial, ' ');
			direction = atoi(partial.c_str());
			getline(iss, partial, ' ');
			value = atof(partial.c_str());
			iss.clear();
			EssentialBC es = EssentialBC(grainId, nodeIdGrainBased, direction, value);
			EssentialBoundaryConditions[index] = es;
		}
		return EssentialBoundaryConditions;
	}


	
//}

#endif  // UTILITIES_H
