/**
 TODO Replace this by your own, correct, triangulation function
 Triangles should be return as arrays of array of indexes
 e.g., [[1,2,3],[2,3,4]] encodes two triangles, where the indices are relative to the array points
**/

DCEL = {
	faces: [],
	vertices: [],
	edges: [],
}

// MODELS:

// face = {
// 	edge: 0
// }
// vertex = {
// 	x: 0.0,
// 	y: 0.0,
// 	edge: 0
// }
// edge = {
// 	vertexBegin: 0,
// 	vertexEnd: 1,
// 	faceLeft: 0,
// 	faceRight: 1,
// 	edgePrevious: 0,
// 	edgeNext: 1
// }

function computeMinMaxX(points){ //Returns min and max points by x coordinate and returns them a little bit further
	var min = points[0];
	var max = points[0];
	
	for(var i = 1; i < points.length; ++i){
		if(points[i].x < min.x) min = points[i];
		else if(points[i].x > max.x) max = points[i];				
	}
	return {min: {x: min.x - (max.x-min.x)/2, y: min.y}, max: {x: max.x + (max.x-min.x)/4, y: max.y}};
}

function computeSupportingLines(p, points){
	var left = points[0];
	var right = points[0];

	for(var i = 1; i < points.length; ++i){
		if(orientationTest(p, left, points[i]) == "left") left = points[i];
		else if(orientationTest(p, right, points[i]) == "right") right = points[i];				
	}
	return {left: left, right: right};
}

function crossProduct(v1, v2){
	return v1[0]*v2[1]-(v2[0]*v1[1]);
}

function orientationTest(p, q, r){
	var v1 = [q.x-p.x, q.y-p.y];
	var v2 = [r.x-p.x, r.y-p.y];
	var det = crossProduct(v1, v2);
	if(det > 0) return "left";
	else if(det < 0) return "right";
	else return "inline";
}

function findTriangle(faceIndex, firstEdgeIndex){ //Through the face index we find all the edges forming the triangle

	var edges = new Array(3); 
	var vertices = new Array(3); 

	if(firstEdgeIndex)
		edges[0] = firstEdgeIndex; //First edge given by input
	else
		edges[0] = DCEL.faces[faceIndex].edge; //First edge referenciat per la face on es troba p
	
	//Edge we are looking at
	var actualEdge = DCEL.edges[edges[0]]; //We take the 1st edge of the triangle from the face
	
	//2nd edge of the triangle and first 2 vertices
	if(actualEdge.faceLeft == faceIndex){
		edges[2] = actualEdge.edgePrevious;
		vertices[0] = actualEdge.vertexBegin;
		vertices[1] = actualEdge.vertexEnd;
	}
	else{
		edges[2] = actualEdge.edgeNext;
		vertices[0] = actualEdge.vertexEnd;
		vertices[1] = actualEdge.vertexBegin;
	}
	actualEdge =  DCEL.edges[edges[2]]; //Look at next edge

	//3rd edge of the triangle and last vertex
	if(actualEdge.faceLeft == faceIndex){
		edges[1] = actualEdge.edgePrevious;
		vertices[2] = actualEdge.vertexBegin;
	}
	else{
		edges[1] = actualEdge.edgeNext;
		vertices[2] = actualEdge.vertexEnd;
	}

	return {edges: edges, vertices: vertices}; //Face, Edges and vertices in CCW and in the same order starting from the same vertex
}


function pointOnEdge(p, edges, f) { //Shall return the number of the edge on which p is placed, or -1 if p is inside the triangle and not on an edge

	//Do all orientation tests of p with the 3 edges just to test if determinant is 0 so edge orientation is not important
	var orientations = edges.map((edge)=>{
		var actualEdge = DCEL.edges[edge];
		return orientationTest(DCEL.vertices[actualEdge.vertexBegin], DCEL.vertices[actualEdge.vertexEnd], p);
	});

	//Check for p on an edge else return -1
	for(var i = 0; i < 3; ++i)
		if(orientations[i] == "inline") return edges[i];
	return -1;
}

function pointInTriangle(p, edges, f) { //Shall return true if point is in triangle or boundary, else returns false

	//Do all orientation tests of p with the 3 edges, all orientation tests are done in counter-clockwise so in the triangle means all edges have p on left
	var orientations = edges.map((edge)=>{
		var actualEdge = DCEL.edges[edge];
		if(actualEdge.faceLeft == f) return orientationTest(DCEL.vertices[actualEdge.vertexBegin], DCEL.vertices[actualEdge.vertexEnd], p);
		return orientationTest(DCEL.vertices[actualEdge.vertexEnd], DCEL.vertices[actualEdge.vertexBegin], p);
	});

	if(orientations[0] == "left" && orientations[1] == "left" && orientations[2] == "left") 
		return true; //Point inside triangle
	else if((orientations[0] == "inline" && orientations[1] == orientations[2]) || (orientations[1] == "inline" && orientations[2] == orientations[0]) || (orientations[2] == "inline" && orientations[0] == orientations[1]))
		return true; //Point on triangle boundary
	else return false; //Point not in triangle
}

function classifyIntersection(a1, a2, b1, b2) {

	var throughVertex = null;
	
	var orientations = [orientationTest(a1, a2, b1), orientationTest(a1, a2, b2),  orientationTest(b1, b2, a1),  orientationTest(b1, b2, a2)];
	var zeros = orientations.reduce((prev, curr)=>{if(curr == "inline") ++prev; return prev;}, 0);

	if(zeros == 1){
		//Endpoint and interior point
		if(orientations[0] == 'inline' && orientations[2] != orientations[3]){
			intersectionType = "INTERIOR POINT";
			throughVertex = 0;
		}
		else if(orientations[1] == 'inline' && orientations[2] != orientations[3]){
			intersectionType = "INTERIOR POINT";
			throughVertex = 1;
		} 
		else{ // No intersection
			intersectionType = "NONE";
		}
	}
	else if(orientations[0] != orientations[1] && orientations[2] != orientations[3]){ // CLASSIC Intersection on interior points
		intersectionType = "CLASSIC";
	}
	else{ //NO INTERSECTION
		intersectionType = "NONE";
	}
		
	// Return object with two fields: a numbered type, and a description
	return {type: intersectionType, vertex: throughVertex} ;
}

function doEnclosingTriangle(points){
	var minMax = computeMinMaxX(points);
	var p1 = minMax.min;
	var supLines = computeSupportingLines(minMax.min, points.slice(1));
	var p2 = {x: minMax.max.x, y: (minMax.max.x-p1.x)*((supLines.right.y-p1.y)/(supLines.right.x-p1.x))+p1.y-(minMax.max.x-minMax.min.x)/8};
	var p3 = {x: minMax.max.x, y: (minMax.max.x-p1.x)*((supLines.left.y-p1.y)/(supLines.left.x-p1.x))+p1.y+(minMax.max.x-minMax.min.x)/8};
	return [p1, p2, p3];
}

function initDCEL(p0, p1, p2){
	//ADD first 2 faces
	DCEL.faces.push(
		{edge: 0}, //First face of all is infinity
		{edge: 0});

	//ADD first 3 edges
	DCEL.edges.push(
		{
			vertexBegin: 0, vertexEnd: 1,
			faceLeft: 1, faceRight: 0,
			edgePrevious: 2, edgeNext: 1
		},
		{
			vertexBegin: 1, vertexEnd: 2,
			faceLeft: 1, faceRight: 0,
			edgePrevious: 0, edgeNext: 2
		},
		{
			vertexBegin: 2, vertexEnd: 0,
			faceLeft: 1, faceRight: 0,
			edgePrevious: 1, edgeNext: 0
		});

	//ADD first 3 vertices
	DCEL.vertices.push(
		{x: p0.x, y: p0.y, edge: 0},
		{x: p1.x, y: p1.y, edge: 1},
		{x: p2.x, y: p2.y, edge: 2});
}

function addToDCEL(p, faceIndex){ //ADDs a point in a known face
	//ADD vertex
	DCEL.vertices.push({x: p.x, y: p.y, edge: DCEL.edges.length}); //We assign the next new edge that we will add since we know it's incident

	//Determine edges of the triangle
	var triangle = findTriangle(faceIndex); //Find edges and points of the triangle in CCW

	//Determine wether point is on boundary or in triangle, if in boundary in what edge
	var onEdge = pointOnEdge(p, triangle.edges, faceIndex);

	if(onEdge == -1){ //Point not on edge
		//ADD 2 new faces
		DCEL.faces.push(
			{edge: triangle.edges[1]}, 
			{edge: triangle.edges[2]});
		
		//ADD triangle faces in triangle in CCW and in the same order for faster understanding and accessibility
		triangle.faces = [faceIndex, DCEL.faces.length-2, DCEL.faces.length-1];

		//ADD what will be the new edges indexs in CCW and in the same order for easy access
		triangle.newEdges = [DCEL.edges.length, DCEL.edges.length+1, DCEL.edges.length+2];

		//ADD 3 new edges in CCW and the same order as faces
		DCEL.edges.push(
		{
			vertexBegin: DCEL.vertices.length-1, vertexEnd: triangle.vertices[0],
			faceLeft: triangle.faces[0], faceRight: triangle.faces[2],
			edgePrevious: DCEL.edges.length+1, edgeNext: triangle.edges[2]
		},
		{
			vertexBegin: DCEL.vertices.length-1, vertexEnd: triangle.vertices[1],
			faceLeft: triangle.faces[1], faceRight: triangle.faces[0],
			edgePrevious: DCEL.edges.length+2, edgeNext: triangle.edges[0]
		},
		{
			vertexBegin: DCEL.vertices.length-1, vertexEnd: triangle.vertices[2],
			faceLeft: triangle.faces[2], faceRight: triangle.faces[1],
			edgePrevious: DCEL.edges.length, edgeNext: triangle.edges[1]
		});

		//UPDATE 3 existing edges
		for(var i = 0; i < 3; ++i){
			if(DCEL.edges[triangle.edges[i]].faceLeft == faceIndex){
				DCEL.edges[triangle.edges[i]].edgePrevious = triangle.newEdges[i];
				DCEL.edges[triangle.edges[i]].faceLeft = triangle.faces[i];
			}
			else{
				DCEL.edges[triangle.edges[i]].edgeNext = triangle.newEdges[i];
				DCEL.edges[triangle.edges[i]].faceRight = triangle.faces[i];
			}
		}
	}
	else{ //Point on returned edge and edge is obliged to be an internal edge

		//Find the 2 affected triangles
		var triangleA = findTriangle(faceIndex, onEdge);
		var triangleB;
		var faceB;
		if(DCEL.edges[onEdge].faceLeft == faceIndex)
			faceB = DCEL.edges[onEdge].faceRight;
		else
			faceB = DCEL.edges[onEdge].faceLeft;

		triangleB = findTriangle(faceB, onEdge);
		
		//Define Quadrilater union of the 2 triangles in CCW order starting by edge pointed by onEdge
		var quadrilater = {
			faces: [faceIndex, DCEL.faces.length, DCEL.faces.length+1, faceB],
			boundaryEdges: [triangleA.edges[1], triangleA.edges[2], triangleB.edges[1], triangleB.edges[2]],
			vertices:[triangleA.vertices[1], triangleA.vertices[2], triangleB.vertices[1], triangleB.vertices[2]],
			internalEdges: [onEdge, DCEL.edges.length, DCEL.edges.length+1, DCEL.edges.length+2]
		}

		//ADD 2 new null faces
		DCEL.faces.push({edge: null},{edge: null})

		//ADD 3 new null edges
		DCEL.edges.push(
		{
			vertexBegin: null, vertexEnd: null,
			faceLeft: null, faceRight: null,
			edgePrevious: null, edgeNext: null
		},
		{
			vertexBegin: null, vertexEnd: null,
			faceLeft: null, faceRight: null,
			edgePrevious: null, edgeNext: null
		},
		{
			vertexBegin: null, vertexEnd: null,
			faceLeft: null, faceRight: null,
			edgePrevious: null, edgeNext: null
		});

		//Overwrite all elements in quadrilater with the corresponding data
		for(var i = 0; i < 4; ++i){
			//Face
			DCEL.faces[quadrilater.faces[i]].edge = quadrilater.boundaryEdges[i];

			//Vertex
			DCEL.vertices[quadrilater.vertices[i]].edge = quadrilater.internalEdges[i];

			//Internal Edge
			DCEL.edges[quadrilater.internalEdges[i]] = {
				vertexBegin: DCEL.vertices.length-1,
				vertexEnd: quadrilater.vertices[i],
				faceLeft: quadrilater.faces[i],
				faceRight: quadrilater.faces[(i+3)%4],
				edgePrevious: quadrilater.internalEdges[(i+1)%4],
				edgeNext: quadrilater.boundaryEdges[(i+3)%4]
			}

			//Boundary Edge only needs to change a face and an edge
			var faceToCheck;
			if(i < 2)
				faceToCheck = faceIndex;
			else
				faceToCheck = faceB;

			if(DCEL.edges[quadrilater.boundaryEdges[i]].faceLeft == faceToCheck){
				DCEL.edges[quadrilater.boundaryEdges[i]].faceLeft = quadrilater.faces[i];
				DCEL.edges[quadrilater.boundaryEdges[i]].edgePrevious = quadrilater.internalEdges[i];
			}
			else{
				DCEL.edges[quadrilater.boundaryEdges[i]].faceRight = quadrilater.faces[i];
				DCEL.edges[quadrilater.boundaryEdges[i]].edgeNext = quadrilater.internalEdges[i];
			}
		}
	}
}



function findAndAddPoint(newPoint, auxiliaryPoint){
	var actualFace = null; //Face we are checking at the moment
	var actualVertex = Object.assign({}, auxiliaryPoint); //If we are now in a vertex then != null
	actualVertex.edges = [];
	actualVertex.faces = [];
	var fromEdge = null; //Edge where we come from
	var visitedFaces = []; //Visited faces along the search to not repeat them in the process

	var containingFace = null; //Face containing the new point

	while(containingFace == null){
		if(actualVertex != null){ //If we currently arrived through a vertex
			//Check all faces and edges incident to vertex
			var initialEdge = actualVertex.vertex.edge;
			var actualEdge;
			if(DCEL.edges[initialEdge].vertexBegin == actualVertex.vertexIndex){ 
				if(visitedFaces.indexOf(DCEL.edges[initialEdge].faceLeft) == -1){
					actualVertex.faces.push(DCEL.edges[initialEdge].faceLeft);
					actualVertex.edges.push(initialEdge);
				}
				actualEdge = DCEL.edges[initialEdge].edgePrevious;
			}
			else{
				if(visitedFaces.indexOf(DCEL.edges[initialEdge].faceRight) == -1){
					actualVertex.faces.push(DCEL.edges[initialEdge].faceRight);
					actualVertex.edges.push(initialEdge);
				}
				actualEdge = DCEL.edges[initialEdge].edgeNext;
			}

			while(actualEdge != initialEdge){ //While we don't check the same edge twice
				if(DCEL.edges[actualEdge].vertexBegin == actualVertex.vertexIndex){ 
					if(visitedFaces.indexOf(DCEL.edges[actualEdge].faceLeft) == -1){
						actualVertex.faces.push(DCEL.edges[actualEdge].faceLeft);
						actualVertex.edges.push(actualEdge);
					}
					actualEdge = DCEL.edges[actualEdge].edgePrevious;
				}
				else{
					if(visitedFaces.indexOf(DCEL.edges[actualEdge].faceRight) == -1){
						actualVertex.faces.push(DCEL.edges[actualEdge].faceRight);
						actualVertex.edges.push(actualEdge);
					}
					actualEdge = DCEL.edges[actualEdge].edgeNext;
				}
			}

			var triangles = new Array(actualVertex.faces.length);
			for(var i = 0; i < actualVertex.faces.length; ++i){
				triangles[i] = findTriangle(actualVertex.faces[i], actualVertex.edges[i]);
				if(pointInTriangle(newPoint, triangles[i].edges, actualVertex.faces[i])){
					containingFace = actualVertex.faces[i];
					break;
				}
			}

			if(containingFace == null){ //Point is not in incident triangles
				visitedFaces = visitedFaces.concat(actualVertex.faces); //We note down the faces we already checked

				//We will only check the opposite edges, those are the only ones necessary
				var edgesToCheck = triangles.map((triangle)=>{ return triangle.edges[1]; });

				for(var i = 0; i < edgesToCheck.length; ++i){ //Check the edges for intersections
					var intersection = classifyIntersection(newPoint, auxiliaryPoint.vertex, DCEL.vertices[DCEL.edges[edgesToCheck[i]].vertexBegin], DCEL.vertices[DCEL.edges[edgesToCheck[i]].vertexEnd]);
					
					if(intersection.type == "CLASSIC"){
						fromEdge = edgesToCheck[i];
						if(DCEL.edges[fromEdge].faceLeft == actualVertex.faces[i])
							actualFace = DCEL.edges[fromEdge].faceRight;
						else
							actualFace = DCEL.edges[fromEdge].faceLeft;
						actualVertex = null;
						break;
					}
					else if(intersection.type == "INTERIOR POINT"){

						if(intersection.vertex == 0){
							actualVertex = {
								vertexIndex: DCEL.edges[edgesToCheck[i]].vertexBegin,
								vertex: DCEL.vertices[DCEL.edges[edgesToCheck[i]].vertexBegin]
							};
						}
						else{
							actualVertex = {
								vertexIndex: DCEL.edges[edgesToCheck[i]].vertexEnd,
								vertex: DCEL.vertices[DCEL.edges[edgesToCheck[i]].vertexEnd]
							};
						}
						break;
					}
				}
			}
		}
		else{ //Check actualFace, simple case
			var triangle = findTriangle(actualFace, fromEdge);
			if(pointInTriangle(newPoint, triangle.edges, actualFace)){
				containingFace = actualFace;
			}
			else{
				visitedFaces.push(actualFace);

				var edgesToCheck = triangle.edges.slice(-2);

				for(var i = 0; i < 2; ++i){ //Check the edges for intersections
					var intersection = classifyIntersection(newPoint, auxiliaryPoint.vertex, DCEL.vertices[DCEL.edges[edgesToCheck[i]].vertexBegin], DCEL.vertices[DCEL.edges[edgesToCheck[i]].vertexEnd]);
					
					if(intersection.type == "CLASSIC"){

						fromEdge = edgesToCheck[i];
						if(DCEL.edges[fromEdge].faceLeft == actualFace)
							actualFace = DCEL.edges[fromEdge].faceRight;
						else
							actualFace = DCEL.edges[fromEdge].faceLeft;
						break;
					}
					else if(intersection.type == "INTERIOR POINT"){

						if(intersection.vertex == 0){
							actualVertex = {
								vertexIndex: DCEL.edges[edgesToCheck[i]].vertexBegin,
								vertex: DCEL.vertices[DCEL.edges[edgesToCheck[i]].vertexBegin]
							};
						}
						else{
							actualVertex = {
								vertexIndex: DCEL.edges[edgesToCheck[i]].vertexEnd,
								vertex: DCEL.vertices[DCEL.edges[edgesToCheck[i]].vertexEnd]
							};
						}
						fromEdge = null;
						actualFace = null;
						break;
					}
				}
			}
		}
	}

	addToDCEL(newPoint, containingFace);
}

function computeTriangulation(points) {

	var enclosingTriangle = doEnclosingTriangle(points);
	// var enclosingTriangle = [{x:-30,y:0}, {x:60,y:-12}, {x:60,y:45}];
	
	//Concatenate the enclosing triangle at the beginning of the points array
	points.unshift(enclosingTriangle[0], enclosingTriangle[1], enclosingTriangle[2]);

	//Initialize DCEL with known enclosing triangle
	initDCEL(enclosingTriangle[0], enclosingTriangle[1], enclosingTriangle[2]);

	//ADD first point
	addToDCEL(points[3], 1);

	
	//Take first point as auxiliary and all i'ts information
	auxiliaryPoint = {
		vertexIndex: 3,
		vertex: DCEL.vertices[3],
	}

	//TODO: find point in space through auxiliary
	// - Look if the point is in any of the incident triangles
	// - If not until we find in which triangle it is located we follow the line between point and auxiliary
	//   checking all the triangles that intersect with the line

	for(var i = 4; i < points.length; ++i){
		findAndAddPoint(points[i], auxiliaryPoint);
	}

	//READ triangles from DCEL
	var outputTriangles = new Array(DCEL.faces.length-1);

	var triangle; 
	for(var i = 1; i < DCEL.faces.length; ++i){
		triangle = findTriangle(i);
		outputTriangles[i-1] = triangle.vertices;
	}

	return outputTriangles;
}



