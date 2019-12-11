DCEL = {
	faces: [],
	vertices: [],
	edges: [],
}

tree = {
	triangle: [], //The vertex of the triangle
	face: 0, //At the beginning it points to infinity
	child: [],


}


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
/*
	Returns the leaf that holds the face
*/
function getLeaf(face, node){
	if(!node)
		node = tree;
	if(node.face == face && node.child.length == 0)
		return node;
	else{
		for(var i=0; i<node.child.length; i++){
			nd = getLeaf(face, node.child[i]);
			if(nd != null)
				return nd;
		}
	}
	return null;
}

function addNode(node,triangles, faceID){
	node.child.push({
		triangle: triangles,
		face: faceID,
		child: []
	});
}

function addToDCEL(p, faceIndex, nodeA){ //ADDs a point in a known face and updates the hierarchy
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
			addNode(nodeA,findTriangle(triangle.faces[i]).vertices,triangle.faces[i])
		}

}	else{ //Point on returned edge and edge is obliged to be an internal edge
		//Find the 2 affected triangles
		var triangleA = findTriangle(faceIndex, onEdge);
		var triangleB;
		var faceB;
		var nodeB;

		if(DCEL.edges[onEdge].faceLeft == faceIndex)
			faceB = DCEL.edges[onEdge].faceRight;
		else
			faceB = DCEL.edges[onEdge].faceLeft;

		triangleB = findTriangle(faceB, onEdge);
		nodeB = getLeaf(faceB);
		//Obtenim el node del arbre de la cara B.
		if(nodeB == null){
			console.log("Se perendio esta vaina");
			return;
		}
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
			if(i < 2){
				faceToCheck = faceIndex;

			}else{
				faceToCheck = faceB;
			}
			if(DCEL.edges[quadrilater.boundaryEdges[i]].faceLeft == faceToCheck){
				DCEL.edges[quadrilater.boundaryEdges[i]].faceLeft = quadrilater.faces[i];
				DCEL.edges[quadrilater.boundaryEdges[i]].edgePrevious = quadrilater.internalEdges[i];
			}
			else{
				DCEL.edges[quadrilater.boundaryEdges[i]].faceRight = quadrilater.faces[i];
				DCEL.edges[quadrilater.boundaryEdges[i]].edgeNext = quadrilater.internalEdges[i];
			}
			if ( i < 2 ) {
				addNode(nodeA,findTriangle(quadrilater.faces[i]).vertices,quadrilater.faces[i])
			}
			else {

				addNode(nodeB,findTriangle(quadrilater.faces[i]).vertices,quadrilater.faces[i])
			}
		}

	}
}

function pointInTriangle(p, triangle) {

	var a = orientationTest(DCEL.vertices[triangle[0]], DCEL.vertices[triangle[1]], p);
	var b = orientationTest(DCEL.vertices[triangle[1]], DCEL.vertices[triangle[2]], p);
	var c = orientationTest(DCEL.vertices[triangle[2]], DCEL.vertices[triangle[0]], p);
	if (a == b && b == c) return true; //Interior point!
	else if ((a == "inline" && (b == c)) || (b == "inline" && (a == c)) || (c=="inline" && (a==b))) return true;
	//Point on one edge.
	// no cal revisar si son dos vertex inline, ja que no es donarà el cas de que un punt
	// estigui just a la mateixa pos que un altre, serien el mateix punt
	return false;
}

/*
	Idea no implementada per fer-lo balancejat:
		donar als nodes del arbre, un id arbitrari i que faci que el arbre sigui balancejat
		també tenir un diccionari que faci referència de ID real a ID falsa
		i un altre diccionari al revés ID falsa a ID real
		per tal de tenir un control.
		No esta implementat, però es una idea futura, ja que aleshores el cost seria
		de nlog(n) en lloc de n² en el pitjor cas.
*/
function findHierarchyTriangle(p, h){
	if(!h) h = tree;

	if(pointInTriangle(p,h.triangle)){
		if(h.child.length == 0){ //It's a leaf
			return h;
		}
		else{
			for(var i=0; i<h.child.length; i++){
				var ch = findHierarchyTriangle(p, h.child[i]);
				if(ch != null) return ch;
			}
		}
	}
	return null;
}

function findAndAddPoint(newPoint){ //Finds the leaf in the triangle hierarchy that contains the point and adds it to the DCEL.
	var h = findHierarchyTriangle(newPoint);
	if(h==null)	{
			DCEL.vertices.push({x: null, y: null, edge: null});
	}

	else	addToDCEL(newPoint, h.face, h);
}

function initTree(p1,p2,p3){
	tree.triangle = [p1,p2,p3];
	tree.face = 1;
}

function computeTriangulation(points) {

	var enclosingTriangle = doEnclosingTriangle(points);
	// var enclosingTriangle = [{x:-30,y:0}, {x:60,y:-12}, {x:60,y:45}];

	//Concatenate the enclosing triangle at the beginning of the points array
	points.unshift(enclosingTriangle[0], enclosingTriangle[1], enclosingTriangle[2]);

	//Initialize DCEL with known enclosing triangle
	initDCEL(enclosingTriangle[0], enclosingTriangle[1], enclosingTriangle[2]);
	//Initialize hierarchical triangle structure
	initTree(0, 1, 2);
	console.log(findTriangle(1).vertices)

	for(var i=3; i<points.length; ++i)
		findAndAddPoint(points[i]);

	var lenOutput = DCEL.faces.length-1;
	var outputTriangles = new Array(lenOutput);

	var triangle;
	for(var i = 0; i < lenOutput; ++i){
		triangle = findTriangle(i+1);
		outputTriangles[i] = triangle.vertices;
	}

	return outputTriangles;
}
