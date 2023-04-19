// Iterate all points inside a 3D polyhedron.

// We use CGAL's Polygon Mesh Processing (PMP):
//   https://doc.cgal.org/latest/Polygon_mesh_processing/
// and/or CGAL's 3D Polyhedral Surface:
//   https://doc.cgal.org/latest/Polyhedron/
// and/or CGAL's 3D Mesh Generation:
//   https://doc.cgal.org/latest/Mesh_3/index.html

//---------------------------------------------
// Read STL file into an imagej-mesh Mesh.

import net.imagej.mesh.naive.NaiveDoubleMesh
import java.io.File
import net.imagej.mesh.io.stl.STLMeshIO

mesh = new NaiveDoubleMesh()
stlFile = new File("/home/curtis/data/3d/t1-head-hull-binary.stl")
new STLMeshIO().read(mesh, stlFile)
println("Loaded mesh: " + mesh)
println("Mesh vertex count: " + mesh.vertices().size())
verts = mesh.vertices()
tris = mesh.triangles()

//---------------------------------------------
// Copy mesh into CGAL data structures.
//
// See:
// https://github.com/CGAL/cgal-swig-bindings/blob/v5.5.2.post202303131425/examples/java/test_pmp.java

import CGAL.Kernel.Point_3
import CGAL.Polygon_mesh_processing.Int_Vector
import CGAL.Polygon_mesh_processing.Point_3_Vector
import CGAL.Polygon_mesh_processing.Polygon_Vector

Point_3_Vector points = new Point_3_Vector();
points.reserve(verts.size());
for (i=0; i<verts.size(); i++) {
	points.add( new Point_3(verts.x(i), verts.y(i), verts.z(i)) );
}
Polygon_Vector polygons = new Polygon_Vector();
for (t=0; t<tris.size(); t++) {
	Int_Vector polygon = new Int_Vector();
	polygon.reserve(3);
	polygon.add((int) tris.vertex0(t));
	polygon.add((int) tris.vertex1(t));
	polygon.add((int) tris.vertex2(t));
	polygons.add(polygon);
}
println("Copied mesh data into CGAL vectors")

// Copy mesh into primitive arrays in CGAL's desired order.
// Unfortunately, there is no orient_polygon_soup for this
// structure, because some points might be duplicated, and
// primitive arrays cannot be dynamically resized.
/*
double[] points = new double[3 * mesh.vertices().size()]
int[] polygons = new int[3 * mesh.triangles().size()]
for (i=0; i<mesh.vertices().size(); i++) {
	points[3*i+0] = mesh.vertices().x(i)
	points[3*i+1] = mesh.vertices().y(i)
	points[3*i+2] = mesh.vertices().z(i)
}
for (t=0; t<mesh.triangles().size(); t++) {
	polygons[3*t+0] = mesh.triangles().vertex0(t)
	polygons[3*t+1] = mesh.triangles().vertex1(t)
	polygons[3*t+2] = mesh.triangles().vertex2(t)
}
println("Copied mesh into arrays")
*/

//---------------------------------------------
// try to consistently orient a soup of polygons in 3D space.
//
// When it is not possible to produce a combinatorial manifold surface, some
// points are duplicated. Because a polygon soup does not have any connectivity
// (each point has as many occurrences as the number of polygons it belongs
// to), duplicating one point (or a pair of points) amounts to duplicate the
// polygon to which it belongs.
//
// These points are either an endpoint of an edge incident to more than two
// polygons, an endpoint of an edge between two polygons with incompatible
// orientations (during the re-orientation process), or more generally a point
// p at which the intersection of an infinitesimally small ball centered at p
// with the polygons incident to it is not a topological disk.
//
// Not sure we really need to do this, but described here:
// https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__orientation__grp.html#gad380465ee62d858d27fab4cfda6c1764

import CGAL.Polygon_mesh_processing.CGAL_Polygon_mesh_processing

println("Executing orient_polygon_soup")
CGAL_Polygon_mesh_processing.orient_polygon_soup(points, polygons)
println("Oriented!")

//---------------------------------------------
// Convert polygon soup to polygon mesh -- i.e. a CGAL Polyhedron
// https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#title43

import CGAL.Polyhedron_3.Polyhedron_3

Polyhedron_3 P = new Polyhedron_3();
println("Executing polygon_soup_to_polygon_mesh")
CGAL_Polygon_mesh_processing.polygon_soup_to_polygon_mesh(points, polygons, P)
println("Constructed CGAL polyhedron: " + P)
println("Polyhedron vertex count: " + P.size_of_vertices())

//---------------------------------------------
// Remove isolated vertices.
// Not sure we really need to do this, but seen here:
// https://github.com/CGAL/cgal-swig-bindings/blob/v5.5.2.post202303131425/examples/java/test_pmp.java

CGAL_Polygon_mesh_processing.remove_isolated_vertices(P);
println("Removed isolated vertices")
println("New vertex count: " + P.size_of_vertices())

//---------------------------------------------
/* NB: Function is not available in Java SWIG bindings.
// Orient to bound the volume.
// Not sure we really need to do this, but seen here:
// https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#title45
println("Orienting to bound a volume")
CGAL_Polygon_mesh_processing.orient_to_bound_a_volume(P)
println("Done!")
*/

//---------------------------------------------
// Sanity check.
println("Is it closed? " + P.is_closed())

// Darn. It needs to be closed, or interior/exterior has no meaning, and CGAL fails.
// It's confusing, because according to the CGAL docs
// (https://doc.cgal.org/latest/Polyhedron/index.html):
// "If the surface is closed we call it a polyhedron"
// But we can construct a Polyhedron_3 that is *not* closed, obviously.

// How to fix it? Let's try triangulating the faces.
// This doesn't actually do anything, because our faces are all already triangles.
// This function is, IIUC, actually just for making non-triangle faces into triangles (e.g. cutting squares in half).
// But I include it here because it's something I tried (unsuccessfully) to make things work.
println("Triangulating faces")
CGAL_Polygon_mesh_processing.triangulate_faces(P);
println("Is it closed now? " + P.is_closed())

// Argh, not sure why this thing isn't closed... it is a saved PLY file of a convex hull computed by imagej-mesh.
// As a workaround, let's recompute the convex hull again, this time using CGAL.

import CGAL.Convex_hull_3.CGAL_Convex_hull_3

println("Recomputing convex hull...")
Polyhedron_3 hull = new Polyhedron_3();
CGAL_Convex_hull_3.convex_hull_3(points.iterator(), hull);
println("convex hull has " + hull.size_of_vertices() + " vertices");
//println("is strongly convex? " + CGAL_Convex_hull_3.is_strongly_convex_3(hull));
println("Is the hull closed? " + hull.is_closed())

//---------------------------------------------
// What's the bounding box of our data?
xMin = yMin = zMin = Double.POSITIVE_INFINITY
xMax = yMax = zMax = Double.NEGATIVE_INFINITY
for (i=0; i<verts.size(); i++) {
	x = verts.x(i)
	y = verts.y(i)
	z = verts.z(i)
	if (x < xMin) xMin = x
	if (y < yMin) yMin = y
	if (z < zMin) zMin = z
	if (x > xMax) xMax = x
	if (y > yMax) yMax = y
	if (z > zMax) zMax = z
}
println("Bounding box: (${xMin}, ${yMin}, ${zMin}) -> (${xMax}, ${yMax}, ${zMax})")

//---------------------------------------------
// Finally, let's determine if a point is inside our polyhedron.
// Code adapted from:
// https://github.com/CGAL/cgal-swig-bindings/blob/v5.5.2.post202303131425/examples/java/test_pmp.java#L334-L345
// See also:
// https://doc.cgal.org/latest/Polygon_mesh_processing/Polygon_mesh_processing_2point_inside_example_8cpp-example.html#_a5

import CGAL.Polygon_mesh_processing.Side_of_triangle_mesh

cx = (xMin + xMax) / 2
cy = (yMin + yMax) / 2
cz = (zMin + zMax) / 2
println("Center point = (${cx}, ${cy}, ${cz})")
Side_of_triangle_mesh sotm = new Side_of_triangle_mesh(hull); // <-- This fails if our original P is used. :-(
bs = sotm.bounded_side(new Point_3(cx, cy, cz))
println("Bounded side = ${bs}")
// The three possibilities are:
// - Bounded_side.ON_BOUNDED_SIDE
// - Bounded_side.ON_BOUNDARY
// - Bounded_side.ON_UNBOUNDED_SIDE

//---------------------------------------------
// Alternate approach, using from:
// https://www.codefull.net/2016/03/cgal-point-in-polyhedron-algorithm/
/*
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef CGAL::Side_of_triangle_mesh<Polyhedron, K> Point_inside;

def pointInside(polyhedron, pt) {
    // Construct AABB tree with a KdTree
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.accelerate_distance_queries();
    // Initialize the point-in-polyhedron tester
    Point_inside inside_tester(tree);

    // Determine the side and return true if inside!
    return inside_tester(query) == CGAL::ON_BOUNDED_SIDE;
}
*/

//---------------------------------------------
// Another potential avenue seen mentioned on Stack Overflow
// (https://stackoverflow.com/a/34806135/1207769) is the
// is_in_complex(Cell_handle c) method of MeshComplex_3InTriangulation_3.
// See:
// https://doc.cgal.org/latest/Mesh_3/classMeshComplex__3InTriangulation__3.html#a7b58574d2c8adfd81fa3b4f6eab1bbb5
// and test_mesh_3 in the Java examples.

// TODO: Try using this.

//---------------------------------------------
// Tangent: what about slicing a polyhedron to the XY plane at a particular Z?
// There might be hints at:
// https://computergraphics.stackexchange.com/questions/7647/projection-of-a-polyhedron-on-xy-plane-with-cgal

//---------------------------------------------
// OK, let's do it! Let's display all points interior to our polyhedron.
// Using the 3D Viewer.

import ij3d.Image3DUniverse;

// First, add the mesh itself.

Image3DUniverse univ = new Image3DUniverse();

//univ.addVoltex(imp, 1);

import customnode.CustomTriangleMesh;
import org.scijava.vecmath.Color3f;

List<Point3f> pts_3dv = new ArrayList<>();
for (t in mesh.triangles()) {
	pts_3dv.add(new Point3f(t.v0xf(), t.v0yf(), t.v0zf()));	
	pts_3dv.add(new Point3f(t.v1xf(), t.v1yf(), t.v1zf()));	
	pts_3dv.add(new Point3f(t.v2xf(), t.v2yf(), t.v2zf()));	
}
CustomTriangleMesh mesh_3dv = new CustomTriangleMesh(pts_3dv);
println("Mesh volume according to 3D Viewer: ${mesh_3dv.getVolume()}");
mesh_3dv.setTransparency(0.5);
mesh_3dv.setColor(new Color3f(java.awt.Color.magenta));
univ.addCustomMesh(mesh_3dv, "Mesh");

// Now, compute and add the points.

import CGAL.Kernel.Bounded_side;
import customnode.CustomPointMesh;
import org.scijava.vecmath.Point3f;

println("Building 3D view of interior integer points")
List<Point3f> vPoints = new ArrayList<>();
int ixMin = Math.ceil(xMin), ixMax = Math.floor(xMax);
int iyMin = Math.ceil(yMin), iyMax = Math.floor(yMax);
int izMin = Math.ceil(zMin), izMax = Math.floor(zMax);
int step = 3;
Point_3 pt = new Point_3()
for (z=izMin; z<=izMax; z += step) {
	println("z -> ${z} of ${izMax - izMin + 1}")
	for (y=iyMin; y<=iyMax; y += step) {
		for (x=ixMin; x<=ixMax; x += step) {
			pt.set_coordinates(x, y, z);
			boolean interior = false;
//			try {
				interior = sotm.bounded_side(pt) != Bounded_side.ON_UNBOUNDED_SIDE;
//			}
//			catch (Exception ex) { } // CGAL.Java.SWIGCGALException can happen.
			if (interior) vPoints.add(new Point3f(x, y, z));
		}
	}
}
univ.addCustomMesh(new CustomPointMesh(vPoints), "Interior Points")
println("Done building points: ${vPoints.size()}")

univ.show();
