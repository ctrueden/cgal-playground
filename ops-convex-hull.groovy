// https://forum.image.sc/t/calculate-mesh-volume-after-convex-hull-3d/63896/4

#@ OpService ops
#@ IOService io
#@ UIService ui
#@ Dataset image
#@ ImagePlus imp
#@ Double(value = 1) isolevel

import java.util.ArrayList;
import java.util.List;

import net.imagej.mesh.Mesh;
import net.imagej.mesh.Meshes;
import net.imagej.mesh.Triangle;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;
import net.imglib2.util.Util;

import org.scijava.vecmath.Point3f;

import customnode.CustomTriangleMesh;
import ij3d.Image3DUniverse;

if (image.getType() instanceof BooleanType) {
	// input image is a binary image
	mask = image;
}
else {
	// convert input image to binary by auto-thresholding
	mask = ops.threshold().otsu(image);
}
println("Mask = " + mask + " [type=" + Util.getTypeFromInterval(mask).getClass().getName() + "]");

Mesh mesh = Meshes.marchingCubes(mask, isolevel);
println("mesh = ${mesh} [${mesh.triangles().size()} triangles, ${mesh.vertices().size()} vertices]")

mesh = Meshes.removeDuplicateVertices(mesh, 0);
println("deduped mesh = ${mesh} [${mesh.triangles().size()} triangles, ${mesh.vertices().size()} vertices]")

mesh = Meshes.simplify(mesh, 0.5, 6);
println("simplified mesh = ${mesh} [${mesh.triangles().size()} triangles, ${mesh.vertices().size()} vertices]")

double meshVolume = ops.geom().size(mesh).getRealDouble();
println("mesh volume = " + meshVolume);

List result = ops.geom().convexHull(mesh);
Mesh hull = (Mesh) result.get(0);
double epsilon = (Double) result.get(1);
println("hull = ${hull} [${hull.triangles().size()} triangles, ${hull.vertices().size()} vertices, epsilon=${epsilon}]")

double hullVolume = ops.geom().size(hull).getRealDouble();
println("hull volume = " + hullVolume);

//outPath = new org.scijava.io.location.FileLocation("/home/curtis/data/3d/t1-head-hull.ply")
//saver = io.getSaver(hull, outPath)
//println("saver = " + saver)

//plyIO = io.getInstances()[1]
//println("plyIO = " + plyIO)
//println(plyIO.supportsSave(hull, outPath))

//io.save(hull, outPath)
//return

// Display original image plus plus mesh in 3D Viewer
List<Point3f> p = new ArrayList<>();
for (Triangle t in mesh.triangles()) {
	p.add(new Point3f(t.v0xf(), t.v0yf(), t.v0zf()));
	p.add(new Point3f(t.v1xf(), t.v1yf(), t.v1zf()));
	p.add(new Point3f(t.v2xf(), t.v2yf(), t.v2zf()));
}
CustomTriangleMesh ctm = new CustomTriangleMesh(p);
println("Volume according to 3D Viewer: ${ctm.getVolume()}");
final Image3DUniverse univ = new Image3DUniverse();
univ.addVoltex(imp, 1);
univ.addCustomMesh(ctm, "Mesh");
univ.show();
