#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <iostream>
#include "Vector.hpp"
#include <string.h>
#include <Eigen/Dense>
#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include <opencv/highgui.h>

#include <octomap/octomap.h>
#include <octomap/OcTree.h>
#include <octomap/ColorOcTree.h>


using namespace std;
using namespace octomap;

void Load(unsigned short *depth, unsigned char *color, const char *depth_fname, const char *color_fname){
	char fname[1000];
	//depth
	sprintf(fname, "depth/%s", depth_fname);
	IplImage *freiburg_depth = cvLoadImage(fname, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR);
	char *ptr = freiburg_depth->imageData;
	for(int i = 0; i < 640*480; i++, ptr += sizeof(unsigned short)){
		depth[i*2] = 0;
		memcpy(&depth[i*2+1], ptr, sizeof(unsigned short));
		depth[i*2+1] /= 5;
	}
	cvRelease( (void **)&freiburg_depth);

	//color
	sprintf(fname, "rgb/%s", color_fname);
	IplImage *freiburg_color = cvLoadImage(fname, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR);
	for(int i = 0; i < 640*480; i++){
		color[i*4] = freiburg_color->imageData[i*3];
		color[i*4+1] = freiburg_color->imageData[i*3+1];
		color[i*4+2] = freiburg_color->imageData[i*3+2];
		color[i*4+3] = 0;
	}
	cvRelease( (void **)&freiburg_color);
}

void CreatePointCloud(const char *depth_fname, const char *color_fname, std::vector<Eigen::Vector3f> &points, std::vector<Eigen::Vector3d> &colors){
	static unsigned short depth[640*480*2];
	static unsigned char color[640*480*4];

	Load(depth, color, depth_fname, color_fname);

	for(int i = 0; i < 480; i++){
		for(int j = 0; j < 640; j++){
			
			int depth_index = (i*640+j)*2+1;
			int color_index = (i*640+j)*4;
			if(depth[depth_index] != 0){
				Eigen::Vector3f point;
				Eigen::Vector3d point_color;
				point.x() = (( (float)j-319.5 )*depth[depth_index]*0.001)/525.0;
				point.y() = (( (float)i-239.5 )*depth[depth_index]*0.001)/525.0;
				point.z() = depth[depth_index]*0.001;

				point_color.x() = color[color_index];
				point_color.y() = color[color_index+1];
				point_color.z() = color[color_index+2];
				
				points.push_back(point);
				colors.push_back(point_color);
			}
		}
	}
}

float bitStack(float sd, int red, int green, int blue, float mu){
	int gray = (red + green + blue)/3;
	int temp = (int)(sd / (2.0*mu / 256.0)) + 128;
	temp > 255 ? 255 : temp;
	temp < 0 ? 0 : temp;

	return (float)(temp << 8 | gray);
}

int getColor(int val){
	return val & 0x000000FF;
}

float getSDF(int val, float mu){
	float sd = (float)((val & 0x0000FF00) >> 8);
	sd = (sd - 128) * (2.0*mu / 256.0);
	return sd;
}
float setSDF(int val, float sd, float mu){
	int temp = (int)(sd / (2.0*mu / 256.0)) + 128;
	temp > 255 ? 255 : temp;
	temp < 0 ? 0 : temp;

	val = (val & 0xFFFF00FF) | temp << 8;
	return val;
}

void computeChildCenter (const unsigned int& pos,
    const float& center_offset,
    const point3d& parent_center,
    point3d& child_center) {
  // x-axis
  if (pos & 1) child_center(0) = parent_center(0) + center_offset;
  else     child_center(0) = parent_center(0) - center_offset;

  // y-axis
  if (pos & 2) child_center(1) = parent_center(1) + center_offset;
  else   child_center(1) = parent_center(1) - center_offset;
  // z-axis
  if (pos & 4) child_center(2) = parent_center(2) + center_offset;
  else   child_center(2) = parent_center(2) - center_offset;
}

/// mimics old deprecated behavior to compare against
void getLeafNodesRecurs(std::list<OcTreeVolume>& voxels,
						std::vector<float> &values,
						std::vector<octomap::ColorOcTreeNode::Color> &colors,
    unsigned int max_depth,
    ColorOcTreeNode* node, unsigned int depth,
    const point3d& parent_center, const point3d& tree_center,
    ColorOcTree* tree, bool occupied)
{
  if ((depth <= max_depth) && (node != NULL) ) {
    if (node->hasChildren() && (depth != max_depth)) {

      double center_offset = tree_center(0) / pow( 2., (double) depth+1);
      point3d search_center;

      for (unsigned int i=0; i<8; i++) {
        if (node->childExists(i)) {
		
          computeChildCenter(i, center_offset, parent_center, search_center);
          getLeafNodesRecurs(voxels, values, colors, max_depth, node->getChild(i), depth+1, search_center, tree_center, tree, occupied);

        } // GetChild
      }
    }
    else {
      if (tree->isNodeOccupied(node) == occupied){
        double voxelSize = tree->getResolution() * pow(2., double(16 - depth));
        voxels.push_back(std::make_pair(parent_center - tree_center, voxelSize));
		values.push_back( node->getValue() );
		colors.push_back( node->getColor() );
      }
    }
  }
}

void PointsFromOctree(octomap::ColorOcTree &tree, const char *fname = "cloud.ply"){
	tree.updateInnerOccupancy();

	const unsigned int tree_max_val(32768);
	point3d tree_center;
	  tree_center(0) = tree_center(1) = tree_center(2)
				  = (float) (((double) tree_max_val) * tree.getResolution());

	  
	std::list<OcTreeVolume> voxels;
	std::vector< float > values;
	std::vector< octomap::ColorOcTreeNode::Color > colors;

	getLeafNodesRecurs(voxels, values, colors, 50, tree.getRoot(), 0, tree_center, tree_center, &tree, true);

	FILE *fptr = fopen(fname,"w");
	fprintf(fptr, "ply\n"
					"format ascii 1.0\n"  
					"element vertex %d\n"
					"property float x\n"    
					"property float y\n"
					"property float z\n"
					"property uchar red\n"
					"property uchar green\n"
					"property uchar blue\n"
					"end_header\n", voxels.size());
	
	std::vector<octomap::ColorOcTreeNode::Color>::iterator cItr = colors.begin();
	for(std::list<OcTreeVolume>::iterator itr = voxels.begin(); itr != voxels.end(); itr++, cItr++){
		//unsigned int val = (unsigned int)(*cItr);
		//fprintf(fptr, "%f %f %f %d %d %d\n", itr->first.x(), itr->first.y(), itr->first.z(), (int)(255.0 * (*cItr)), (int)(255.0 * (*cItr)), (int)(255.0 * (*cItr)));
		//int red, green, blue;
		
		//int gray = getColor(val);
		fprintf(fptr, "%f %f %f %d %d %d\n", itr->first.x(), itr->first.y(), itr->first.z(),  cItr->b, cItr->g, cItr->r);
	}

	fclose(fptr);

}

void CreateOctree(const char *depth_fname, const char *color_fname, octomap::ColorOcTree &tree, Eigen::Matrix<float, 3, 3, Eigen::RowMajor> &R, Eigen::Vector3f &T){
	static unsigned short depth[640*480*2];
	static unsigned char color[640*480*4];

	Load(depth, color, depth_fname, color_fname);

	for(int i = 0; i < 480; i ++){
		for(int j = 0; j < 640; j ++){
			
			int depth_index = (i*640+j)*2+1;
			int color_index = (i*640+j)*4;
			if(depth[depth_index] != 0){
				Eigen::Vector3f point;
				Eigen::Vector3d point_color;
				point.x() = (( (float)j-319.5 )*depth[depth_index]*0.001)/525.0;
				point.y() = (( (float)i-239.5 )*depth[depth_index]*0.001)/525.0;
				point.z() = depth[depth_index]*0.001;

				point = R*point + T;

				point_color.x() = color[color_index];
				point_color.y() = color[color_index+1];
				point_color.z() = color[color_index+2];

				int val = (int)(point_color.z()) << 16 | (int)(point_color.y()) << 8 | (int)point_color.x();
				ColorOcTreeNode *node = tree.updateNode(octomap::point3d (point.x(),point.y(),point.z()), true);
				node->setValue(bitStack(0,point_color.x(),point_color.y(),point_color.z(),0));
				node->setColor(point_color.x(), point_color.y(), point_color.z());
				
			}
		}
	}
}

void SavePointCloud(const std::vector<Eigen::Vector3f> &points, const std::vector<Eigen::Vector3d> &colors){
	printf("%d %d\n", points.size(), colors.size());
	FILE *fptr = fopen("output.ply","w");
	fprintf(fptr,	"ply\n"
					"format ascii 1.0\n"
					"element vertex %d\n"
					"property float x\n"
					"property float y\n"
					"property float z\n"
					"property uchar red\n"
					"property uchar green\n"
					"property uchar blue\n"
					"property uchar alpha\n"
					"end_header\n", points.size());
	for(int i = 0; i < points.size(); i++){
		fprintf(fptr, "%f %f %f %d %d %d %d\n", points[i].x(), points[i].y(), points[i].z(), (int)colors[i].x(),(int)colors[i].y(),(int)colors[i].z(),0);
	}
	fclose(fptr);
}




void CreateTSDF(const char *depth_fname, const char *color_fname, octomap::ColorOcTree &tree, Eigen::Matrix<float, 3, 3, Eigen::RowMajor> &R, Eigen::Vector3f &T){
	int size = 512;
	int cellsPerMeter = 128;
	float step = (1.0 / (float)cellsPerMeter) ;

	static unsigned short depth[640*480*2];
	static unsigned char color[640*480*4];
	Load(depth, color, depth_fname, color_fname);

	std::vector<Eigen::Vector3f> points;
	std::vector<Eigen::Vector3d> colors;

	Eigen::Vector3f point;
	Eigen::Vector3f lambda;
	
	float mu = 0.01;
	for(int row = 0; row < size; row ++){
		for(int col = 0; col < size; col ++){
			for(int z = 0; z < size; z ++){

				//row = 256;
				//col = 256;

				point.x() = (float)(col - size/2)*step;
				point.y() = (float)(row - size/2)*step;
				point.z() = (float)(z+cellsPerMeter/2)*step;

				
				float x = point.x()*525 / point.z() + 319.5;
				float y = point.y()*525 / point.z() + 239.5;

				if(x < 0 || x >= 640 || y < 0 || y >= 480)
					continue;

				int index = ((int)y * 640 + (int)x)*2 + 1;
				int color_index = ( (int)y*640 + (int)x )*4;
				float val = (float)depth[index]*0.001;

				lambda.x() = (x-319.5) / 525;
				lambda.y() = (y-239.5) / 525;
				lambda.z() = 1;
				

				float sd = val - 1.0 / (lambda.norm()) * point.norm();
	
				sd > mu ? mu : sd;
				sd < -mu ? -mu : sd;

				point = R*point + T;
				//tree.coordToKey( octomap::point3d (point.x(),point.y(),point.z()) );
				
				ColorOcTreeNode *node = tree.search( octomap::point3d (point.x(),point.y(),point.z()), 0 );
				if(node == NULL){
					if(fabs(sd) < mu){
						int temp = (int)(sd / (2.0*mu / 256.0)) + 128;
						temp > 255 ? 255 : temp;
						temp < 0 ? 0 : temp;

						int red = (int)color[color_index+2];
						int green = (int)color[color_index+1];
						int blue = (int)color[color_index];
						
						octomap::ColorOcTreeNode *new_node = tree.updateNode(octomap::point3d (point.x(),point.y(),point.z()), true);//->setValue(bitStack(sd,red,green,blue,mu));
						new_node->setValue(sd);
						new_node->setColor(blue,green,red);
					}
				}else if(sd > 0){
					//float newVal = (sd + getSDF(node->getValue(), mu))/2;
					float newVal = (sd + node->getValue())/2;
					if( fabs(newVal) < mu ){
						node->setValue(newVal);
					}
					else
						tree.deleteNode(octomap::point3d (point.x(),point.y(),point.z()));
				}
				

				
			}
		}
	}
	//SavePointCloud(points,colors);
	//PointsFromOctree(tree);
}

int main(int argc, char *argv[]){
	if(atoi(argv[1]) == -1){
		system("del  associatedImages.txt");
		system("del  imagesWithTraj.txt");
		system("associate.py depth.txt rgb.txt > associatedImages.txt");
		system("associate.py associatedImages.txt groundTruth.txt > imagesWithTraj.txt");
		exit(1);
	}

	FILE *fptr = fopen("imagesWithTraj.txt","r");
	
	double timestamp;
	char depth_name[1000];
	char color_name[1000];
	float x, y, z;
	float i, j, k, w;
	FILE *trajFptr = fopen("trajectory.txt","w");
	BPVector::Quaternion q_inv;
	int counter = 0;
	Eigen::Vector3f center;

	std::vector<Eigen::Vector3f> points;
	std::vector<Eigen::Vector3d> colors;
	octomap::ColorOcTree tree(0.02);
	octomap::ColorOcTree tsdf_tree(0.02);
	
	float tsdf[512*512*512];

	
	while(1){
		if(fscanf(fptr, "%lf %s %lf %s %lf %f %f %f %f %f %f %f", &timestamp, &depth_name[0], &timestamp, &color_name[0], &timestamp, &x, &y, &z, &i, &j, &k, &w) < 1)
			break;
		BPVector::Quaternion q(w,i,j,k);
		if(counter == 0){
			q_inv = q.inverse();
			center.x() = x;
			center.y() = y;
			center.z() = z;
		
		}
		Eigen::Matrix<float, 3, 3, Eigen::RowMajor> R = (q).RotationMatrix();
		
		printf("%s\n", depth_name);
		fprintf(trajFptr, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", R(0,0), R(0,1), R(0,2), 0,
																			   R(1,0), R(1,1), R(1,2), 0,
																			   R(2,0), R(2,1), R(2,2), 0,
																			   x ,y,z,1.0);
		char *token;
		char *depth_fname;
		char *color_fname;
		char *search = "/";

		// Token will point to "SEVERAL".
		token = strtok(depth_name, search);
		depth_fname = strtok(NULL, search);

		token = strtok(color_name, search);
		color_fname = strtok(NULL, search);
	
		//CreateOctree(depth_fname, color_fname, tree, R, Eigen::Vector3f(x,y,z));
		counter++;
		printf("%d\n", counter);
		//if(counter > 50)
			//break;
		//CreatePointCloud(depth_fname, color_fname, points, colors);
		//SavePointCloud(points, colors);
		CreateTSDF(depth_fname, color_fname, tsdf_tree, R, Eigen::Vector3f(x,y,z));
		if(counter == atoi(argv[1]))
			break;
	}

	FILE *size = fopen("sizes.txt","w");
	fprintf(size, "OCC Tree = %d\nTSDF Tree = %d\n", tree.memoryUsage(), tsdf_tree.memoryUsage());
	fclose(size);

	//PointsFromOctree(tree, "occCloud.ply");
	PointsFromOctree(tsdf_tree, "tsdfCloud.ply");
	fclose(trajFptr);
}