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



class TSDFOctree{
public:
	TSDFOctree(float resolution){
		tsdf = new octomap::ColorOcTree(resolution);
	}

	void IntegrateData(unsigned short *depth, unsigned char *color, Eigen::Matrix<float, 3, 3, Eigen::RowMajor> &R, Eigen::Vector3f &T){
		int size = 512;
		int cellsPerMeter = 128;
		float step = (1.0 / (float)cellsPerMeter) ;

		std::vector<Eigen::Vector3f> points;
		std::vector<Eigen::Vector3d> colors;

		Eigen::Vector3f point;
		Eigen::Vector3f lambda;
	
		float mu = 0.01;
		for(int row = 0; row < size; row ++){
			for(int col = 0; col < size; col ++){
				for(int z = 0; z < size; z ++){
					//Convert from grid coordinate to 3D point
						point.x() = (float)(col - size/2)*step;
						point.y() = (float)(row - size/2)*step;
						point.z() = (float)(z+cellsPerMeter/2)*step;

					//Project point onto image plane
						float x = point.x()*525 / point.z() + 319.5;
						float y = point.y()*525 / point.z() + 239.5;

						if(x < 0 || x >= 640 || y < 0 || y >= 480)
							continue;


					/*int index = ((int)y * 640 + (int)x)*2 + 1;
					int color_index = ( (int)y*640 + (int)x )*4;
*/
					int index = ((int)y * 640 + (int)x);
					int color_index = ( (int)y*640 + (int)x )*3;
					float val = (float)depth[index]*0.001;

					lambda.x() = (x-319.5) / 525;
					lambda.y() = (y-239.5) / 525;
					lambda.z() = 1;
				

					float sd = val - 1.0 / (lambda.norm()) * point.norm();
	
					sd > mu ? mu : sd;
					sd < -mu ? -mu : sd;

					point = R*point + T;
				
					//search for a node at the specified point
					ColorOcTreeNode *node = tsdf->search( octomap::point3d (point.x(),point.y(),point.z()), 0 );
					if(node == NULL){
						if(fabs(sd) < mu){ //no node and signed distance is below truncation
							int red = (int)color[color_index+2];
							int green = (int)color[color_index+1];
							int blue = (int)color[color_index];
						
							octomap::ColorOcTreeNode *new_node = tsdf->updateNode(octomap::point3d (point.x(),point.y(),point.z()), true);
							new_node->setValue(sd);
							new_node->setColor(blue,green,red);
						}
					}else if(sd > 0){ //node exists and is "behind" the surface
						float newVal = (sd + node->getValue())/2; //compute new value for node
						if( fabs(newVal) < mu ){ //update node if it still below truncation, otherwise delete node
							node->setValue(newVal);
						}
						else
							tsdf->deleteNode(octomap::point3d (point.x(),point.y(),point.z()));
					}			
				}
			}
		}
	}

	void SurfaceEstimation(Eigen::Matrix<float, 3, 3, Eigen::RowMajor> &R, Eigen::Vector3f &T){
		octomap::KeyRay keys;

		octomap::point3d end;
		
		//octomap::point3d direction(0,0,1);
		octomap::point3d start(T.x(), T.y(), T.z());
		octomap::point3d coord;
		octomap::point3d zeroCrossing;
		IplImage *img = cvCreateImage(cvSize(640,480), IPL_DEPTH_8U, 1);

		printf("Starting estimation\n");
		for(int row = 0; row < 480; row++){
			for(int col = 0; col < 640; col++){

				//Compute direction of ray passing through pixel (col, row)
				float X = ( (float)col - 319.5) / 525.f;
				float Y = ( (float)row - 239.5) / 525.f;

				octomap::point3d direction(X,Y,1);
				direction.normalize();
			
				if( tsdf->castRay( octomap::point3d(T.x(), T.y(), T.z()), direction, end, true) ){ //ray hit an occupied cell
					tsdf->computeRayKeys(end, end + direction*tsdf->getResolution()*10, keys); //trace nodes near first hit cell
					float prev_val;


					for(std::vector<octomap::OcTreeKey>::iterator key_itr = keys.begin(); key_itr < keys.end(); key_itr++){ //find zero crossing
						octomap::ColorOcTreeNode *node = tsdf->search(*key_itr);
				
						if(node != '\0'){
							octomap::point3d temp = tsdf->keyToCoord(*key_itr);
							//printf("%f\t(%f %f %f)\n", node->getValue(), temp.x(), temp.y(), temp.z());

							if(key_itr == keys.begin()){
								coord = tsdf->keyToCoord(*key_itr);
								prev_val = node->getValue();
								continue;
							}
							else{
								if(prev_val / fabs(prev_val) != node->getValue() / fabs(node->getValue()) ){ //changed sign
									//printf("SIGN CHANGE\n");
									//printf("prev = %f, curr = %f\n", prev_val, node->getValue());
									//printf("prev = (%f %f %f)\tcurr = (%f %f %f)\n", coord.x(), coord.y(), coord.z(), temp.x(), temp.y(), temp.z());

									float distance = fabs(node->getValue() - prev_val);
									float alpha = fabs(prev_val) / distance;
									zeroCrossing = coord + (temp - coord) * alpha;

									//printf("alpha = %f\n", alpha);
									//printf("zero crossing = (%f %f %f)\n", zeroCrossing.x(), zeroCrossing.y(), zeroCrossing.z());

									break;
								}else{ //didnt change sign
									coord = tsdf->keyToCoord(*key_itr);
									prev_val = node->getValue();
								}
							}

					
						}
				
					}

					cvSet2D(img, row, col, cvScalar(zeroCrossing.z()/.0078));
				}else{
					cvSet2D(img, row, col, cvScalar(0));
				}
			}
		}
		printf("done\n");

		cvNamedWindow("surface", CV_WINDOW_AUTOSIZE);
		cvShowImage("surface", img);
		cvWaitKey(0);
		/*printf("%f %f %f\n", end.x(), end.y(), end.z());

		tsdf->computeRayKeys(octomap::point3d(T.x(), T.y(), T.z()), octomap::point3d(T.x(), T.y(), T.z()+10), keys);
		
		octomap::KeyRay::iterator key_itr = keys.begin();
		printf("%d\n", keys.size());
		for(; key_itr != keys.end(); key_itr++){
			
		}*/


	}

	octomap::ColorOcTree &GetTSDF(){ return *tsdf; }

private:
	octomap::ColorOcTree *tsdf;
};

void CreateTSDF(const char *depth_fname, const char *color_fname, octomap::ColorOcTree &tree, Eigen::Matrix<float, 3, 3, Eigen::RowMajor> &R, Eigen::Vector3f &T){
	int size = 512;
	int cellsPerMeter = 128;
	float step = (1.0 / (float)cellsPerMeter) ;

	static unsigned short depth[640*480];
	static unsigned char color[640*480*3];
	Load(depth, color, depth_fname, color_fname);

	std::vector<Eigen::Vector3f> points;
	std::vector<Eigen::Vector3d> colors;

	Eigen::Vector3f point;
	Eigen::Vector3f lambda;
	
	float mu = 0.01;
	for(int row = 0; row < size; row ++){
		for(int col = 0; col < size; col ++){
			for(int z = 0; z < size; z ++){
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
				
				ColorOcTreeNode *node = tree.search( octomap::point3d (point.x(),point.y(),point.z()), 0 );
				if(node == NULL){
					if(fabs(sd) < mu){
						int temp = (int)(sd / (2.0*mu / 256.0)) + 128;
						temp > 255 ? 255 : temp;
						temp < 0 ? 0 : temp;

						int red = (int)color[color_index+2];
						int green = (int)color[color_index+1];
						int blue = (int)color[color_index];
						
						octomap::ColorOcTreeNode *new_node = tree.updateNode(octomap::point3d (point.x(),point.y(),point.z()), true);
						new_node->setValue(sd);
						new_node->setColor(blue,green,red);
					}
				}else if(sd > 0){
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
}



void Load2(unsigned short *depth, unsigned char *color, const char *depth_fname, const char *color_fname){
	char fname[1000];
	sprintf(fname, "depth/%s", depth_fname);
	IplImage *freiburg_depth = cvLoadImage(fname, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR);
	
	char *ptr = freiburg_depth->imageData;
	unsigned short temp[640*480];

	cvNamedWindow("depth", CV_WINDOW_AUTOSIZE);
	cvShowImage("depth", freiburg_depth);
	cvWaitKey(25);
	for(int i = 0; i < 640*480; i++, ptr += sizeof(unsigned short)){
		memcpy(&(depth[i]), ptr, sizeof(unsigned short));
		depth[i] /= 5;
	}
	cvRelease( (void **)&freiburg_depth);

	//color
	sprintf(fname, "rgb/%s", color_fname);
	IplImage *freiburg_color = cvLoadImage(fname, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR);
	for(int i = 0; i < 640*480; i++){
		color[i*3] = freiburg_color->imageData[i*3];
		color[i*3+1] = freiburg_color->imageData[i*3+1];
		color[i*3+2] = freiburg_color->imageData[i*3+2];
	}
	cvRelease( (void **)&freiburg_color);
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

	unsigned short *depth = (unsigned short *)malloc(640*480*sizeof(unsigned short));
	unsigned char *color = (unsigned char *)malloc(640*480*3*sizeof(unsigned short));

	TSDFOctree tsdf(0.01);

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
	
		counter++;
		printf("%d.) ", counter);
		printf("%s %s\n", depth_fname, color_fname);
		Load2(depth, color, depth_fname, color_fname);
		cout << R << endl << endl;
		printf("%f %f %f\n", x, y, z);

		R(0,0) = 1;	R(0,1) = 0;	R(0,2) = 0;
		R(1,0) = 0;	R(1,1) = 1;	R(1,2) = 0;
		R(2,0) = 0;	R(2,1) = 0;	R(2,2) = 1;

		tsdf.IntegrateData(depth, color, R, Eigen::Vector3f(0,0,0));
		tsdf.SurfaceEstimation(R, Eigen::Vector3f(0,0,0));

		//CreateTSDF(depth_fname, color_fname, tsdf_tree, R, Eigen::Vector3f(x,y,z));
		if(counter == atoi(argv[1]))
			break;
	}

	//PointsFromOctree(tree, "occCloud.ply");
	//PointsFromOctree(tsdf.GetTSDF(), "tsdfCloud.ply");
	fclose(trajFptr);
}