/******************************************************************************** 
** @Copyright(c) $year$ $registered organization$ All Rights Reserved.
** @auth： taify
** @date： 2021/01/13
** @desc： MyICP调用demo
** @Ver : V1.0.0
*********************************************************************************/

#include "omp.h"
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h> 
#include <pcl/registration/icp.h>
#include <pcl/visualization/pcl_visualizer.h>
#include "preprocess.h"
#include <sensor_msgs/PointCloud2.h>
#define SKEW_SYM_MATRX(v) 0.0,-v[2],v[1],v[2],0.0,-v[0],-v[1],v[0],0.0

std::mutex mtx_buffer;
std::shared_ptr<Preprocess> p_pre(new Preprocess());
std::deque<PointCloudXYZI::Ptr>        lidar_buffer;
Eigen::Matrix4d transformation_matrix = Eigen::Matrix4d::Identity();
boost::shared_ptr<PointCloudXYZI>  source_cloud(new PointCloudXYZI());
boost::shared_ptr<PointCloudXYZI>  target_cloud(new PointCloudXYZI());
boost::shared_ptr<PointCloudXYZI>  source_cloud_mid(new PointCloudXYZI());
boost::shared_ptr<PointCloudXYZI>  target_cloud_mid(new PointCloudXYZI());
pcl::search::KdTree<PointType>::Ptr kdtree(new pcl::search::KdTree<PointType>);

Eigen::Matrix<double,4,Eigen::Dynamic> Points;
std::vector<Eigen::Matrix<double,4,Eigen::Dynamic>,Eigen::aligned_allocator<Eigen::Matrix<double,4,Eigen::Dynamic>>> result;
bool first= true;
int idx =0;
void icp_(const boost::shared_ptr<PointCloudXYZI>  &source_cloud,const boost::shared_ptr<PointCloudXYZI>  &target_cloud);
void icp_process(const boost::shared_ptr<PointCloudXYZI>  &source_cloud,const boost::shared_ptr<PointCloudXYZI>  &target_cloud);
template <typename S>
void FitLine(const Eigen::Matrix<S,4,Eigen::Dynamic> &re);

void standard_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg)
{
    std::unique_lock<std::mutex> ul(mtx_buffer, std::defer_lock);

	
    PointCloudXYZI::Ptr  ptr(new PointCloudXYZI());
    p_pre->process(msg, ptr);
    ul.lock();
    lidar_buffer.push_back(ptr);
	std::cout <<"lidar_buffer.push_back(ptr);"<< lidar_buffer.size() << "\n";
		// if (lidar_buffer.size() >= 359)
		if (first){
			//do nothing
			auto size =  ptr->points.size();

			Points.resize(4,size/5 +1);
			Points.setZero();

			for (size_t i = 0; i < size; i++)
			{
				if (i % 5 != 0) continue;
				auto &point = ptr->points[i];
				Eigen::Vector4d tmp{point.x,point.y,point.z,1};
				Points.col(idx) = tmp;

				idx++;
			}
			idx = 0;
			first =false;
		}else{
			// for (std::deque<PointCloudXYZI::Ptr> ::iterator it = lidar_buffer.begin(); it!=lidar_buffer.end() -1; ++it)
			std::deque<PointCloudXYZI::Ptr> ::iterator it =  lidar_buffer.begin() + idx;
			{
				target_cloud = *(it);
				source_cloud =  *(it+1);
				icp_(source_cloud,target_cloud);
			}
			idx++;
		}
    ul.unlock();
	if (lidar_buffer.size() ==359)
	{
		Eigen::Matrix<double,4,Eigen::Dynamic> re;
		re.resize(4,Points.cols());
		re.setZero();
		for (auto it : result)
		{
			re += it;
		}
		re = re/(result.size());
		FitLine(re);
	}
	
}

template <typename S>
void FitLine(const Eigen::Matrix<S,4,Eigen::Dynamic> &re)
{
    int num = re.cols();
	if (num < 2) {
		std::cout <<"num less than 2"  << "\n";
		
    }
	Eigen::Matrix<S, 4, 1> dir;
    Eigen::Matrix<S,4,1> points;
    Eigen::Matrix<S,4,4> point_cov;
    Eigen::Matrix<S,4,1> line_center;
	point_cov.setZero();

    for (int j = 0; j < num; j++)
    {

		points = re.col(j); //world point 

        point_cov += points * points.transpose();        
        line_center += points;
    }


	line_center = line_center / num;
    point_cov = point_cov / num - line_center * line_center.transpose();
	std::cout <<"line_center ="<< line_center << "\n";
	Eigen::Matrix<S,3,1> C {line_center(0),line_center(1),line_center(2)};
    Eigen::JacobiSVD<Eigen::Matrix<S,4,4>> svd(point_cov, Eigen::ComputeFullV);
	Eigen::Matrix<S,4,1> Sigma = svd.singularValues();
	std::cout <<"Sigma ="<< Sigma << "\n";

    if (std::isnan(Sigma(0)))
    {
    }
    if (Sigma(0) > 3*Sigma(1))
    {
        dir = svd.matrixV().col(0);
		std::cout <<"dir = "<< dir << "\n";
		Eigen::Matrix<S,3,1> di{dir(0),dir(1),dir(2)};
		// Eigen::Matrix<S,3,3> crossmat;
		// crossmat << SKEW_SYM_MATRX(di);
		for (int i = 0; i < 360; i++)
		{
			std::cout <<"theta  = "<< i  << "\n";
			
			Eigen::AngleAxis<S> rotation_vector(i,di);
			Eigen::Matrix<S,3,3> rotation_matrix = rotation_vector.matrix();
			std::cout <<"rotation_matrix  = "<< rotation_matrix  << "\n";

			Eigen::Matrix<S,3,1> T ;
			T = (Eigen::Matrix<S,3,3>::Identity() - rotation_matrix) * C;
			std::cout <<"T  = "<< T  << "\n";
		}
    }

    
    // for (int j = 0; j < NUM_MATCH_POINTS; j++)
    // {
    //     points(0) = point[j].x; //world point 
    //     points(1) = point[j].y; //world point
    //     points(2) = point[j].z; //world point
    //     if (dir.template cross(points - line_center).template squaredNorm() > eps)
    //     {
    //         return false;
    //     }
    // }
	// std::cout << " Sigma line = " << Sigma.transpose()<<"\n";

}


void icp_process(const boost::shared_ptr<PointCloudXYZI>  &source_cloud,const boost::shared_ptr<PointCloudXYZI>  &target_cloud){

	// boost::shared_ptr<PointCloudXYZI> cloud_icp (new PointCloudXYZI());
	pcl::IterativeClosestPoint<PointType, PointType> icp;
	icp.setMaximumIterations (1);
	icp.setInputSource(source_cloud);
		std::cout << "\n  	source_cloud.size() is " << (*source_cloud).size()  << std::endl;
		std::cout << "\n  	target_cloud.size() is " << (*target_cloud).size() << std::endl;


	icp.setInputTarget(target_cloud);
	// icp.align(*cloud_icp);

	icp.setMaximumIterations (500);
	// icp.setTransformationEpsilon (1e-9);
	// icp.setMaxCorrespondenceDistance (0.05);
	// icp.setEuclideanFitnessEpsilon (1);
	// icp.setRANSACOutlierRejectionThreshold (1.5);


	if (icp.hasConverged ())
	{
		std::cout << "\nICP has converged, score is " << icp.getFitnessScore () << std::endl;
		// std::cout << "\nICP transformation " << iterations << " : cloud_icp -> cloud_in" << std::endl;
		transformation_matrix = icp.getFinalTransformation ().cast<double>();
		std::cout << "\nICP transformation_matrix =" << transformation_matrix << " " << std::endl;

	}
	else
	{
		PCL_ERROR ("\nICP has not converged.\n");
	}

}

void icp_(const boost::shared_ptr<PointCloudXYZI>  &source_cloud,const boost::shared_ptr<PointCloudXYZI>  &target_cloud)
{
	//建立kd树
	kdtree->setInputCloud(target_cloud);

	//初始化变换矩阵等参数
	int iters = 0;
	double error = std::numeric_limits<double>::infinity();
	Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
	Eigen::Vector3d T = Eigen::Vector3d::Zero();
	Eigen::Matrix4d H = Eigen::Matrix4d::Identity();
	Eigen::Matrix4d H_final = H;
	Eigen::Matrix4d transformation_matrix = Eigen::Matrix4d::Identity();
	pcl::copyPointCloud(*source_cloud, *source_cloud_mid);

	//开始迭代，直到满足条件
	clock_t start = clock();
	#pragma omp parallel for
	for (;error > 0.0001 && iters < 100;iters++)
	{
		
		double last_error = error;
		double err = 0.0;
		pcl::transformPointCloud(*source_cloud_mid, *source_cloud_mid, H);
		std::vector<int>indexs(source_cloud_mid->size());
		#pragma omp parallel for
		for (int i = 0; i < source_cloud_mid->size(); ++i)
		{
			std::vector<int>index(1);
			std::vector<float>distance(1);
			kdtree->nearestKSearch(source_cloud_mid->points[i], 1, index, distance);
			err = err + sqrt(distance[0]);
			indexs[i] = index[0];
		}
		pcl::copyPointCloud(*target_cloud, indexs, *target_cloud_mid);
		error = err / source_cloud->size();
		// std::cout << "iters:" << iters << "  " << "error:" << error << std::endl;

		if (fabs(last_error - error) < 1e-6)
			break;
		//计算点云中心坐标
		Eigen::Vector4d source_centroid, target_centroid_mid;
		pcl::compute3DCentroid(*source_cloud_mid, source_centroid);
		pcl::compute3DCentroid(*target_cloud_mid, target_centroid_mid);

		//去中心化
		Eigen::MatrixXd souce_cloud_demean, target_cloud_demean;
		pcl::demeanPointCloud(*source_cloud_mid, source_centroid, souce_cloud_demean);
		pcl::demeanPointCloud(*target_cloud_mid, target_centroid_mid, target_cloud_demean);

		//计算W=q1*q2^T
		Eigen::MatrixXd W = (souce_cloud_demean * target_cloud_demean.transpose()).topLeftCorner(3, 3);

		//SVD分解得到新的旋转矩阵和平移矩阵
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(W, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::MatrixXd U = svd.matrixU();
		Eigen::MatrixXd V = svd.matrixV();

		if (U.determinant() * V.determinant() < 0)
		{
			for (int x = 0; x < 3; ++x)
				V(x, 2) *= -1;
		}

		R = V * U.transpose();
		T = target_centroid_mid.head(3) - R * source_centroid.head(3);
		H << R, T, 0, 0, 0, 1;
		H_final = H * H_final; //更新变换矩阵	
	}
	transformation_matrix << H_final;
	Eigen::Matrix<double,4,Eigen::Dynamic> R_P = transformation_matrix * Points;
	Points = R_P;
	result.push_back(R_P);
	std::cout << "transformation_matrix:" << transformation_matrix   << std::endl;

}





int main(int argc, char** argv)
{
	ros::init(argc, argv, "myicp");
    ros::NodeHandle nh;
	ros::Subscriber sub_pcl =  nh.subscribe("/os_cloud_node/points", 1000, &standard_pcl_cbk );
	omp_set_num_threads(8);

	ros::spin();

	return 0;
}