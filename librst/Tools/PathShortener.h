/**
 * @file Smoothie.h
 */

#include <iostream>
#include <Tools/Robot.h>
#include <ctime>
#include <list>
#include <vector>


#ifndef _SMOOTHIE_
#define _SMOOTHIE_

#define RAND12(N1,N2) N1 + ((N2-N1) * ((double)rand() / ((double)RAND_MAX + 1))) // random # between N&M

/**
 * @function Smoothie
 */
class PathShortener
{
   public:
     EIGEN_MAKE_ALIGNED_OPERATOR_NEW

     PathShortener();
     ~PathShortener();
     void Initialize( World &world,
                      int robotId,
                      std::vector<int> linksId );
     void BruteForce( std::list< Eigen::VectorXd > &rawPath, double stepSize = 0.1);
     bool CheckPathSegment( Eigen::VectorXd config1, Eigen::VectorXd config2 ) const;

     std::list< Eigen::VectorXd > mBruteForcePath; 
   private:
     World* mWorld;
     int mRobotId;     
     std::vector<int> mLinksId;

     Eigen::VectorXd mStartConfig;
     Eigen::VectorXd mGoalConfig;   
     Eigen::Vector3d mStartPos;
     Eigen::Vector3d mGoalPos;   

     double mStepSize;

     std::list< Eigen::VectorXd > mRawPath;
     std::list< Eigen::VectorXd > mSmoothiePath; 
  

     // Heuristic stuff
     Eigen::Vector3d mX1;
     Eigen::Vector3d mX2;
     Eigen::Vector3d mX2_X1;
     float mX2_X1Norm;   

};

#endif /** _SMOOTHIE_ */
