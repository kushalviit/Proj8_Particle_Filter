/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>


#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

   default_random_engine gen;
        normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);
        num_particles=300;
        for(int i=0;i< num_particles;i++)
        {
          Particle temp_p;
          temp_p.id=i;
          temp_p.x=dist_x(gen);
          temp_p.y=dist_y(gen);
          temp_p.theta=dist_psi(gen);
          temp_p.weight=1.0;
          particles.push_back(temp_p);
          weights.push_back(1.0);
         }
   is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  // define normal distributions for sensor noise
  default_random_engine gen;
         normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_psi(0, std_pos[2]);
         
         for(int i=0;i< num_particles;i++)
         {
          if(fabs(yaw_rate) > 0.00001)
           {
           
           particles[i].x+=((velocity/yaw_rate)*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta)))+dist_x(gen);
           particles[i].y+=((velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t)))+dist_y(gen);
           particles[i].theta=particles[i].theta+yaw_rate*delta_t+dist_psi(gen);
           }
          else
           {
             particles[i].x+=((velocity*delta_t)*(cos(particles[i].theta)))+dist_x(gen);
             particles[i].y+=((velocity*delta_t)*(sin(particles[i].theta)))+dist_y(gen);
             particles[i].theta+=dist_psi(gen);
           }
         }
  
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.
  for(int i=0;i<observations.size();i++)
              { 
                  double observations_x=observations[i].x;
                  double observations_y=observations[i].y; 
                  double min_dis=dist(observations_x,observations_y, predicted[0].x, predicted[0].y);
                   observations[i].id=predicted[0].id;
                   for(int j=1;j<predicted.size();j++)
                       { 
                          double dis=dist(observations_x, observations_y, predicted[j].x, predicted[j].y);
                          if(dis<min_dis)
                             {
                               observations[i].id=predicted[j].id;
                               min_dis=dis;
                              }
         
                        }  
               
              }
          
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
    std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation 
  //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
  //   for the fact that the map's y-axis actually points downwards.)
  //   http://planning.cs.uiuc.edu/node99.html

  // for each particle...
   double prob_sum=0.0;
         weights.clear();
         for(int i=0;i<num_particles;i++)
         {
         //transformation
             std::vector<LandmarkObs> mp_ldmarks;
            std::vector<LandmarkObs> map_coor_observations=observations;
              for(int j=0;j<observations.size();j++)
              {
                 map_coor_observations[j].x=observations[j].x*cos(particles[i].theta)-observations[j].y*sin(particles[i].theta)+particles[i].x;
                 map_coor_observations[j].y=observations[j].x*sin(particles[i].theta)+observations[j].y*cos(particles[i].theta)+particles[i].y;
              }
             
             
            for(int j=0;j<map_landmarks.landmark_list.size();j++)
             {    
                  LandmarkObs temp;
                  temp.id=map_landmarks.landmark_list[j].id_i;
                  temp.x=map_landmarks.landmark_list[j].x_f;
                  temp.y=map_landmarks.landmark_list[j].y_f;
                  //double dis=dist(particles[i].x,particles[i].y,temp.x,temp.y);
                if(fabs(particles[i].x - temp.x) <= sensor_range && fabs(particles[i].y - temp.y) <= sensor_range)
                  mp_ldmarks.push_back(temp);
                  
              }
             
             
         //Association
         
         dataAssociation(mp_ldmarks,map_coor_observations);
             std::vector<int> associations;
             std::vector<double> sense_x;
             std::vector<double> sense_y;
              for(int j=0;j<map_coor_observations.size();j++)
                {
                  associations.push_back(map_coor_observations[j].id);
                  sense_x.push_back(map_coor_observations[j].x);
                  sense_y.push_back(map_coor_observations[j].y);
                 }

          particles[i]=SetAssociations(particles[i], associations, sense_x, sense_y);
        
         //assigning weights for each particle
          
               double temp_prob=1.0;
               for(int j=0;j<particles[i].associations.size();j++)
               {
                  double s_x=particles[i].sense_x[j];
                  double s_y=particles[i].sense_y[j];
                  int pos;
                     for(int k=0;k<mp_ldmarks.size();k++)
                      {
                          if(particles[i].associations[j]==mp_ldmarks[k].id)
                             {
                             pos=k;
                             break;
                             }
                       }
                   double m_x=mp_ldmarks[pos].x;
                   double m_y=mp_ldmarks[pos].y;
                   double var_2x=2*std_landmark[0]*std_landmark[0];
                   double var_2y=2*std_landmark[1]*std_landmark[1];
                   double parm_x=(s_x-m_x)*(s_x-m_x);
                   double parm_y=(s_y-m_y)*(s_y-m_y);
                   double fac=1/(2*M_PI*std_landmark[0]*std_landmark[1]);
                   temp_prob*=(fac*exp(-1*(parm_x+parm_y)));
               }
               
              weights.push_back(temp_prob);
              prob_sum+=temp_prob;
             }
          for(int i=0;i<num_particles;i++)
           { 
            
            particles[i].weight=weights[i]/prob_sum;
             weights[i]=particles[i].weight;
            }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
      
  std::default_random_engine generator;
      
      std::uniform_int_distribution<int> index_dist(0,num_particles-1);
      std::vector<Particle> temp_particles;
      double beta=0.0;
      int index=index_dist(generator);
      double max_weight=*std::max_element(weights.begin(),weights.end());
      std::uniform_real_distribution<double> weight_dist(0,max_weight);
      for(int i=0;i<num_particles;i++)
      {
         double random=weight_dist(generator);
         beta += random * 2.0 * max_weight;
          while(beta>weights[index])
          {
           beta -= weights[index];
            index = (index + 1) % num_particles;
           }
          temp_particles.push_back(particles[index]);
      }
      
      particles=temp_particles;
     for(int i=0;i<num_particles;i++)
      { 
        weights[i]=particles[i].weight;
      }
}


Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
