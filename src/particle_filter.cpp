/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::default_random_engine;
using std::discrete_distribution;
using std::uniform_real_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  this->num_particles = 1000;  // TODO: Set the number of particles
  default_random_engine gen;
  normal_distribution<double> distr_x(x,std[0]);
  normal_distribution<double> distr_y(y,std[1]);
  normal_distribution<double> distr_theta(theta,std[2]);

  for(int i=0;i<num_particles;i++){
    Particle tmp;
    tmp.id = i;
    tmp.x = distr_x(gen);
    tmp.y = distr_y(gen);
    tmp.theta = distr_theta(gen);
    tmp.weight = 1.0;
    particles.push_back(tmp);
  }

  this->is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  default_random_engine gen;
  normal_distribution<double> distr_x(0,std_pos[0]);
  normal_distribution<double> distr_y(0,std_pos[1]);
  normal_distribution<double> distr_theta(0,std_pos[2]);

  std::vector<Particle>::iterator itr = particles.begin();
  for(;itr!=particles.end();++itr){
    Move((*itr),velocity,yaw_rate,delta_t,distr_x,distr_y,distr_theta,gen);
  }
}

void ParticleFilter::Move(Particle& particle, double vel, double yawRate, double dt,
            std::normal_distribution<double>& distr_x,
            std::normal_distribution<double>& distr_y,
            std::normal_distribution<double>& distr_theta,
            std::default_random_engine& gen)
{
  double dx = 0, dy = 0, dtheta = 0;
  if(yawRate == 0){
    dx = vel*dt*cos(particle.theta);
    dy = vel*dt*sin(particle.theta);
  }
  else{
    dx = vel/yawRate*(sin(particle.theta+yawRate*dt)-sin(particle.theta));
    dy = vel/yawRate*(cos(particle.theta)-cos(particle.theta+yawRate*dt));
    dtheta = yawRate*dt;
  }
  
  std::normal_distribution<double>::param_type xParam(dx, distr_x.param().stddev());
  std::normal_distribution<double>::param_type yParam(dy, distr_y.param().stddev());
  std::normal_distribution<double>::param_type thetaParam(dtheta, distr_theta.param().stddev());
  particle.x += distr_x(gen,xParam);
  particle.y += distr_y(gen,yParam);
  particle.theta += distr_theta(gen,thetaParam);
}

void ParticleFilter::dataAssociation(Particle& particle,
                                     const vector<LandmarkObs>& predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  std::vector<int> asso;
  std::vector<double> assoX;
  std::vector<double> assoY;
  //TODO: re-write using KD Tree
  for(auto itr = observations.begin(); itr != observations.end(); ++itr){
    //iterate through the copy of observation, do coord trans into WCS
    double tmpX = itr->x, tmpY = itr->y;
    itr->x = particle.x + tmpX*cos(particle.theta) - tmpY*sin(particle.theta);
    itr->y = particle.y + tmpX*sin(particle.theta) + tmpY*cos(particle.theta);
    double minDist = 100000;
    //compare and associate with list of landmarks in range in WCS
    for(auto & landmark : predicted){
      double distance = dist(itr->x,itr->y,landmark.x,landmark.y);
      if(distance<minDist){
        minDist = distance;
        itr->id = landmark.id;
      }
    }
    //in this particls's copy of observations, coord is now WCS, ids are their matched LM's ids
    //Now set association
    asso.push_back(itr->id);
    assoX.push_back(itr->x);
    assoY.push_back(itr->y);
  }
  SetAssociations(particle,asso,assoX,assoY);

}

void ParticleFilter::Observe(const Particle& particle,
              const Map& map,
              double sensor_range,
              std::vector<LandmarkObs>& LandMarksInRange)
{
  for(auto & landmark : map.landmark_list){
    if(dist(landmark.x_f,landmark.y_f,particle.x,particle.y)<=sensor_range){
      LandmarkObs tmpObs;
      tmpObs.id = landmark.id_i;
      tmpObs.x = landmark.x_f;
      tmpObs.y = landmark.y_f;
      LandMarksInRange.push_back(tmpObs);
    }
  }
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  double normalizer = 0;
  for(auto itr = particles.begin();itr!=particles.end();++itr){
    std::vector<LandmarkObs> landMarksInRange;
    std::vector<LandmarkObs> tmpObservations = observations;
    Observe((*itr),map_landmarks,sensor_range,landMarksInRange);
    dataAssociation((*itr),landMarksInRange,tmpObservations);
    //update weight for this particle
    for(unsigned int i=0;i<itr->associations.size();i++){
      //Get pseudo ranges as map landmark's x and y
      double pseudo_x, pseudo_y;
      for(auto& Lmk_gt : landMarksInRange){
        if(itr->associations.at(i)==Lmk_gt.id){
          pseudo_x = Lmk_gt.x;
          pseudo_y = Lmk_gt.y;
          break;
        }
      }
      itr->weight *= multiv_prob(std_landmark[0],std_landmark[1],
                                 itr->sense_x.at(i),itr->sense_y.at(i),
                                 pseudo_x,pseudo_y);
      //end multiplying weight for one measurement
    }
    //end updating weight for one particle
    normalizer += itr->weight;
  }
  if(normalizer>0){
    for(auto & ptcl : particles){
      ptcl.weight /= normalizer;
    }
  }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  default_random_engine gen;
  discrete_distribution<int> distr_i(0,num_particles);
  std::vector<Particle> tmpParticles;
  //TODO: re-write Particle to compare
  //auto itr_max = std::max_element(particles.begin(),particles.end());
  double max_w = -1;
  double beta = 0;
  for(const auto& particle : particles){
    max_w = (particle.weight>max_w) ? particle.weight : max_w;
  }
  uniform_real_distribution<double> distr_w(0.0,max_w*2.0);
  for(int i=0; i<num_particles; i++){
    unsigned int ind = distr_i(gen);
    beta += distr_w(gen);
    while(beta>particles.at(ind).weight){
      beta -= particles.at(ind).weight;
      ind = (ind+1)%num_particles;
    }
    tmpParticles.push_back(particles.at(ind));
  }
  particles = tmpParticles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}