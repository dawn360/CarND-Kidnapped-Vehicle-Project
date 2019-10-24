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
#include <random>
#include <string>
#include <vector>
#include <limits>       // std::numeric_limits


#include "helper_functions.h"

using std::string;
using std::vector;

using std::normal_distribution;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::numeric_limits;

double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;

  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta

  num_particles = 200;  // TODO: Set the number of particles Done
  for (int i = 0; i < num_particles; ++i) {
    Particle np;
    np.id = i;
    np.x = x;
    np.y = y;
    np.theta = theta;
    np.weight = 1;
    particles.push_back(np);
  } 
  // set Standard Deviations
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];
  // create a normal (Gaussian) distribution for x,y and theta
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  for (int i = 0; i < num_particles; ++i) {    
    
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }

  is_initialized = true;

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
  std::default_random_engine gen;
  double x,y,theta,std_x, std_y, std_theta,
        xf,yf,thetaf;  // Standard deviations for x, y, and theta

// set Standard Deviations
  std_x = std_pos[0];
  std_y = std_pos[1];
  std_theta = std_pos[2];

   for (int i = 0; i < num_particles; ++i) { 
     x = particles[i].x;
     y = particles[i].y;
     theta = particles[i].theta;

    if (yaw_rate == 0) {
      xf = x + ((velocity*(delta_t))*cos(theta));
      yf = y + ((velocity*(delta_t))*sin(theta));
      thetaf = theta;

    } else {
      xf = x + (velocity/yaw_rate)*(sin(theta+(yaw_rate*(delta_t))) - sin(theta));
      yf = y + (velocity/yaw_rate)*(cos(theta)-cos(theta+(yaw_rate*(delta_t))));
      thetaf = theta+(yaw_rate*(delta_t));
    }
    
    // create a normal (Gaussian) distribution for x,y and theta
    normal_distribution<double> dist_x(xf, std_x);
    normal_distribution<double> dist_y(yf, std_y);
    normal_distribution<double> dist_theta(thetaf, std_theta);
  
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }


  
  // // create a normal (Gaussian) distribution for x,y and theta
  // normal_distribution<double> dist_x(x, std_x);
  // normal_distribution<double> dist_y(y, std_y);
  // normal_distribution<double> dist_theta(theta, std_theta);
  
  // for (int i = 0; i < num_particles; ++i) {    
    
  //   particles[i].x = dist_x(gen);
  //   particles[i].y = dist_y(gen);
  //   particles[i].theta = dist_theta(gen);
  // }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   // define coordinates and theta
  
  for(auto &obs:observations){
    double shortest = numeric_limits<double>::max();
    for (auto &pred:predicted){
      double distance = dist(obs.x,obs.y,pred.x,pred.y);
      if( distance <= shortest )
      {
        shortest = distance;
        obs.id = pred.id; // assign predicted identifier to obs
      }
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

   // for each particle 
  // 1. based on the sensor_range predict a set of observations from map
  // 2. Find the predicted measurement that is closest to each 
  //    observed measurement and assign the observed measurement to this 
  //    particular landmark.
  // 3.
  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];
  for (unsigned int i = 0; i < particles.size(); ++i) {
    //get predicted based on sensor range
    Particle p = particles[i];
    vector<LandmarkObs> predicted_f;
    for (auto &lmark:map_landmarks.landmark_list) {
      // add landmark if in sensor range
      if (dist(p.x,p.y,lmark.x_f,lmark.y_f) <= sensor_range) {
        predicted_f.push_back(LandmarkObs{lmark.id_i,lmark.x_f,lmark.y_f});
      }
    }

    // create a new vector of observations transformed 
    // from world coordinates to map coordinates

    vector<LandmarkObs> transObservations;
    for (auto &obs:observations) {
      double x_map,y_map;

      // transform to map x coordinate
      x_map = p.x + (cos(p.theta) * obs.x) - (sin(p.theta) * obs.y);

      // transform to map y coordinate
      y_map = p.y + (sin(p.theta) * obs.x) + (cos(p.theta) * obs.y);
      transObservations.push_back(LandmarkObs{obs.id,x_map,y_map});
    }
  // associate closest predicted landmark to observation
  dataAssociation(predicted_f,transObservations);

  // update weights
    particles[i].weight = 1;
    for (auto &obs:transObservations){
      for (auto &mu:predicted_f){
        if (obs.id == mu.id){
          // Calculate OBS1 weight
          particles[i].weight *= multiv_prob(sig_x, sig_y, obs.x, obs.y, mu.x, mu.y);
        }
      }
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
   //Get max weight.
  std::default_random_engine gen;
  double maxWeight = numeric_limits<double>::min();
  for(int i = 0; i < num_particles; ++i) {
    if(particles[i].weight > maxWeight) {
      maxWeight = particles[i].weight;
    }
  }

  uniform_real_distribution<double> distDouble(0.0, maxWeight);
  uniform_int_distribution<int> distInt(0, num_particles-1);
  int index = distInt(gen);
  double beta = 0.0;
  vector<Particle> resampledParticles;
  for(int i = 0; i < num_particles; ++i) {
    beta += distDouble(gen) * 2.0;
    while(beta > particles[index].weight) {
      beta -= particles[index].weight; //weights[index];
      index = (index + 1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }

  particles = resampledParticles;

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

