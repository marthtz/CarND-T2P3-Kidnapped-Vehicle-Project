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
#include <map>

#include "particle_filter.h"

using namespace std;


#define DEBUG_PRINTS 0



/**
 * init Initializes particle filter by initializing particles to Gaussian
 *   distribution around first position and all the weights to 1.
 * @param x Initial x position [m] (simulated estimate from GPS)
 * @param y Initial y position [m]
 * @param theta Initial orientation [rad]
 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 10;

  // from L14, 4/5
  default_random_engine gen;

  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i=0; i<num_particles; i++)
  {
    Particle particle;

    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1;

    particles.push_back(particle);
    weights.push_back(1);

#if DEBUG_PRINTS
    cout << "Particle init " << i << endl;
    cout << "\tx, y, theta: " << x << " " << y << " " << theta << endl;
    cout << "\tparticle x, y, theta: " << particle.x << " " << particle.y << " " << particle.theta << endl;
#endif
  }

  is_initialized = true;
}




/**
 * prediction Predicts the state for the next time step
 *   using the process model.
 * @param delta_t Time between time step t and t+1 in measurements [s]
 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 * @param velocity Velocity of car from t to t+1 [m/s]
 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;
    double x_pred;
    double y_pred;
    double theta_pred;

  // L12 3 and L14 6/7/8 (org from fusion)
  for (int i=0; i<num_particles; i++)
  {

#if DEBUG_PRINTS
    cout << "Particle predict " << i << endl;
    cout << "\tVelocity:"<<velocity<<";Yawrate:"<<yaw_rate<<std::endl;
    cout << "\tparticle before x, y, theta: " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
#endif

    if (yaw_rate == 0)
    {
      // xf = x0 + v * dt * cos(theta)
      // yf = y0 + v * dt * sin(theta)
      x_pred = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      y_pred = particles[i].y + velocity * delta_t * sin(particles[i].theta);
      theta_pred = particles[i].theta;
    }
    else
    {
      // xf = x0 + v/yaw_rate * (sin(yaw + yaw_rate * dt) - sin(jaw))
      // yf = y0 + v/yaw_rate * (cos(yaw) - cos(yaw + yaw_rate * dt))
      x_pred = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      y_pred = particles[i].y + velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      theta_pred = particles[i].theta + yaw_rate * delta_t;
    }

    // Add noise
    normal_distribution<double> dist_x(x_pred, std_pos[0]);
    normal_distribution<double> dist_y(y_pred, std_pos[1]);
    normal_distribution<double> dist_theta(theta_pred, std_pos[2]);

    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);

#if DEBUG_PRINTS
    cout << "\tx, y, theta: " << x_pred << " " << y_pred << " " << theta_pred << endl;
    cout << "\tparticle after x, y, theta: " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
#endif
  }

}




/**
 * dataAssociation Finds which observations correspond to which landmarks (likely by using
 *   a nearest-neighbors data association).
 * @param predicted Vector of predicted landmark observations
 * @param observations Vector of landmark observations
 */
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase. }


  // predicted = map
  // task: compute nearest neighbour for each obs (using dist) and map
  //       corresponding ID from map to obs
}





/**
 * updateWeights Updates the weights for each particle based on the likelihood of the 
 *   observed measurements. 
 * @param sensor_range Range [m] of sensor
 * @param std_landmark[] Array of dimension 2 [standard deviation of range [m], standard deviation of bearing [rad]]
 * @param observations Vector of landmark observations
 * @param map Map class containing map landmarks
 */
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
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html

  // Loop counter
  int i, j, k;
  double theta;
  //vector<LandmarkObs> obs_glob_vec;
  LandmarkObs obs_glob;

  int num_obs = observations.size();
  int num_lm = map_landmarks.landmark_list.size();

  double ox = std_landmark[0];
  double oy = std_landmark[1];
  double sin_theta;
  double cos_theta;

  // for each particle
  for (i=0; i<num_particles; i++)
  {
    double final_weight = 1;
    double weight;

#if DEBUG_PRINTS
    cout << "Particle update " << i << endl;
    cout << "\tParticle (x,y,theta): " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << endl;
#endif

    // Steps from L14 S13
    sin_theta = sin(particles[i].theta);
    cos_theta = cos(particles[i].theta);
    for (j=0; j<num_obs; j++)
    {
      // step 1: transform obs to global map coordinates
      //            obs_glob_x = particel_x + (obs_x * cos(theta) - obs_y * sin(theta))
      //            obs_glob_y = particel_x + (obs_x * sin(theta) + obs_y * cos(theta))
      obs_glob.x = particles[i].x + (observations[j].x * cos_theta - observations[j].y * sin_theta);
      obs_glob.y = particles[i].y + (observations[j].x * sin_theta + observations[j].y * cos_theta);

      // step 2: dataAssociation
      //dataAssociation(map_landmarks.landmark_list, obs_glob_vec);
      // Assign obs to landmark
      double min_distance=999999.9;
      double distance;
      for (k=0; k<num_lm; k++)
      {
        distance = dist(obs_glob.x,
                        obs_glob.y,
                        map_landmarks.landmark_list[k].x_f,
                        map_landmarks.landmark_list[k].y_f);
        if (distance < min_distance)
        {
          min_distance = distance;
          obs_glob.id = map_landmarks.landmark_list[k].id_i-1;
        }
      }

#if DEBUG_PRINTS
      cout << "\tGlob obs " << j << endl;
      cout << "\t\tObs (x,y): " << observations[j].x << " " << observations[j].y << endl;
      cout << "\t\tObs global (x,y): " << obs_glob.x << " " << obs_glob.y << endl;
      cout << "\t\tAssoc. landmark (x, y, id): " << map_landmarks.landmark_list[obs_glob.id].x_f << " " << map_landmarks.landmark_list[obs_glob.id].y_f << " " << obs_glob.id+1 << endl;
#endif

      // step 3: calc each obs weight:
      // P(x,y) = 1/(2*PI*ox*oy) * exp(-1 * ( (obs_glob_x - pred_x)^2/(2*ox^2) + (obs_glob_y - pred_y)^2/(2*oy^2)))
      double factor1;
      double factor2;

      factor1 = (obs_glob.x - map_landmarks.landmark_list[obs_glob.id].x_f);
      factor1 = factor1 * factor1;
      factor1 = factor1 / (2.0*(ox*ox));

      factor2 = (obs_glob.y - map_landmarks.landmark_list[obs_glob.id].y_f);
      factor2 = factor2 * factor2;
      factor2 = factor2 / (2*(oy*oy));

      weight = 1.0/(2.0*M_PI*ox*oy) * exp( -1 * (factor1 + factor2) );
      final_weight *= weight;

#if DEBUG_PRINTS
      cout << "\t\tWeight: " << weight << endl;
#endif
    }

#if DEBUG_PRINTS
    cout << "\tFinal weight: " << final_weight << endl;
#endif

    // step 4: calc final weight by multiplying all obs weights
    particles[i].weight = final_weight;
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  int i;
  std::vector<Particle> particles_new;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());
  std::map<int, int> m;
  for(i=0; i<num_particles; i++)
  {
    particles_new.push_back(particles[d(gen)]);
  }

  for(i=0; i<num_particles; i++)
  {
    particles[i] = particles_new[i];
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
