/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *  Modified on: 09 Aug 2017
 *      Author: Krishtof Korda
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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  
  num_particles = 20;
  
  // Distribution variables for adding random noise
  normal_distribution<double> x_dist(x, std[0]);
  normal_distribution<double> y_dist(y, std[1]);
  normal_distribution<double> theta_dist(theta, std[2]);
  
  // Declare new particle variable
  Particle p;
  
  for (int i=0; i<num_particles; i++) {
    
    // Assign values based on GPS with random Gaussian noise
    p.id = i;
    p.x = x_dist(gen);
    p.y = y_dist(gen);
    p.theta = theta_dist(gen);
    p.weight = 1.0;
    
    // Store new particle in the vector of particles
    particles.push_back(p);
  }
  
  // Set initialization as complete
  is_initialized = true;
  
  cout << "||||||||number of particles after initialization = " << particles.size() << "||||||||\n\n";
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  
  cout << "||||||||Predicting particle positions after one time step||||||||\n\n";
  
  // declare readability variables
  double x_0;
  double y_0;
  double theta_0;
  double x_f;
  double y_f;
  double theta_f;
  
  for (int i=0; i<particles.size(); i++) {
    
    // extract variables for readability
    x_0 = particles[i].x; //cout << "x = " << x_0 << "\n\n";
    y_0 = particles[i].y; //cout << "y = " << y_0 << "\n\n";
    theta_0 = particles[i].theta; //cout << "theta = " << theta_0 << "\n\n";
    
    
    //// Prediction equations ////
    // Motion prediction equations for yaw_rate ~zero
    if (fabs(yaw_rate) < .001) {
      
      x_f =  x_0 + velocity * cos(theta_0) * delta_t;
      y_f = y_0 + velocity * sin(theta_0) * delta_t;
      theta_f = theta_0;
    }
    
    // Motion prediction equations for non-zero yaw_rate
    else {
      
      x_f =  x_0 + velocity * (sin(theta_0 + yaw_rate * delta_t) - sin(theta_0)) / yaw_rate;
      y_f = y_0 + velocity * (-cos(theta_0 + yaw_rate * delta_t) + cos(theta_0)) / yaw_rate;
      theta_f = theta_0 + yaw_rate * delta_t;
    }
    
    // Noise generators with mean of current GPS position
    normal_distribution<double> x_f_dist(x_f, std_pos[0]);
    normal_distribution<double> y_f_dist(y_f, std_pos[1]);
    normal_distribution<double> theta_f_dist(theta_f, std_pos[2]);
    
    // Assigning predicted noisy position and orientation
    particles[i].x = x_f_dist(gen);
    particles[i].y = y_f_dist(gen);
    particles[i].theta = theta_f_dist(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  
  // Declare index for minimum distance
  int min_pos;
  
  // Declare vectors to store distances
  std::vector<double> distances;
  std::vector<int> ids;
  double dist_i;
  
  // Loop through observations
  for (int i=0; i < observations.size(); i++) {
    
    // Loop through predicted landmark and calculate distance between predicted and observation
    for (int j=0; j < predicted.size(); j++) {
      
      // Calculate distance between current observation  and landmark
      dist_i = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
      
      // Store distance of jth landmark to ith observation
      distances.push_back(dist_i);
    }
    
    // Find the index of the min distance for current observation and landmark
    min_pos = min_element(distances.begin(), distances.end()) - distances.begin();
    
    // Assign closest landmark id to observation
    observations[i].id = predicted[min_pos].id;
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
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  
  cout << "||||||||Updating particle weights||||||||\n\n";
  
  // Pull standard devation values
  const double std_x = std_landmark[0];
  const double std_y = std_landmark[1];
  
  // Declare variable for particle position and yaw
  double px;
  double py;
  double pth;
  
  // Declare variable for landmark position, id, and distance
  double x_f;
  double y_f;
  double id_f;
  double dist_j;
  
  // Declare landmark object variable for transformed observation
  LandmarkObs o_map;
  
  // Declare landmark object variable for current predicted landmark
  LandmarkObs pred_lm;
  
  // Declare variables for Multi-variate probability calculations
  double xdiff;
  double ydiff;
  double P;
  
  // Clear weights variable that is about to be filled
  weights.clear();
  
  // Initialize weight sum variable for normalization
  double weight_sum = 0.0;
  
  // Loop through each particle
  for (int i=0; i<particles.size(); i++) {
    
    // Pull position and yaw info
    px = particles[i].x;
    py = particles[i].y;
    pth = particles[i].theta;

    // Declare Vector of predicted landmarks
    std::vector<LandmarkObs> predicted_landmarks;
    
    // Loop through landmarks to find landmarks in range of each particle.
    for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
      
      // Pull landmark position and id from the map
      x_f = map_landmarks.landmark_list[j].x_f;
      y_f = map_landmarks.landmark_list[j].y_f;
      id_f = map_landmarks.landmark_list[j].id_i;
      
      // Calculate distance of the landmark to observation
      dist_j = dist(px,  py,  x_f,  y_f);
      
      // If landmark is in range append to vector of predicted landmarks
      if (fabs(dist_j) <= sensor_range) {
        
        // Declare landmark object to store in vector of in range landmarks
        LandmarkObs in_range;
        
        // Id and position of matching landmark in the map.
        in_range.id = id_f;
        in_range.x = x_f;
        in_range.y = y_f;
        
        //Store in range landmark
        predicted_landmarks.push_back(in_range);
      }
    }
    
    ///// Transform observations into map coordinates////
    // Declare vector and variable containers for tranformations
    std::vector<LandmarkObs> observations_map_coord;
    LandmarkObs obs_map_coord;
    
    // Loop through observations to be transformed
    for (int k=0; k < observations.size(); k++) {
      
      // Equations for 2D transformation rotation and translation
      obs_map_coord.x = px + observations[k].x * cos(pth) - observations[k].y * sin(pth);
      obs_map_coord.y = py + observations[k].x * sin(pth) + observations[k].y * cos(pth);
      
      // Copy ids from observations to transformed
      obs_map_coord.id = observations[k].id;
      
      // Build the vector of transformed observations
      observations_map_coord.push_back(obs_map_coord);
    }
    
    // Add association landmarks to each particle
    dataAssociation(predicted_landmarks, observations_map_coord);
    
    // Re-initialize particle weight so as not to carry previous probabilities from resample
    particles[i].weight = 1.0;
    
    // Declare constants for Multi-variate Gaussian prob
    const double sig_y2 = 2*std_y*std_y;
    const double sig_x2 = 2*std_x*std_x;
    const double pi_sig2 = 2*M_PI*std_x*std_y;
    
    //// Loop through observations to do Mulit-variate Gaussian probability ////
    for (int k=0; k < observations_map_coord.size(); k++) {
      
      // Pull current observation, k
      o_map = observations_map_coord[k];

      // Loop through landmarks found to be in range
      for (int j=0; j < predicted_landmarks.size(); j++) {
        
        // Pull current landmark, j
        pred_lm = predicted_landmarks[j];
        
        // Check if observation is matched to current landmark
        if (o_map.id == pred_lm.id) {
          
          // Calculate Multi-variate Gaussian probability of current observation being correct
          xdiff = pow((o_map.x - pred_lm.x), 2) / (sig_x2);
          ydiff = pow((o_map.y - pred_lm.y), 2) / (sig_y2);
          P = exp(-(xdiff + ydiff)) / (pi_sig2);
       
          // Multiply current probability with new probability and assign to particle weight
          particles[i].weight *= P;
        }
      }
    }
    
    // Build weights vector for use in resampling
    weights.push_back(particles[i].weight);
    
    // Accumulate weight sum for normalization
    weight_sum += particles[i].weight;
  }
  
  //Normalize weights to a sum of 1
  for (int n=0; n<weights.size(); n++) {
    weights[n] /= weight_sum;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  cout << "||||||||Resampling particles based on weight||||||||\n\n";
  
  //// Determine the max weight of all particles for use in gaussian noise ////
  // Declare and initialize max weight
  double max_weight = *max_element(weights.begin(), weights.end());

  // Normal integer distribution for index of particles
  std::uniform_int_distribution<unsigned int> index_dist(0, num_particles-1);
  unsigned int index = index_dist(gen);

  // Delcare storage variable for newly sampled particles
  std::vector<Particle> resampled_particles;
  
  // Declare and initialize normal distribution of weights
  double beta=0.0;
  double mw2 = 2*max_weight;
  std::uniform_real_distribution<double> weight_dist(0, mw2);
  
  // Resample particles with higher probability of selection based on particle weight
  // This is a resampling wheel implementation learned in lesson
  for (int r=0; r < num_particles; r++) {
    
    beta += weight_dist(gen);

    while ( weights[index] < beta ) {
      
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    
    // Add selected sample to new particles
    resampled_particles.push_back(particles[index]);
  }
  
  // Assign new set of particles to particles vector
  particles = resampled_particles;
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
