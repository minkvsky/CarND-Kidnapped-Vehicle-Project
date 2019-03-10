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

#include "helper_functions.h"

// using std::string;
// using std::vector;
using namespace std;

std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 10;  // TODO ajust: Set the number of particles
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++){
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;

    particles.push_back(p);
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
  for (int i = 0; i < num_particles; i++) {

    // Add measurements to each particle
    if (fabs(yaw_rate) < 0.00001) {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }
    else {
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }

    // add random Gaussian noise
    // necessary 0?
    normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs> &observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
   for (unsigned int i = 0; i < observations.size(); i++){
     LandmarkObs o = observations[i];
     // base on sensor_range
     double min_dist = 50 * 3;
     int min_id;
     // double min_x;
     // double min_y;
     for (unsigned int j = 0; j < predicted.size(); j++){
       LandmarkObs p = predicted[j];
       double op_dist = dist(o.x, o.y, p.x, p.y);
       // find min
       if (op_dist < min_dist){
         min_dist = op_dist;
         min_id = p.id;
         // min_x = p.x;
         // min_y = p.y;
       }
     }
     // assign the observed measurement to this particular landmark
     observations[i].id = min_id;
     // observations[i].x = min_x;
     // observations[i].y = min_y;
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
   // TODO Update the weights of each particle using a mult-variate Gaussian
   // TODO get predicted
  for (int i = 0; i < num_particles; i++) {

    // get pridictions
    vector<LandmarkObs> predictions;
    // for each map landmark...
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      int lm_id = map_landmarks.landmark_list[j].id_i;

      if (fabs(lm_x - particles[i].x) <= sensor_range && fabs(lm_y - particles[i].y) <= sensor_range) {
        predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
      }
    }
    // transform
    vector<LandmarkObs> trans_obs;
    for (unsigned int j = 0; j < observations.size(); j++) {
      trans_obs[j] = transform_obs(particles[i], observations[j]);
    }
    // assocations
    dataAssociation(predictions, trans_obs);
    particles[i].weight = 1.0;
    // update weights
    for (unsigned int j = 0; j < trans_obs.size(); j++) {

      // placeholders for observation and associated prediction coordinates
      double pr_x, pr_y;
      // get the x,y coordinates of the prediction associated with the current observation
      for (unsigned int k = 0; k < predictions.size(); k++) {
        if (predictions[k].id == trans_obs[j].id) {
          pr_x = predictions[k].x;
          pr_y = predictions[k].y;
        }
      }

      // calculate weight for this observation with multivariate Gaussian
      double obs_w = multiv_prob(std_landmark[0], std_landmark[1], pr_x, pr_y, trans_obs[j].x, trans_obs[j].y);
      // product of weights
      particles[i].weight *= obs_w;
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


   // get all of the current weights
   vector<double> weights;
   for (int i = 0; i < num_particles; i++) {
     weights.push_back(particles[i].weight);
   }
   // get max weight
   double max_weight = *max_element(weights.begin(), weights.end());

   // generate random starting index for resampling wheel
   uniform_int_distribution<int> uniintdist(0, num_particles-1);
   auto index = uniintdist(gen);

   // uniform random distribution [0.0, max_weight)
   uniform_real_distribution<double> unirealdist(0.0, max_weight);

   // beta?
   double beta = 0.0;
   vector<Particle> tmp_particles;
   // spin the resample wheel!
   for (int i = 0; i < num_particles; i++) {
     beta += unirealdist(gen) * 2.0;
     while (beta > weights[index]) {
       beta -= weights[index];
       index = (index + 1) % num_particles;
     }
     tmp_particles.push_back(particles[index]);
   }

   particles = tmp_particles;
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
