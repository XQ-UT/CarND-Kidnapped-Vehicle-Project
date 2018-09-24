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
#include <limits>
#include <unordered_map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	
	// Set the number of particles 
	num_particles = 10;

	for(int i = 0; i < num_particles; ++i){
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;

		particles.push_back(p);
		weights.push_back(1.0);
	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];
	default_random_engine gen;

	for(Particle& p : particles){
		if( yaw_rate == 0.0){
			// deal the case when yaw_rate is zero.
			p.x += velocity * delta_t * cos(p.theta);
			p.y += velocity * delta_t * sin(p.theta);
		}else{
			p.x += velocity * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta)) / yaw_rate;
			p.y += velocity * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t)) / yaw_rate;
			p.theta += yaw_rate * delta_t;
		}

		normal_distribution<double> dist_x(p.x, std_x);
		normal_distribution<double> dist_y(p.y, std_y);
		normal_distribution<double> dist_theta(p.theta, std_theta);
		
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for(LandmarkObs& ob : observations){
		double min_distance = numeric_limits<double>::max();
		for(const LandmarkObs& predicted_ob: predicted){
			double distance = dist(predicted_ob.x, predicted_ob.y, ob.x, ob.y);
			if(distance <= min_distance){
				min_distance = distance;
				ob.id = predicted_ob.id;
			}	
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	for(int i = 0; i < num_particles; ++i){
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		// convert observation from vehicle coordicates to map coordinates system.
		std::vector<LandmarkObs> converted_observations;
		for(const LandmarkObs& ob : observations){
			LandmarkObs converted_ob;
			converted_ob.id = ob.id;
			converted_ob.x = x + ob.x * cos(theta) - ob.y * sin(theta);
			converted_ob.y = y + ob.x * sin(theta) + ob.y * cos(theta);
			converted_observations.push_back(std::move(converted_ob));
		}

		// fetch in-range landmarks from map.
		std::vector<LandmarkObs> predicted_landmarks;
		std::unordered_map<int, LandmarkObs> id_to_landmark;
		for(const Map::single_landmark_s& landmark: map_landmarks.landmark_list){
			if(dist(x, y, landmark.x_f, landmark.y_f) >= sensor_range){
				continue;
			}

			LandmarkObs inrange_landmark;
			inrange_landmark.id = landmark.id_i;
			inrange_landmark.x = landmark.x_f;
			inrange_landmark.y = landmark.y_f;
			predicted_landmarks.push_back(inrange_landmark);
			id_to_landmark[inrange_landmark.id] = inrange_landmark;
		}

		dataAssociation(predicted_landmarks, converted_observations);
		double new_weight = 1.0;
		for(const LandmarkObs& c_ob :converted_observations){
			LandmarkObs& p_ob = id_to_landmark[c_ob.id];
			double index = 	pow(p_ob.x - c_ob.x, 2) / (2 * std_x * std_x) +
							pow(p_ob.y - c_ob.y, 2) / (2 * std_y * std_y);
			double prob = exp(-index) / (2 * M_PI * std_x * std_y);
			new_weight *= prob;
		}

		particles[i].weight = new_weight;
		weights[i] = new_weight;
	}		
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	double max_weight = *std::max_element(weights.begin(), weights.end());
	std::discrete_distribution<int> d_index(0, num_particles);
	std::uniform_real_distribution<double> d_prob(0.0, 2 * max_weight);
	default_random_engine gen;
	int index = d_index(gen);
	double beta = 0.0;
	std::vector<Particle> new_particles;
	std::vector<double> new_weights;

	for(int i = 0; i < num_particles; ++i){
		beta += d_prob(gen);
		while(weights[index] < beta){
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		new_particles.push_back(particles[index]);
		new_weights.push_back(particles[index].weight);
	}
	
	particles = new_particles;
	weights = new_weights;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
