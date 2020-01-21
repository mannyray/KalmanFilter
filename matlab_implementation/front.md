# Kalman Filter

This website has various implementations in Matlab(/Octave) of different flavours of the Kalman filter. Check out the tabs of this website to explore the filters. The repository is located at

**[https://github.com/mannyray/KalmanFilter](https://github.com/mannyray/KalmanFilter)**


The insipiration to create this repository is rlabbe's \[1] great introduction to the Kalman filter. This repository aims to provide users a basic and ready to use arsenal to use in exploring filtering. I spent some time working with the Kalman Filter as part of my thesis \[2] where I coded up `CD-EKF` and `DD-EKF`. After coding up the two filters, I decided to keep things interesting and added other filters as well. In order to test my implementations, I used the filters in various contexts as well as checked if the steady state covariances match (see Example 2 of `DD-EKF`, `UKF`, `CD-EKF`, `Particle`, `ENKF-stochastic` and the first example of `CC-KF` and `SS-KF`) as explained in \[3].


If you have further questions or want to contribute, then feel free to reach me at `contact at zonov dot ca`. This repository will grow to be a great collection of various filter implementations for anyone to use and refer to. For contribution, I don't just need help with implementing filters - most importantly I need examples in order to make the code more accessible for newcomers. **If you end up using the code and have a code sample you would like to share then please let me know - I'll feature it here!** If you _do_ have additional filter implementation then I'll gladly make you a contributor to the repository.


\[1] [https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python](https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python)

\[2] [https://uwspace.uwaterloo.ca/bitstream/handle/10012/14740/Zonov_Stanislav.pdf](https://uwspace.uwaterloo.ca/bitstream/handle/10012/14740/Zonov_Stanislav.pdf) (see chapter 3 for a good introduction to the Kalman Filter) 

\[3] Smith, M. W. A., and A. P. Roberts. "An exact equivalence between the discrete-and continuous-time formulations of the Kalman filter." Mathematics and Computers in Simulation 20.2 (1978): 102-109.
