# Kalman Filter

This is a Java implementation of the discrete discrete Kalman filter or most commonly referred to as the regular Kalman Filter. In case of a linear model, the filter will serve as a regular Kalman filter and in the case of a non linear model the code will serve as an extended Kalman Filter. The Kalman filter is of square root flavour. Square root does not change the math behind the algorithm but ensures that numerical rounding errors that are natural to a computer do not prevent the computed covariance matrices maintaining their property of being positive semi definite.
 
The Java implementation here is standalone meaning that you do not need to import any other extra libraries. All the code is contained here including the basic matrix operations, QR decompositions, basic solvers and so forth. The basic matrix operations were implemented based on `Numerical Linear Algebra` by `Trefethen and Bau`.

The code here was verified against already existing and tested Python implementation in [examples directory](https://github.com/mannyray/KalmanFilter/tree/master/java_implementation/examples) by running `./script.sh`. The output is:

```
difference between Java and Python ddekf implementation:
max relative error for estimates: 4.5871854250104974e-12
max relative error for covariances: 2.0008819506531873e-13
```

The max relative error between the computed estimates and covariances of the Python and Java implementations is very small meaning that the Java code was accurately implemented with respect to the Python code. 


### How to use

The explanation here will be succinct and assumes the reader has a familiarity with the Kalman filter. 

The example `test2.java` in [examples directory](https://github.com/mannyray/KalmanFilter/tree/master/java_implementation/examples) is a good example to see how to implement code for your own context. Just copy the `java` files in `java_implementation` directory to your project.

You first need a model which is defined by class `Function` (or its subclass in the case of `test2.java`). In this class you need to define a `next` method which computes `x_{k+1} = next(x_k,k)` and another method that computes the Jacobian. You can use `Function`'s constructor to define constants/matrices that might be needed for your computations. All Matrix based operations can be found in the custom `Matrix.java`.

Your model requires to define the process noise matrix Q and its sensor noise matrix R. However, as seen in the example you need to provide the root of the matrices (Cholesky decomposition transposed) to the `ddekf` function. Unfortunately, Cholesky decomposition was not included in the code here so in `test2.java` these root matrices are entered manually. The Cholesky should not be too hard to implement if really needed.

In addition, you need to provide observation matrix C, the initial estimate and initial covariance estimate (its root). Finally, you need to provide the array of measurements. The output will be a two dimensional array where the first row is the estimates with column n of array representing the nth timestep's estimate while the second row are the covariances.

To run the Kalman filter in a realtime setting see the implementation of `ddekf` and extract the calls of `updatePhase` and `predictPhase` to your relevant context.
