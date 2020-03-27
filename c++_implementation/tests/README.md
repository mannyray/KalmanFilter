## Tests
To compile tests you need GTest. Run

```
$ mkdir build
$ cd build
$ cmake ..
$ make
$./runTests
$./runTests2
$./runTests3
$./runTests4
```

 - `runTests` tests continuos solvers
 - `runTests2` tests `predict` method of discrete Kalman filter and then the update phase method
 - `runTests3` compare the output of discrete Kalman filter using three different types for template. Results should all be the same.
 - `runTests4` tests 'randomVector' functions
