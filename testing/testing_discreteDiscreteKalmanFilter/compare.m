%in seperate terminal go to c++_implementation and run 'make'
%go to bin/test5 and run './test5'
% now comment out 'return' below and run the rest of the script
%return 

estimates_c = csvread('output_data_filtered_data.txt');
estimates_c = estimates_c';
covariances_c = csvread('output_data_covariances.txt');
covariances_cpp = cell(length(covariances)-1,1);
for i=1:length(covariances_cpp)
  covariances_cpp{i} = covariances_c(...
  ((i-1)*state_count+1):((i)*state_count),:);
end
%make sure the relative error is really small between c++ and matlab
relative_tolerance = 0.00001;
for i=1:length(covariances_cpp)
  assert(norm(estimates_c(:,i) - estimates(:,i))./...
  norm(estimates(:,i)) < relative_tolerance);
  assert(   (norm(covariances{i} - covariances_cpp{i})...
  /norm(covariances{i}) ) < relative_tolerance);
  assert( (trace(covariances{i}) - trace(covariances_cpp{i}))./trace(covariances{i}) < relative_tolerance);
end
disp('c++ discrete discrete Kalman filter matches matlab implementation')