clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

logistic;
close all;

particle_count = 30; 
pkg load statistics;
particle1 = mvnrnd(x_0-15,P_0,particle_count/2)';
particle2 = mvnrnd(x_0+15,P_0,particle_count/2)';
particle = [particle1,particle2];
[estimates_particle, covariances_particle, particles] = particle_filter(next_func,jacobian_func,dt,t_start,state_count,sensor_count,...
  outputs,particle_count,particle,C,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);
  
  
h = figure;
hold on;
lim = 975;
plot(0:dt:(outputs-lim-1)*dt,measurements(1:end-lim),'LineWidth',3); 
plot(0:dt:(outputs-lim-1)*dt,process_noise_data(1:end-lim),'LineWidth',3);
plot(0:dt:(outputs-lim-1)*dt,estimates_particle(1:end-lim-1),'LineWidth',3);

for ii=1:(outputs-lim)
	[xUnique, ignore, ixs] = unique(particles{ii});
	counts = zeros(size(xUnique,1),1);
	for ix = 1:size(counts,1);
		counts(ix) = sum(ixs == ix);
	end
	scatter(dt*(ii-1)*ones(length(xUnique),1),xUnique,counts*20,'g','filled');
end
legend('Measurement','Real data','Estimate');
xlabel('Time')
ylabel('Population')
