#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include "model.h"
#include "filterModel.h"
#include <iostream>

/*
 * Abstract class KalmanFilter that describes a general Kalman filter 
 * with its properties. Any kalman filter has a filterModel that describes
 * the model. The protected methods are the interface layer between the 
 * filterModel and KalmanFilter subclasses such as discreteDiscreteKalmanFilter.
 * The public getter methods are for the subclass users. The pure virtual methods
 * are the predict and update phase as that depends on the specific filter type.
 */
template<class VECTOR,class MATRIX, class T1>
class KalmanFilter{

	filterModel<VECTOR,MATRIX,T1> *model;

	T1 currentTime;

	VECTOR currentEstimate;
	MATRIX currentCovariance;
	
	protected:

		int getStateCount(){
			return model->getStateCount();
		}

		int getSensorCount(){
			return model->getSensorCount();
		}

		MATRIX getTransitionJacobian(VECTOR v, T1 t){
			return model->transitionJacobian(v,t);
		}	

		MATRIX getMeasurementJacobian(VECTOR v, T1 t){
			return model->measurementJacobian(v,t);
		}

		VECTOR transition(VECTOR t, T1 time){
			return model->transition(t,time);
		}		

		VECTOR measure(VECTOR t, T1 time){
			return model->measurement(t,time);
		}


		void setCurrentCovariance(MATRIX t){
			currentCovariance = t;
		}

		void setCurrentTime(T1 t){
			currentTime=t;
		}

		void setCurrentEstimate(VECTOR t){
			currentEstimate = t;
		}

		MATRIX getProcessNoiseCovariance(VECTOR v, T1 t){
			return model->getProcessNoiseCovariance(v,t);
		}

		MATRIX getSensorNoiseCovariance(VECTOR v, T1 t){
			return model->getSensorNoiseCovariance(v,t);			
		}

		MATRIX getProcessNoiseCovarianceSqrt(VECTOR v, T1 t){
			return model->getProcessNoiseCovarianceSqrt(v,t);
		}

		MATRIX getSensorNoiseCovarianceSqrt(VECTOR v, T1 t){
			return model->getSensorNoiseCovarianceSqrt(v,t);			
		}

	public:

		MATRIX getCurrentCovariance(){
			return currentCovariance;
		}	
		
		VECTOR getCurrentEstimate(){
			return currentEstimate;
		}

		T1 getCurrentTime(){
			return currentTime;
		}

		KalmanFilter(
			T1 initialTime,
			VECTOR initialEstimate,
			MATRIX initialCovariance,
			filterModel<VECTOR,MATRIX,T1> *model):
			model(model){
				currentTime = initialTime;
				currentEstimate = initialEstimate;
				currentCovariance = initialCovariance;
			};

		virtual void predict(T1 timeUnitsForward) = 0;
		virtual void update(VECTOR measurement) = 0;
};


/*
 * 
 */
template<class VECTOR, class MATRIX>
class discreteDiscreteKalmanFilter: public KalmanFilter<VECTOR,MATRIX,int>{
	using KF = KalmanFilter<VECTOR,MATRIX,int>;
	
	public:
		discreteDiscreteKalmanFilter(
			int initialTime,
			VECTOR initialEstimate,
			MATRIX initialCovariance,
			filterModel<VECTOR,MATRIX,int> *model):
			KF::KalmanFilter(initialTime,
			initialEstimate,initialCovariance,model){}

		void predict(int timeUnitsForward){//TODO - have to run things time steps forward
			VECTOR est = KF::getCurrentEstimate();
			int time = KF::getCurrentTime();  
			MATRIX jacob = KF::getTransitionJacobian(est,KF::getCurrentTime());
			MATRIX covariance = KF::getCurrentCovariance();

			KF::setCurrentEstimate(KF::transition(est,time));
			KF::setCurrentCovariance(jacob*covariance*jacob.transpose() + KF::getProcessNoiseCovariance(est,time));
			KF::setCurrentTime(time + 1);
		}

		void update(VECTOR measurement){
			VECTOR est = KF::getCurrentEstimate();
			MATRIX covariance = KF::getCurrentCovariance();
			int time = KF::getCurrentTime();  
			//TODO: prevent/fix inverse
			MATRIX measurementJacobian = KF::getMeasurementJacobian(est,time);
			MATRIX Sk = KF::getSensorNoiseCovariance(est,time) + measurementJacobian*covariance*measurementJacobian.transpose();
			MATRIX KalmanGain = covariance*measurementJacobian.transpose()*Sk.inverse();
			KF::setCurrentEstimate(est + KalmanGain*(measurement - KF::measure(est,time)));
			KF::setCurrentCovariance((MATRIX::Identity(KF::getStateCount(),KF::getStateCount())-KalmanGain*measurementJacobian)*covariance);
		}
};


template<class VECTOR, class MATRIX>
class discreteDiscreteKalmanFilterSquareRoot: public KalmanFilter<VECTOR,MATRIX,int>{
	using KF = KalmanFilter<VECTOR,MATRIX,int>;
	
	
	MATRIX getCurrentCovariance(){
		return currentCovarianceSQRT*currentCovarianceSQRT.transpose();
	}	

	private:
		MATRIX currentCovarianceSQRT;

		MATRIX getCurrentCovarianceSQRT(){
			return currentCovarianceSQRT;
		}

	public:
		discreteDiscreteKalmanFilterSquareRoot(
			int initialTime,
			VECTOR initialEstimate,
			MATRIX initialCovariance,
			filterModel<VECTOR,MATRIX,int> *model):
			KF::KalmanFilter(initialTime,
			initialEstimate,initialCovariance,model){
				currentCovarianceSQRT = initialCovariance;
				MATRIX::lowerTriangulate(currentCovarianceSQRT);
			}

		void predict(int timeUnitsForward){
			int stateCount = KF::getStateCount();
			int sensorCount = KF::getSensorCount();
			VECTOR est = KF::getCurrentEstimate();
			int time = KF::getCurrentTime();  
			MATRIX jacob = KF::getTransitionJacobian(est,KF::getCurrentTime());
			MATRIX covarianceSQRT = getCurrentCovarianceSQRT();

			MATRIX tmp(stateCount,stateCount*2);//TODO: empty constructor
			tmp.block(0,0,stateCount,stateCount) = jacob*covarianceSQRT;
			tmp.block(0,stateCount,stateCount,stateCount) = KF::getProcessNoiseCovarianceSqrt(est,time);

			MATRIX::lowerTriangulate(tmp);

			KF::setCurrentEstimate(KF::transition(est,time));
			currentCovarianceSQRT = tmp.block(0,0,stateCount,stateCount);
			KF::setCurrentTime(time + 1);
		}

		void update(VECTOR measurement){
			int stateCount = KF::getStateCount();
			int sensorCount = KF::getSensorCount();
			VECTOR est = KF::getCurrentEstimate();
			MATRIX covariance = KF::getCurrentCovariance();
			int time = KF::getCurrentTime();  
			MATRIX measurementJacobian = KF::getMeasurementJacobian(est,time);

			MATRIX tmp(stateCount + sensorCount, stateCount + sensorCount);
			tmp.block(0,0,sensorCount,sensorCount) = KF::getSensorNoiseCovarianceSqrt(est,time); 
			tmp.block(0,sensorCount,sensorCount,stateCount) = measurementJacobian*currentCovarianceSQRT;
			tmp.block(sensorCount,sensorCount,stateCount,stateCount) = currentCovarianceSQRT;

			MATRIX::lowerTriangulate(tmp);
			MATRIX tmp2 = tmp.block(0,0,sensorCount,sensorCount);
			VECTOR ek = VECTOR::solve(tmp2,(measurement - KF::measure(est,time)));//TODO - is this correct?
			KF::setCurrentEstimate(est+tmp.block(sensorCount,0,stateCount,sensorCount)*ek);
			currentCovarianceSQRT = tmp.block(sensorCount,sensorCount,stateCount,stateCount);
		}
};
#endif
