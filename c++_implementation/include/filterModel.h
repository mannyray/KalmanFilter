#ifndef FILTERMODEL
#define FILTERMODEL
#include "model.h"
#include <random>

template<class VECTOR, class MATRIX, class T1>
class filterModel{
        public:
                virtual VECTOR transition(VECTOR t, T1 timeUnit) = 0;
                virtual MATRIX transitionJacobian(VECTOR t, T1 timeUnit) = 0;
                virtual VECTOR measurement(VECTOR t, T1 timeUnit) = 0;
                virtual MATRIX measurementJacobian(VECTOR t, T1 timeUnit) = 0;
                virtual MATRIX getProcessNoiseCovariance(VECTOR t,T1 time) = 0;
                virtual MATRIX getSensorNoiseCovariance(VECTOR t,T1 time) = 0;
                virtual MATRIX getProcessNoiseCovarianceSqrt(VECTOR t,T1 time) = 0;
                virtual MATRIX getSensorNoiseCovarianceSqrt(VECTOR t,T1 time) = 0;
		virtual int getSensorCount() = 0;
		virtual int getStateCount() = 0;
};

template<class VECTOR,class MATRIX>
class discreteDiscreteFilterModel: public filterModel<VECTOR,MATRIX,int>{
	private:
		discreteModel<VECTOR> *transitionModel;
		jacobianDiscrete<VECTOR,MATRIX> * transitionJac;
		discreteNoiseCovariance<VECTOR,MATRIX> *processNoise;
		discreteModel<VECTOR> *measurementModel;
		jacobianDiscrete<VECTOR,MATRIX> * measurementJac;
		discreteNoiseCovariance<VECTOR,MATRIX> *sensorNoise;
		int stateCount;
		int sensorCount;

	public:

		discreteDiscreteFilterModel(discreteModel<VECTOR> *transitionModel,
		discreteModel<VECTOR> *measurementModel,
		discreteNoiseCovariance<VECTOR,MATRIX>* processNoise,
		discreteNoiseCovariance<VECTOR,MATRIX>* sensorNoise,
		jacobianDiscrete<VECTOR,MATRIX>* transitionJac,
		jacobianDiscrete<VECTOR,MATRIX>* measurementJac,
		int stateCount,
		int sensorCount):
		transitionModel(transitionModel),measurementModel(measurementModel),
		processNoise(processNoise),sensorNoise(sensorNoise),
		transitionJac(transitionJac),measurementJac(measurementJac),
		stateCount(stateCount),sensorCount(sensorCount){}
		
		int getStateCount(){
			return stateCount;
		}

		int getSensorCount(){
			return sensorCount;
		}

		VECTOR transition(VECTOR t, int time){
			return transitionModel->function(t,time);
		}

		MATRIX transitionJacobian(VECTOR t, int time){
			return transitionJac->function(t,time);
		}

		VECTOR measurement(VECTOR t, int time){
			return measurementModel->function(t,time);
		}

		MATRIX measurementJacobian(VECTOR t, int time){
			return measurementJac->function(t,time);
		}

		MATRIX getProcessNoiseCovariance(VECTOR t, int time){
			return processNoise->function(t,time);
		}


		MATRIX getProcessNoiseCovarianceSqrt(VECTOR t, int time){
			return processNoise->sqrt(t,time);
		}

		MATRIX getSensorNoiseCovariance(VECTOR t, int time){
			return sensorNoise->function(t,time);
		}

		MATRIX getSensorNoiseCovarianceSqrt(VECTOR t, int time){
			return sensorNoise->sqrt(t,time);
		}

};




template<class VECTOR,class MATRIX>
class discreteDiscreteFilterSolver{
	private:
		discreteDiscreteFilterModel<VECTOR,MATRIX> * model;
		VECTOR currentState;
		VECTOR measurement;
		int time;

		VECTOR * solvedMeasurements;
		VECTOR * solvedStates;
		bool usedMeasurementsPointer;
		bool usedStatesPointer;
		uint32_t seed;


	public:
		discreteDiscreteFilterSolver(discreteDiscreteFilterModel<VECTOR,MATRIX> *model):model(model){
			solvedMeasurements = NULL;
			solvedStates = NULL;
			usedMeasurementsPointer = true;
			usedStatesPointer = true;
		}

		discreteDiscreteFilterSolver(discreteDiscreteFilterModel<VECTOR,MATRIX> * model,int initialTime, VECTOR initialCondition)
		:model(model),time(initialTime),currentState(initialCondition){
			uint32_t seed = std::random_device()();
			solvedMeasurements = NULL;
			solvedStates = NULL;
			usedMeasurementsPointer = true;
			usedStatesPointer = true;
		}

		VECTOR * getSolvedMeasurements(){
			usedMeasurementsPointer = true;//users responsibility to delete now
			return solvedMeasurements;
		}

		VECTOR * getSolvedStates(){
			usedStatesPointer = true;
			return solvedStates;
		}
		
		void setSeed(uint32_t s){
			seed = s;
		}

		uint32_t getSeed(){
			return seed;
		}

		VECTOR * solve(int steps){
			int stateCount = model->getStateCount();
			int sensorCount = model->getSensorCount();
			if(usedMeasurementsPointer==false){
				delete [] solvedMeasurements;//TODO mention it is users responsibility to delete 
			}
			if(usedStatesPointer==false){
				delete [] solvedStates;
			}
			usedMeasurementsPointer = false;
			usedStatesPointer = false;

			solvedMeasurements = new VECTOR[steps];

			solvedStates = new VECTOR[steps];


			for(int i = 0; i < steps; i++){
				currentState = model->transition(currentState,time + i);
				currentState = currentState +
					model->getProcessNoiseCovarianceSqrt(currentState,time + i)*VECTOR::randomVector(stateCount,seed+i);
				measurement = model->measurement(currentState,time + i);
				measurement = measurement + 
					model->getSensorNoiseCovarianceSqrt(currentState,time)*VECTOR::randomVector(sensorCount,seed-i);
				
				solvedStates[i] = currentState;
				solvedMeasurements[i] = measurement;
			}
			
		}
};
#endif
