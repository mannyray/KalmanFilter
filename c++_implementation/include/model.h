#ifndef MODEL_H
#define MODEL_H

#include <assert.h>
#include <stdexcept>

/*
 * Empty class designed for representing a model. Inherited types
 *  include discrete and continuos version. K represents time unit. 
 *  The reason for various sub classes is for more user to clearly
 *  understand what models they are working with.
 */
template <class K>
class model{};

/**
 * Abstract continuousModel inherits from model<double> and thus specifies
 * that time is continuous. Class only declares single pure virtual
 * method 'function' that describes a ODE where returned value is
 * is the derivative of T evaluated at a specific time and value T.
 * T is the data type whether double, Eigen or custom implemented.
 *
 * Following two examples implement the same differential equation 
 * but using different T's. 
 *
 * \section ex1 Example using double
 * \snippet model/continuousModel1.cpp example1
 * \section ex2 Example using Eigen
 * \snippet model/continuousModel2.cpp example1
 *
 */
template <typename T> class continuousModel: public model<double>{
	public:
                /** 
		 * Pure virtual function. User is responsible for implementing 
		 * this method in inherited class. The function returns
		 * rhs (right hand side) where rhs = diff(T,t) = function(T,t)
		 * where rhs is the derivative of variable T at time t.
		 *
		 * User is also responsible for specifying data type T. 
		 */  
		virtual T function(const T & val, const double time) const = 0;
};

/**
 * Abstract discreteModel inherits from model<int> and thus specifies
 * that time is 'discrete'. Class only declares single pure virtual
 * method 'function' that describes an incremental function where
 * the current state at time k that is T_k is given by function(T_{k-1},k-1). 
 * T is the data type whether double, Eigen or custom implemented.
 *
 * \section ex1 Example using double
 * \snippet model/discreteModel1.cpp example1
 * \section ex2 Example using Eigen
 * \snippet model/discreteModel2.cpp example1
 */
template <typename T> class discreteModel: public model<int>{
        public:
                /** 
		 * Pure virtual function. User is responsible for implementing 
		 * this method in inherited class. The function returns
		 * rhs (right hand side) where rhs = T{k} = function(T_{k-1},k-1)
		 * where rhs is the next model state after state being T_{k-1}
		 * at time k-1.
		 *
		 * User is also responsible for specifying data type T. 
		 */  
                virtual T function(const T& val, const int index) const = 0;
};

/*
 * Abstract jacobian that defines a general Jacobian class. Class declares
 * a single method `function` that returns the jacobian of a variable val
 * evaluated at time T. The function with which the Jacobian is associated 
 * here is not defined here and it is the users responsibility that the Jacobian
 * is correct.  
 * This is akin to the model<K> type class from which continuousModel<G> and 
 * discreteModel<G> inherit from.
 * The reason for various sub classes is for more user to clearly
 * understand what jacobians they are working with.
 */
template<class VECTOR, class MATRIX, class T>
class jacobian{
	public:
		///Pure virtual method returning Jacobian of some function.
		virtual MATRIX function(const VECTOR& val, const T t) = 0;
};


/*
 * Abstract jacobianContinuous. Class declares
 * a single method `function` that returns the jacobian of a variable val
 * evaluated at time T. The function with which the Jacobian is associated 
 * here is not defined here and it is the users responsibility that the Jacobian
 * is correct. The jacobianContinuous in other classes is associated with a continuousModel.  
 * See continuousModel<K> for an example of what the jacobianContinuous<K> implementation
 * would look like.
 */
template <class VECTOR, class MATRIX>
class jacobianContinuous: public jacobian<VECTOR,MATRIX,double>{
	public:
		///Pure virtual method returning Jacobian of some function.
		virtual MATRIX function(const VECTOR & val, double t) = 0; 
};



/*
 * Abstract jacobianDiscrete. Class declares
 * a single method `function` that returns the jacobian of a variable val
 * evaluated at time T. The function with which the Jacobian is associated 
 * here is not defined here and it is the users responsibility that the Jacobian
 * is correct. The jacobianDiscrete in other classes is associated with a discreteModel.  
 * See continuousModel<K> for an example of what the jacobianContinuous<K> implementation
 * would look like.
 */
template <class VECTOR, class MATRIX>
class jacobianDiscrete: public jacobian<VECTOR,MATRIX,int>{
	public:
		///Pure virtual method returning Jacobian of some function.
		virtual MATRIX function(const VECTOR & val, int t) = 0; 
};


/*
 * Abstract class akin to continuousModel<K> and discreteModel<K> 
 * where 'function' only takes time T argument and returns type M
 * argument.
 */
template <class M, class T1>
class timeBasedFunction{
	public:
		//Pure virtual method.
		virtual M const function(T1 t) = 0;	
};


template<class VECTOR, class MATRIX, class T>
class noiseCovariance{
	public:
		///Pure virtual method returning Jacobian of some function.
		virtual MATRIX function(const VECTOR& val, const T t) = 0;
		virtual MATRIX sqrt(const VECTOR & val, int t) = 0;
};


template <class VECTOR, class MATRIX>
class discreteNoiseCovariance: public noiseCovariance<VECTOR,MATRIX,int>{
	public:
		///Pure virtual method returning Jacobian of some function.
		virtual MATRIX function(const VECTOR & val, int t) = 0; 
		virtual MATRIX sqrt(const VECTOR & val, int t) = 0;
};


/*
 * Non abstract solver for solving discrete models. 
 * Unlike continuousSolver<T> there is no abstract base class
 * for the discrete case as there is no variance for
 * discrete solvers since all the hard work is encoded in the 
 * discreteModel. Despite this, both discreteSolver<T>
 * and continuousSolver<T> have very similar methods and variables
 * raising the question as to why they do not inherit from the same method.
 * This was an initial 'mistake' and is now the way the code is.
 */
template<class T> class discreteSolver{
	private:
		///discreteSolver solves a discreteModel
		const discreteModel<T> *M;

		///Current time/index in model.
		int currentIndex;

		///Current state of model at time currentIndex.
		T currentState;

	public:
		/*
		 * \param startIndex initial time 
		 * \param startState initial state of model at startIndex
		 * Constructor.
		 */
                discreteSolver(const discreteModel<T> *M, const int startIndex, const T& startState)
			:M(M),currentIndex(startIndex),currentState(startState){}

		/*
		 * Constructor for when no model provided. No check is conducted to
		 * see if user ran setInitialConditions before running solve - in the
		 * case when user does that the behaviour is undefined. 
		 */
		discreteSolver(discreteModel<T> *M):M(M){}

                ///solve timeForward time units forward from currentTime
                void solve(int stepsForward){
			if(stepsForward < 0){
				throw std::out_of_range("");
			}
			for(int i = 0; i < stepsForward; i++){
				currentState = M->function(currentState,currentIndex);
				currentIndex++;
			}
		}

		/* 
		 * Setter of time and current model state included in order for easier solver reuse.
		 */
		void setInitialConditions(const T & newState, int newIndex){
			currentState = newState;
			currentIndex = newIndex;
		}

		///Getter for current state.
                T getCurrentState(){
                        return currentState;
                }

		///Getter for current time.
                double getCurrentIndex(){
                        return currentIndex;
                }

                /*
		 * Interface for interacting with model. Does not affect
		 * the solver's stored current time and current state.
		 */
                T evaluateModel(const T & t,const double time){
                        return M->function(t,time);
                }
};

/*
 * Abstract numerical solver of ODEs represented by continuousModel<T>
 */
template<class T>class continuousSolver{
        private:
		///continuousSolver solves a continuousModel
                const continuousModel<T> *M;
		
		///current time in the model
                double currentTime;
		
		///current state of model
                T currentState;

		//TODO:destructor/constructor? 
        public:
		/*
		 * \param startTime initial time 
		 * \param startState initial state of model at initial time
		 * \param continuousModel encoding ODE 
		 * Constructor.
		 */
                continuousSolver(const continuousModel<T> *M, const double startTime, const T& startState)
			:M(M),currentTime(startTime),currentState(startState){};

                /*
		 * Solve timeForward time units forward from currentTime. 
		 * This pure virtual method leaves the detail of solving the 
		 * continuousModel up to the inheriting class - be it Eulers, RK4
		 * or whatever.
		 */
                virtual void solve(double timeForward) = 0;

		/* 
		 * Setter of time and current model state included in order for easier solver reuse.
		 */
		void setInitialConditions(const T & newState, double newTime){
			currentState = newState;
			currentTime = newTime;
		}
		
		///Getter for state
                T getCurrentState(){
                        return currentState;
                }

		///Getter for current time.
                double getCurrentTime(){
                        return currentTime;
                }

                /*
		 * Interface for interacting with model. Does not affect
		 * the solver's stored current time and current state.
		 */
                T evaluateModel(const T & t,const double time) const{
                        return M->function(t,time);
                }
};


/*
 * Non abstract RK4 solver inheriting from continuousSolver<T>
 */
template <typename T> class continuousSolverRK4: public continuousSolver<T>{

        private:
		///RK4 has a fixed timestep set in the constructor
                double timeStep;
        public:
		/*
		 * \param startTime initial time 
		 * \param startState initial state of model at initial time
		 * \param continuousModel encoding ODE 
		 * \param timeStep the timestep RK4 will use
		 * Constructor.
		 */
                continuousSolverRK4( continuousModel<T>* M,const double startTime,const T& startState, double timeStep):
                        continuousSolver<T>(M,startTime,startState),timeStep(timeStep){
                }

		/*
		 * Implementation of solve method that is pure virtual in superclass continuousSolver. 
		 * Propagate the current state and current time according to the ODE described in continuousModel
		 * provided in the constructor. 
		 */
                void solve(double timeForward){
			
			if(timeForward < 0){
				throw std::out_of_range("");
			}

                        double localTime = continuousSolver<T>::getCurrentTime();
                        T tmp = continuousSolver<T>::getCurrentState();

                        double startTime = localTime;
                        double localTimeStep = timeStep;
                        double finalTime = startTime + timeForward;

			while(localTime < finalTime){
				//deal with unequal timestep not dividing evently into timeForward
                                if(localTime + localTimeStep > finalTime){
                                        localTimeStep = finalTime - localTime;
                                }

                                T k1 = localTimeStep*continuousSolver<T>::evaluateModel(tmp         , localTime                    );
                                T k2 = localTimeStep*continuousSolver<T>::evaluateModel(tmp + 0.5*k1, localTime + 0.5*localTimeStep);
                                T k3 = localTimeStep*continuousSolver<T>::evaluateModel(tmp + 0.5*k2, localTime + 0.5*localTimeStep);
                                T k4 = localTimeStep*continuousSolver<T>::evaluateModel(tmp +     k3, localTime +     localTimeStep);

                                tmp = tmp + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
                                localTime+=localTimeStep;
                        }
                        continuousSolver<T>::setInitialConditions(tmp,localTime);
                }
};
#endif
