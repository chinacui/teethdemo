#pragma  once
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>

using std::vector;
using std::map;
using std::string;
#define  SAFE_PROB 1  //using logistical arithmetic for avoiding small probs
#define  NEGA_MIN_PORB -99999999999
template<class StateType>
class viterbi_triple_t{
public:
	double prob;
	double vprob;
	std::vector<StateType> vpath;

	viterbi_triple_t(){};
	viterbi_triple_t(double _total, std::vector<StateType> & _vpath, double _valmax)
	{
		prob = _total;
		vpath.clear();	vpath.resize(_vpath.size());
		std::copy(_vpath.begin(),_vpath.end(),vpath.begin());
		vprob = _valmax;
	}

	void print_path() const
	{
		size_t np = vpath.size();
		for (unsigned i=0; i<np-1; ++i)
			cout<<vpath[i]<<"->";
		cout<<vpath.back()<<". "<<"\n"
			<<"Observation sequence prob :"<<prob<<"\n"
			<<"hidden state sequence prob :"<<vprob<<"\n\n";
	}
};


template<class StateType>
class HMM
{
public:
	typedef typename map<StateType, double>	Map_SD;
	typedef typename map<StateType, Map_SD>	Map_SMapSD;
	typedef typename std::vector<std::vector<double> >	PM;  //prob matrix
public:
	//init probabilities
	size_t m_nState, m_nOb;
	std::vector<StateType> m_state_set;
	std::vector<double> m_init_prob;
	Map_SD m_start_prob;
	//probabilities matrix
	Map_SMapSD transition_prob, emission_prob;
	//observations
	std::vector<StateType> m_ob_set;

	//for accelaration
	PM	*m_pTP,*m_pEP;
	
public:	
	HMM(){};
	bool Init(std::vector<StateType> &state_set, std::vector<double> &init_prob,	std::vector<StateType> &ob_set, 
		std::vector<double> &tr_pro, std::vector<double> &em_pro)
	{
		//given init probabilities
		if(state_set.size() != init_prob.size()){
			std::cout<<"init states are not same size with init prob \n"; 
			return false;
		}


		if( (m_nState = state_set.size())!=0 ){
			m_state_set.clear();		m_state_set.resize(m_nState);
			std::copy(state_set.begin(), state_set.end(), m_state_set.begin());
			m_init_prob.clear();		m_init_prob.resize(m_nState);
			std::copy(init_prob.begin(), init_prob.end(), m_init_prob.begin());
		}
		else{
			std::cout<<"init state and init prob wrong \n"; 
			return false;
		}

		for ( unsigned i=0 ; i<m_nState ; ++i )
			m_start_prob.insert(Map_SD::value_type(m_state_set[i],m_init_prob[i]));

		//given observations
		if( (m_nOb = ob_set.size())!= 0 ){
			m_ob_set.clear();		m_ob_set.resize(m_nOb);
			std::copy(ob_set.begin(), ob_set.end(), m_ob_set.begin());				
		}
		else{
			std::cout<<"init observation wrong \n"; 
			return false;
		}		
		//give transition matrix
		if(tr_pro.size() == m_nState*m_nState){
			for (unsigned i=0; i<m_nState; ++i){
				Map_SD state_i_trans_probs;
				for (unsigned j=0; j<m_nState; ++j)
					state_i_trans_probs.insert( Map_SD::value_type(m_state_set[j],tr_pro[i*m_nState+j]) );
				transition_prob.insert( Map_SMapSD::value_type(m_state_set[i],state_i_trans_probs) );
			}
		}
		else{
			std::cout<<"init transition matrix wrong \n"; 
			return false;
		}

		//give emission matrix
		if(em_pro.size() == m_nState*m_nOb){
			for (unsigned i=0; i<m_nState; ++i){
				Map_SD state_i_ems_probs;
				for (unsigned j=0; j<m_nOb; ++j)
					state_i_ems_probs.insert( Map_SD::value_type(m_ob_set[j],em_pro[i*m_nOb+j]) );
				emission_prob.insert( Map_SMapSD::value_type(m_state_set[i],state_i_ems_probs) );
			}
		}
		else{
			std::cout<<"init emission matrix wrong \n"; 
			return false;
		}

		return true;
	}
	
	bool Init(int nstate, int nob, std::vector<double> &init_prob, 
			  PM *transProb, PM *emisProb)
	{
		m_nState = nstate;	m_nOb = nob;
		if(m_nState != init_prob.size()){
			std::cout<<"init init prob wrong \n";
			return false;
		}
		m_init_prob.clear();		m_init_prob.resize(m_nState);
		std::copy(init_prob.begin(), init_prob.end(), m_init_prob.begin());
		if(transProb->size()!=m_nState || emisProb->size()!=m_nState ){
			std::cout<<"prob matrix sizes are wrong \n";
			return false;
		}
		else{
			m_pTP = transProb;			m_pEP = emisProb;
		}
		return true;
	}



    /* it's a much forward version of Viterbi */
	bool ViterbiRun(viterbi_triple_t<StateType> &vtriple, std::vector<StateType> &ob_sequence)
	{
		//alias
		Map_SD& sp = m_start_prob;
		Map_SMapSD &tp = transition_prob;
		Map_SMapSD &ep = emission_prob;

		//init,
		map<StateType, viterbi_triple_t<StateType> > P;  //will hold all the optimal path for each of states
		for (std::vector<StateType>::iterator it = m_state_set.begin(); it != m_state_set.end(); ++it)
		{
			viterbi_triple_t<StateType> foo;
			foo.prob = sp[*it];
			foo.vpath.push_back(*it);
			foo.vprob = sp[*it];
			P[*it] = foo;
		}

		map<StateType, viterbi_triple_t<StateType> > U;

		double total = 0.0;
		std::vector<StateType> argmax;
		double valmax = 0;
		double p = 0;

		//for the given observation sequence , NOTE: observation set is not the observation sequence!
		//of course one could simply take the observation set as the ob-sequence, while it's not always the case
		if(ob_sequence.size()<=0) return false;
		for (std::vector<StateType>::const_iterator itob = ob_sequence.begin(); itob != ob_sequence.end(); ++itob)
		{
			cout<<"observation=" << *itob << endl;
			U.clear();

			for (std::vector<StateType>::iterator itNextState = m_state_set.begin();itNextState != m_state_set.end();++itNextState)
			{
				cout<<"\tnext_state="<< *itNextState << endl;
				total = 0;
				argmax.clear();
				valmax = 0;
				for( std::vector<StateType>::iterator itSrcState = m_state_set.begin();itSrcState != m_state_set.end();++itSrcState)
				{
					cout <<"\t\tstate=" << *itSrcState << endl;
					viterbi_triple_t<StateType> foo = P[*itSrcState];
					p = ep[*itSrcState][*itob] * tp[*itSrcState][*itNextState];
					cout << "\t\t\tp=" << p << endl;
					foo.prob *= p;
					foo.vprob *= p;
					//cout <<"\t\t\t triple="<< foo << endl;
					total += foo.prob;

					if (foo.vprob > valmax){
						foo.vpath.push_back(*itNextState);
						argmax = foo.vpath;
						valmax = foo.vprob;
					}
				}

				U[*itNextState] = viterbi_triple_t<StateType>(total, argmax, valmax);
				//cout <<"\tUpdate U["<< *itNextState << "]=" << U[*itNextState] << " " << endl;
			}
			P.swap(U);
		}

		total = 0;
		argmax.clear();
		valmax = 0;

		for (std::vector<StateType>::iterator itState = m_state_set.begin(); itState != m_state_set.end(); ++itState)
		{
			viterbi_triple_t<StateType> foo = P[*itState];
			total += foo.prob;
			if (foo.vprob > valmax){
				argmax.swap(foo.vpath);
				valmax = foo.vprob;
			}
		}

		vtriple.prob = total;
		vtriple.vpath = argmax;
		vtriple.vprob = valmax;

		cout << "final triple  : " <<endl;
		vtriple.print_path();

		return true;
	}


	bool ViterbiEvolve(viterbi_triple_t<StateType> &vtriple, std::vector<StateType> &ob_sequence, std::vector<double> ob_evolve,double delta_C )
	{
		if(ob_sequence.size()<=0) {
			std::cout<<"zero observation sequence.\n"; 
			return false;
		}
		if(ob_evolve.size() != ob_sequence.size()) {
			std::cout<<"evolve sequence number is wrong.\n"; 
			return false;
		}

		//alias
		Map_SD& sp = m_start_prob;
		Map_SMapSD &tp = transition_prob;
		Map_SMapSD &ep = emission_prob;

		//init,
		map<StateType, viterbi_triple_t<StateType> > P;  //will hold all the optimal path for each of states
		for (std::vector<StateType>::iterator it = m_state_set.begin(); it != m_state_set.end(); ++it)
		{
			viterbi_triple_t<StateType> foo;
			foo.prob =sp[*it];
			foo.vprob = sp[*it];
			P[*it] = foo;
		}

		map<StateType, viterbi_triple_t<StateType> > U;

		double total = 0.0;
		std::vector<StateType> argmax;
		double valmax = 0;
		double p = 0;

		//for the given observation sequence , NOTE: observation set is not the observation sequence!
		//of course one could simply take the observation set as the ob-sequence, while it's not always the case

		for (std::vector<StateType>::const_iterator itob = ob_sequence.begin(); itob != ob_sequence.end(); ++itob)
		{
			U.clear();
		

			for (std::vector<StateType>::iterator itNextState = m_state_set.begin();itNextState != m_state_set.end();++itNextState)
			{
				total = 0;
				argmax.clear();
				valmax = 0;
				for(std::vector<StateType>::iterator itSrcState = m_state_set.begin();itSrcState != m_state_set.end();++itSrcState)
				{
					viterbi_triple_t<StateType> foo = P[*itSrcState];
					//recify the transition probabilities based on the current evolve..
					double dC_1= (tp[*itSrcState][*itNextState]/ob_evolve[*itob]-1.0);
					double evolved_tp = exp(-0.5*dC_1*dC_1/delta_C/delta_C);
					//double evolved_tp = -dC_1*dC_1/delta_C/delta_C;
					p = ep[*itNextState][*itob] * evolved_tp ;

					foo.prob  *= p;
					foo.vprob *= p;
					total *= foo.prob;

					if (foo.vprob > valmax){
						foo.vpath.push_back(*itNextState);
						argmax = foo.vpath;
						valmax = foo.vprob;
					}
				}

				U[*itNextState] = viterbi_triple_t<StateType>(total, argmax, valmax);
			}
			P.swap(U); //update all the optimal path for every sources at this observation stage
		}

		total = 0;
		argmax.clear();
		valmax = 0;

		for (std::vector<StateType>::iterator itState = m_state_set.begin(); itState != m_state_set.end(); ++itState)
		{
			viterbi_triple_t<StateType> foo = P[*itState];
			total += foo.prob;
			if (foo.vprob > valmax){
				argmax.swap(foo.vpath);
				valmax = foo.vprob;
			}
		}

		vtriple.prob = total;
		vtriple.vpath = argmax;
		vtriple.vprob = valmax;

		return true;
	}
	


	bool ViterbiEvolve_Modified(viterbi_triple_t<int> &vtriple, std::vector<int> &ob_sequence, std::vector<double> &ob_evolve, double delta_C )
	{
		std::cerr << "start viterbi" << std::endl;
		if(ob_sequence.size()<=0) {
			std::cout<<"zero observation sequence.\n"; 
			return false;
		}
		if(ob_evolve.size() != ob_sequence.size()) {
			std::cout<<"evolve sequence number is wrong.\n"; 
			return false;
		}

		PM &transition_probM = *m_pTP;
		PM &emission_prob = *m_pEP;
		size_t nOb = ob_sequence.size();
		//init,
		std::vector< viterbi_triple_t<int> > P(m_nState);  //will hold all the optimal path for each of states

		for (unsigned i=0; i<m_nState; ++i)
		{
			viterbi_triple_t<int> foo;
			foo.prob = m_init_prob[i];
			foo.vprob = m_init_prob[i];
			P[i] = foo;
		}

		std::vector< viterbi_triple_t<int> > U(m_nState);

		double total = 0.0;
		std::vector<int> argmax;
		#if SAFE_PROB
		double valmax = NEGA_MIN_PORB;
		#else
		double valmax = 0.0;
		#endif
		double p = 0;

		//for the given observation sequence , NOTE: observation set is not the observation sequence!
		//of course one could simply take the observation set as the ob-sequence, while it's not always the case
		std::cerr << "size nob "<<nOb << std::endl;
		std::cerr << "size state " << m_nState << std::endl;
		for (size_t i_ob=0; i_ob<nOb; ++i_ob)			
		{
			U.clear();U.resize(m_nState);
			for (unsigned i_nextstate=0; i_nextstate<m_nState; ++i_nextstate)
			{
				total = 0;
				argmax.clear();
				#if SAFE_PROB
				valmax = NEGA_MIN_PORB;
				#else
				valmax = 0.0;
				#endif
				for (unsigned i_srcstate=0; i_srcstate<m_nState; ++i_srcstate)
				{
					viterbi_triple_t<int> foo = P[i_srcstate];
					//recify the transition probabilities based on the current evolve..
					double dC_1 = (transition_probM[i_srcstate][i_nextstate] / ob_evolve[ob_sequence[i_ob]] - 1.0);
					#if SAFE_PROB
					double evolved_tp = -dC_1*dC_1/delta_C/delta_C ;
					p = emission_prob[i_nextstate][ob_sequence[i_ob]] + evolved_tp;
					foo.prob += p;
					foo.vprob += p;
					total += foo.prob;
					#else
					double evolved_tp = exp(-0.5*dC_1*dC_1/delta_C/delta_C) ;
					p = emission_prob[i_nextstate][ob_sequence[i_ob]] *
						transition_probM[i_srcstate][i_nextstate];
					foo.prob *= p;
					foo.vprob *= p;
					total += foo.prob;
					#endif

					if (foo.vprob > valmax){
						foo.vpath.push_back(i_nextstate);
						argmax = foo.vpath;
						valmax = foo.vprob;
					}
				}

				U[i_nextstate] = viterbi_triple_t<int>(total, argmax, valmax);
			}
			P.swap(U); //update all the optimal path for every sources at this observation stage
		}

		total = 0;
		argmax.clear();
		#if SAFE_PROB
		valmax = NEGA_MIN_PORB;
		#else
		valmax = 0.0;
		#endif

		for (unsigned i_state=0; i_state<m_nState; ++i_state)
		{
			viterbi_triple_t<int> foo = P[i_state];
			total += foo.prob;
			if (foo.vprob > valmax){
				argmax.swap(foo.vpath);
				valmax = foo.vprob;
			}
		}

		vtriple.prob = total;
		vtriple.vpath = argmax;
		vtriple.vprob = valmax;
		std::cerr << "end viterbi" << std::endl;
		return true;
	}




};