#ifndef _CM_SERIALISE_
#define _CM_SERIALISE_

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp> 
#include <boost/core/nvp.hpp>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <stdexcept>
#include "cm.h"
#include "rnarep.h"
#include "rnarep_serialise.h"

//declare to split serilize to save/load
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::CompartPool)
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::Compart)

//serialisers
namespace boost { namespace serialization {

	//CompartPool
	template<class Archive>
	void save(Archive & ar, const cmdv::CompartPool & sim, unsigned int){
		/* 
		 * savedir <- should be loaded from pars, than -> ...
		 * output <- inic_output
		 * also dont forget to update static no_alive!!!
		 * rnarep::CellContetn::no_replicators also needs to be updated
		 * */

		//save parameters
		ar 	<< BOOST_SERIALIZATION_NVP(par_noEA)
		 	<< BOOST_SERIALIZATION_NVP(par_quit)
		 	<< BOOST_SERIALIZATION_NVP(par_maxtime)
		 	<< BOOST_SERIALIZATION_NVP(par_poolsize)
		 	<< BOOST_SERIALIZATION_NVP(par_output_interval)
		 	<< BOOST_SERIALIZATION_NVP(par_save_interval)
		 	<< BOOST_SERIALIZATION_NVP(par_claimNorep)
		 	<< BOOST_SERIALIZATION_NVP(par_splitfrom)
		 	<< BOOST_SERIALIZATION_NVP(par_num_input_content)
		 	<< BOOST_SERIALIZATION_NVP(par_ID)
		 	<< BOOST_SERIALIZATION_NVP(par_str_pool)
		 	<< BOOST_SERIALIZATION_NVP(par_outdir)
		 	<< BOOST_SERIALIZATION_NVP(par_output_filename)
		 	<< BOOST_SERIALIZATION_NVP(par_savedir)
		 	<< BOOST_SERIALIZATION_NVP(par_load)
		 	<< BOOST_SERIALIZATION_NVP(par_seed_file)
		 	<< BOOST_SERIALIZATION_NVP(par_bubbles)
		 	<< BOOST_SERIALIZATION_NVP(par_init_grid)
		 	<< BOOST_SERIALIZATION_NVP(par_ll)
		 	<< BOOST_SERIALIZATION_NVP(par_sigma)
		 	<< BOOST_SERIALIZATION_NVP(par_substitution)
		 	<< BOOST_SERIALIZATION_NVP(par_insertion)
		 	<< BOOST_SERIALIZATION_NVP(par_deletion)
		 	<< BOOST_SERIALIZATION_NVP(par_g)
		 	<< BOOST_SERIALIZATION_NVP(par_b1)
		 	<< BOOST_SERIALIZATION_NVP(par_b2)
		 	<< BOOST_SERIALIZATION_NVP(par_c)
		 	<< BOOST_SERIALIZATION_NVP(par_Emin)
		 	<< BOOST_SERIALIZATION_NVP(par_gc_bonus)
		 	<< BOOST_SERIALIZATION_NVP(par_rangePdeg)
		 	<< BOOST_SERIALIZATION_NVP(par_minPdeg)
		 	<< BOOST_SERIALIZATION_NVP(par_flexPdeg);

		//save general properties
		ar	<< boost::serialization::make_nvp("size", sim.size)
			<< boost::serialization::make_nvp("time", sim.time);

		//add comparts
		ar	<< boost::serialization::make_nvp("cells", boost::serialization::make_array(sim.comparts, sim.size));

		//save counters
		ar	<< boost::serialization::make_nvp("no_last_splits", sim.no_last_splits);
		ar	<< boost::serialization::make_nvp("no_last_replicates", sim.no_last_replicates);
		ar	<< boost::serialization::make_nvp("no_last_deaths", sim.no_last_deaths);
		ar	<< boost::serialization::make_nvp("savedir", sim.savedir);



	}

	template<class Archive>
	void load(Archive & ar, cmdv::CompartPool & sim, unsigned int){
		//load parameters
		ar 	>> BOOST_SERIALIZATION_NVP(par_noEA)
		 	>> BOOST_SERIALIZATION_NVP(par_quit)
		 	>> BOOST_SERIALIZATION_NVP(par_maxtime)
		 	>> BOOST_SERIALIZATION_NVP(par_poolsize)
		 	>> BOOST_SERIALIZATION_NVP(par_output_interval)
		 	>> BOOST_SERIALIZATION_NVP(par_save_interval)
		 	>> BOOST_SERIALIZATION_NVP(par_claimNorep)
		 	>> BOOST_SERIALIZATION_NVP(par_splitfrom)
		 	>> BOOST_SERIALIZATION_NVP(par_num_input_content);

		ar	>> par_ID;
		strcat(par_ID, "_cont\0");

		ar 	>> BOOST_SERIALIZATION_NVP(par_str_pool)
		 	>> BOOST_SERIALIZATION_NVP(par_outdir)
		 	>> BOOST_SERIALIZATION_NVP(par_output_filename)
		 	>> BOOST_SERIALIZATION_NVP(par_savedir)
		 	>> BOOST_SERIALIZATION_NVP(par_load)
		 	>> BOOST_SERIALIZATION_NVP(par_seed_file)
		 	>> BOOST_SERIALIZATION_NVP(par_bubbles)
		 	>> BOOST_SERIALIZATION_NVP(par_init_grid)
		 	>> BOOST_SERIALIZATION_NVP(par_ll)
		 	>> BOOST_SERIALIZATION_NVP(par_sigma)
		 	>> BOOST_SERIALIZATION_NVP(par_substitution)
		 	>> BOOST_SERIALIZATION_NVP(par_insertion)
		 	>> BOOST_SERIALIZATION_NVP(par_deletion)
		 	>> BOOST_SERIALIZATION_NVP(par_g)
		 	>> BOOST_SERIALIZATION_NVP(par_b1)
		 	>> BOOST_SERIALIZATION_NVP(par_b2)
		 	>> BOOST_SERIALIZATION_NVP(par_c)
		 	>> BOOST_SERIALIZATION_NVP(par_Emin)
		 	>> BOOST_SERIALIZATION_NVP(par_gc_bonus)
		 	>> BOOST_SERIALIZATION_NVP(par_rangePdeg)
		 	>> BOOST_SERIALIZATION_NVP(par_minPdeg)
		 	>> BOOST_SERIALIZATION_NVP(par_flexPdeg);

		// general properties
		ar	>> boost::serialization::make_nvp("size", sim.size)
			>> boost::serialization::make_nvp("time", sim.time);

		//reserve space for comparts
		sim.comparts = new cmdv::Compart [sim.size];

		//reserve space for replicators
		replicators = new class Compart::ScmRep[sim.size*(par_splitfrom-1)+1];
		
		//add comparts
		ar	>> boost::serialization::make_nvp("cells", boost::serialization::make_array(sim.comparts, sim.size));
		for(cmdv::Compart *comp = sim.comparts, *endcomp = comp+sim.size; comp != endcomp; comp++){
			comp->parent = &sim;
		}

		//set counters
		ar	>> boost::serialization::make_nvp("no_last_splits", sim.no_last_splits);
		ar	>> boost::serialization::make_nvp("no_last_replicates", sim.no_last_replicates);
		ar	>> boost::serialization::make_nvp("no_last_deaths", sim.no_last_deaths);

	}

	//Compart 
	template<class Archive>
	void save(Archive & ar, const cmdv::Compart & cell, unsigned int){
		double metabolism = cell.get_M();
		bool is_alive = const_cast<cmdv::Compart &>(cell).alive();

		ar	<< boost::serialization::make_nvp("reps", cell.reps);
		ar	<< boost::serialization::make_nvp("alive", is_alive )
			<< boost::serialization::make_nvp("reciproc_noEA", cell.reciproc_noEA)
			<< boost::serialization::make_nvp("M", metabolism);

	}

	template<class Archive>
	void load(Archive & ar, cmdv::Compart & cell, unsigned int ver){
		if(ver < 5) throw std::runtime_error("cmdv::Compart is out of date, can not load it.\n");

		//need to save them to parent->replicators, but it is not set yet, then get a pointer an put it here. Do not forget to assign(Compart)
		ar	>> boost::serialization::make_nvp("reps", );

		bool is_alive;
		ar	>> boost::serialization::make_nvp("alive", is_alive )
			>> boost::serialization::make_nvp("reciproc_noEA", cell.reciproc_noEA);
	}

	// replicators
	
	template<class Archive>
	void save(Archive & ar, const cmdv::Compart::ScmRep * repl, unsigned int i){
		save(ar, static_cast<const rnarep::CellContent&>(*repl), i);
	}
	
}} //namespace boost::serialize

//declare version
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::Compart::ScmRep)
BOOST_CLASS_VERSION(cmdv::Compart, 5)
BOOST_CLASS_VERSION(cmdv::CompartPool, 6)
BOOST_CLASS_VERSION(cmdv::Compart::ScmRep, 3)

#endif

