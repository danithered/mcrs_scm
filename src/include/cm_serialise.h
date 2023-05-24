#ifndef _CM_SERIALISE_
#define _CM_SERIALISE_

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp> 
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/list.hpp>
#include "rnarep_serialise.h"

//declare to split serilize to save/load
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::CompartPool)

//serialisers
namespace boost { namespace serialization {

	//CellContent 
	template<class Archive>
	void serialize(Archive & ar, cmdv::Compart & cell, unsigned int){
		double metabolism = cell.M();
		bool is_alive = cell.alive();

		ar	& BOOST_SERIALIZATION_NVP(cell.reps);
		ar	//& BOOST_SERIALIZATION_NVP(cell.updateable)
			& BOOST_SERIALIZATION_NVP( is_alive )
			& BOOST_SERIALIZATION_NVP(cell.reciproc_noEA)
			& BOOST_SERIALIZATION_NVP(metabolism);

	}

	//CompartPool
	template<class Archive>
	void save(Archive & ar, const cmdv::CompartPool & sim, unsigned int){
		/* 
		 * savedir <- should be loaded from pars, than -> ...
		 * output <- inic_output
		 * also dont forget to update static no_alive!!!
		 * rnarep::CellContetn::no_replicators also needs to be updated
		 * */


		ar	<< BOOST_SERIALIZATION_NVP(sim.size)
			<< boost::serialization::make_nvp("time", sim.time);

		//add comparts
		ar	<< boost::serialization::make_nvp("cells", boost::serialization::make_array(sim.comparts, sim.size));

		ar	<< BOOST_SERIALIZATION_NVP(sim.no_last_splits);

	}

	template<class Archive>
	void load(Archive & ar, cmdv::CompartPool & sim, unsigned int){
		ar	>> BOOST_SERIALIZATION_NVP(sim.size)
			>> BOOST_SERIALIZATION_NVP(sim.time);

		//add comparts
		sim.comparts = new cmdv::Compart [sim.size];
		ar	>> boost::serialization::make_nvp("cells", boost::serialization::make_array(sim.comparts, sim.size));
		for(cmdv::Compart *comp = sim.comparts, *endcomp = comp+sim.size; comp != endcomp; comp++){
			comp->parent = &sim;
		}

	}
	
}} //namespace boost::serialize

//declare version
BOOST_CLASS_VERSION(cmdv::Compart, 2)
BOOST_CLASS_VERSION(cmdv::CompartPool, 3)

#endif

