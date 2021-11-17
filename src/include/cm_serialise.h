#ifndef _RNAREP_SERIALISE_
#define _RNAREP_SERIALISE_

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp> 
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include "rnarep_serialise.h"

//declare to split serilize to save/load
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::CompartPool)

//declare version
BOOST_CLASS_VERSION(cmdv::Compart, 1)
BOOST_CLASS_VERSION(cmdv::CompartPool, 1)

//serialisers
namespace boost { namespace serialization {

	//CellContent 
	template<class Archive>
	void serialize(Archive & ar, cmdv::Compart & cell, unsigned int){
		ar	& BOOST_SERIALIZATION_NVP(cell.reps);
		ar	& BOOST_SERIALIZATION_NVP(cell.updateable)
			& BOOST_SERIALIZATION_NVP(cell.alive)
			& BOOST_SERIALIZATION_NVP(cell.reciproc_noEA)
			& BOOST_SERIALIZATION_NVP(cell.leftover);

		//save metabolism for analytical purposes
		{
			double metabolism = cell.M();
			ar & BOOST_SERIALIZATION_NVP(metabolism);
		}

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
			<< BOOST_SERIALIZATION_NVP(sim.time);

		//add comparts
		ar	<< boost::serialization::make_nvp("cells", boost::serialization::make_array(sim.comparts, sim.size));

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

#endif

