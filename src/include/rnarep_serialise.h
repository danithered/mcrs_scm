#ifndef _RNAREP_SERIALISE_
#define _RNAREP_SERIALISE_

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp> 
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/array_wrapper.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/list.hpp>
#include "rnarep.h"

//serialisers
namespace boost { namespace serialization {

	//CellContent 
	template<class Archive>
	void save(Archive & ar, const rnarep::CellContent & repl, unsigned int){
		ar << BOOST_SERIALIZATION_NVP(repl.seq);

		//adding str (char array)
		{
			std::string str (repl.str);
			ar << BOOST_SERIALIZATION_NVP(str); 
		}
		
		ar	<< BOOST_SERIALIZATION_NVP(repl.mfe)
			<< BOOST_SERIALIZATION_NVP(repl.Pfold)
			<< BOOST_SERIALIZATION_NVP(repl.Pdeg)
			<< BOOST_SERIALIZATION_NVP(repl.R)
			<< BOOST_SERIALIZATION_NVP(repl.empty)
			<< BOOST_SERIALIZATION_NVP(repl.no_sites)
			<< BOOST_SERIALIZATION_NVP(repl.no_acts)
			<< BOOST_SERIALIZATION_NVP(repl.type)
			<< BOOST_SERIALIZATION_NVP(repl.prev_type)
			<< BOOST_SERIALIZATION_NVP(repl.annot_level);

		//adding activities
		ar << boost::serialization::make_nvp("activities", boost::serialization::make_array(repl.a, par_noEA));

	}

	template<class Archive>
	void load(Archive & ar, rnarep::CellContent & repl, unsigned int){
		ar >> BOOST_SERIALIZATION_NVP(repl.seq);

		//loading str (-> char*)
		{
			std::string str;
			ar >> BOOST_SERIALIZATION_NVP(str);
			std::strcpy(repl.str, str.c_str()); //no need to allocate, constructor did it
		}

		//other properties
		ar	>> BOOST_SERIALIZATION_NVP(repl.mfe)
			>> BOOST_SERIALIZATION_NVP(repl.Pfold)
			>> BOOST_SERIALIZATION_NVP(repl.Pdeg)
			>> BOOST_SERIALIZATION_NVP(repl.R)
			>> BOOST_SERIALIZATION_NVP(repl.empty)
			>> BOOST_SERIALIZATION_NVP(repl.no_sites)
			>> BOOST_SERIALIZATION_NVP(repl.no_acts)
			>> BOOST_SERIALIZATION_NVP(repl.type)
			>> BOOST_SERIALIZATION_NVP(repl.prev_type)
			>> BOOST_SERIALIZATION_NVP(repl.annot_level);

		//loading activities
		ar >> boost::serialization::make_nvp("activities", boost::serialization::make_array(repl.a, par_noEA));
	}

}} // namespace boost::serialization

//declare that split serialize to save/load
BOOST_SERIALIZATION_SPLIT_FREE(rnarep::CellContent)

//declare version
BOOST_CLASS_VERSION(rnarep::CellContent, 1)

#endif

