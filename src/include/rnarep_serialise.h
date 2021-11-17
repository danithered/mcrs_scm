#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp> 
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>

//declare that split serialize to save/load
BOOST_SERIALIZATION_SPLIT_FREE(rnarep::CellContent)

//declare version
BOOST_CLASS_VERSION(rnarep::CellContent, 1)

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
		
		ar	<< BOOST_SERIALIZATION_NVP(mfe)
			<< BOOST_SERIALIZATION_NVP(Pfold)
			<< BOOST_SERIALIZATION_NVP(Pdeg)
			<< BOOST_SERIALIZATION_NVP(R)
			<< BOOST_SERIALIZATION_NVP(empty)
			<< BOOST_SERIALIZATION_NVP(no_sites)
			<< BOOST_SERIALIZATION_NVP(no_acts)
			<< BOOST_SERIALIZATION_NVP(type)
			<< BOOST_SERIALIZATION_NVP(prev_type)
			<< BOOST_SERIALIZATION_NVP(annot_level);

		//adding activities
		ar << boost::serialization::make_nvp("activities", boost::serialization::make_array(a, par_noEA));

	}

	template<class Archive>
	void load(Archive & ar, rnarep::CellContent & repl, unsigned int){
		ar >> BOOST_SERIALIZATION_NVP(repl.seq);

		//loading str (-> char*)
		{
			str::string str;
			ar >> BOOST_SERIALIZATION_NVP(str);
			std::strcpy(repl.str, str.c_str()); //no need to allocate, constructor did it
		}

		//other properties
		ar	>> BOOST_SERIALIZATION_NVP(mfe)
			>> BOOST_SERIALIZATION_NVP(Pfold)
			>> BOOST_SERIALIZATION_NVP(Pdeg)
			>> BOOST_SERIALIZATION_NVP(R)
			>> BOOST_SERIALIZATION_NVP(empty)
			>> BOOST_SERIALIZATION_NVP(no_sites)
			>> BOOST_SERIALIZATION_NVP(no_acts)
			>> BOOST_SERIALIZATION_NVP(type)
			>> BOOST_SERIALIZATION_NVP(prev_type)
			>> BOOST_SERIALIZATION_NVP(annot_level);

		//loading activities
		ar >> boost::serialization::make_nvp("activities", boost::serialization::make_array(a, par_noEA));
	}

}} // namespace boost::serialization

