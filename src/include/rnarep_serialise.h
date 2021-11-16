#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp> 
#include <boost/serialization/nvp.hpp>

namespace boost {
namespace serialization {

	//CellContent 
	template<class Archive>
	void serialize(Archive & ar, rnarep::CellContent & r, const unsigned int version)
	{
	    ar & BOOST_SERIALIZATION_NVP(r.seq);
	    ar & BOOST_SERIALIZATION_NVP(r.); //char* cant be imported easyly -> has to be spllitted to save/load
	}

} // namespace serialization
} // namespace boost
