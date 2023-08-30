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
#include <boost/serialization/vector.hpp>
#include <stdexcept>
#include "cm.h"
#include "rnarep.h"
#include "rnarep_serialise.h"

//declare to split serilize to save/load
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::CompartPool)
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::Compart)
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::Compart*)
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::Compart**)

//serialisers
namespace boost { namespace serialization {

	//CompartPool
/*	template<class Archive>
	void serialize(Archive & ar, const cmdv::CompartPool & sim, unsigned int){
		std::string ID(par_ID), str_pool(par_str_pool), outdir(par_outdir), output_filename(par_output_filename);
		std::string savedir(par_savedir) ;
		std::string load(par_load) ;
		std::string seed_file(par_seed_file) ;
		std::string bubbles(par_bubbles) ;

		//save parameters
		ar 	 & BOOST_SERIALIZATION_NVP(par_noEA)
		 	 & BOOST_SERIALIZATION_NVP(par_quit)
		 	 & BOOST_SERIALIZATION_NVP(par_maxtime)
		 	 & BOOST_SERIALIZATION_NVP(par_poolsize)
		 	 & BOOST_SERIALIZATION_NVP(par_output_interval)
		 	 & BOOST_SERIALIZATION_NVP(par_save_interval)
		 	 & BOOST_SERIALIZATION_NVP(par_claimNorep)
		 	 & BOOST_SERIALIZATION_NVP(par_splitfrom)
		 	 & BOOST_SERIALIZATION_NVP(par_num_input_content)
		 	 & boost::serialization::make_nvp("par_ID", ID)
		 	 & boost::serialization::make_nvp("par_str_pool", str_pool)
		 	 & boost::serialization::make_nvp("par_outdir", outdir)
		 	 & boost::serialization::make_nvp("par_output_filename", output_filename)
		 	 & boost::serialization::make_nvp("par_savedir", savedir)
		 	 & boost::serialization::make_nvp("par_load", load)
		 	 & boost::serialization::make_nvp("par_seed_file", seed_file)
		 	 & boost::serialization::make_nvp("par_bubbles", bubbles)
		 	 & BOOST_SERIALIZATION_NVP(par_init_grid)
		 	 & BOOST_SERIALIZATION_NVP(par_ll)
		 	 & BOOST_SERIALIZATION_NVP(par_sigma)
		 	 & BOOST_SERIALIZATION_NVP(par_substitution)
		 	 & BOOST_SERIALIZATION_NVP(par_insertion)
		 	 & BOOST_SERIALIZATION_NVP(par_deletion)
		 	 & BOOST_SERIALIZATION_NVP(par_g)
		 	 & BOOST_SERIALIZATION_NVP(par_b1)
		 	 & BOOST_SERIALIZATION_NVP(par_b2)
		 	 & BOOST_SERIALIZATION_NVP(par_c)
		 	 & BOOST_SERIALIZATION_NVP(par_Emin)
		 	 & BOOST_SERIALIZATION_NVP(par_gc_bonus)
		 	 & BOOST_SERIALIZATION_NVP(par_rangePdeg)
		 	 & BOOST_SERIALIZATION_NVP(par_minPdeg)
		 	 & BOOST_SERIALIZATION_NVP(par_flexPdeg);

		strcpy(par_ID, ID.c_str());
		strcat(par_ID, "_cont\0");
		strcpy(par_str_pool, str_pool.c_str());
		strcpy(par_outdir, outdir.c_str());
		strcpy(par_output_filename, output_filename.c_str());
		strcpy(par_savedir, savedir.c_str());
		strcpy(par_load, load.c_str());
		strcpy(par_seed_file, seed_file.c_str());
		strcpy(par_bubbles, bubbles.c_str());

		//save general properties
		ar	& boost::serialization::make_nvp("size", sim.size)
			& boost::serialization::make_nvp("time", sim.time);

		//add comparts
		ar	& boost::serialization::make_nvp("cells", boost::serialization::make_array(sim.comparts, sim.size));

		//save counters
		ar	& boost::serialization::make_nvp("no_last_splits", sim.no_last_splits);
		ar	& boost::serialization::make_nvp("no_last_replicates", sim.no_last_replicates);
		ar	& boost::serialization::make_nvp("no_last_deaths", sim.no_last_deaths);
		ar	& boost::serialization::make_nvp("savedir", sim.savedir);

		ar 	& boost::serialization::make_nvp("replicators", boost::serialization::make_array(sim.replicators, rnarep::CellContent::no_replicators) );
	}
*/
	template<class Archive>
	void save(Archive & ar, const cmdv::CompartPool & sim, unsigned int){
		std::string ID(par_ID), str_pool(par_str_pool), outdir(par_outdir), output_filename(par_output_filename);
		std::string savedir(par_savedir) ;
		std::string load(par_load) ;
		std::string seed_file(par_seed_file) ;
		std::string bubbles(par_bubbles) ;


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
		 	<< boost::serialization::make_nvp("par_ID", ID)
		 	<< boost::serialization::make_nvp("par_str_pool", str_pool)
		 	<< boost::serialization::make_nvp("par_outdir", outdir)
		 	<< boost::serialization::make_nvp("par_output_filename", output_filename)
		 	<< boost::serialization::make_nvp("par_savedir", savedir)
		 	<< boost::serialization::make_nvp("par_load", load)
		 	<< boost::serialization::make_nvp("par_seed_file", seed_file)
		 	<< boost::serialization::make_nvp("par_bubbles", bubbles)
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

		ar	<< boost::serialization::make_nvp("no_replicators", rnarep::CellContent::no_replicators);

		ar 	<< boost::serialization::make_nvp("replicators", boost::serialization::make_array(sim.replicators, sim.size*(par_splitfrom-1)+1) );

		//add comparts
		ar	<< boost::serialization::make_nvp("cells", boost::serialization::make_array(sim.comparts, sim.size));

		//save counters
		ar	<< boost::serialization::make_nvp("no_last_splits", sim.no_last_splits);
		ar	<< boost::serialization::make_nvp("no_last_replicates", sim.no_last_replicates);
		ar	<< boost::serialization::make_nvp("no_last_deaths", sim.no_last_deaths);
		ar	<< boost::serialization::make_nvp("savedir", sim.savedir);

		//ar 	<< boost::serialization::make_nvp("rep_stack", sim.rep_stack); // solved in load!!
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

		{
		std::string ID;
		ar	>> boost::serialization::make_nvp("par_ID", ID);
		strcpy(par_ID, ID.c_str());
		strcat(par_ID, "_cont\0");
		}

		{
		std::string str_pool;
		ar	>> boost::serialization::make_nvp("par_str_pool", str_pool);
		strcpy(par_str_pool, str_pool.c_str());
		}

		{
		std::string outdir;
		ar	>> boost::serialization::make_nvp("par_outdir", outdir);
		strcpy(par_outdir, outdir.c_str());
		}

		{
		std::string output_filename;
		ar	>> boost::serialization::make_nvp("par_output_filename", output_filename);
		strcpy(par_output_filename, output_filename.c_str());
		}

		{
		std::string savedir;
		ar	>> boost::serialization::make_nvp("par_savedir", savedir);
		strcpy(par_savedir, savedir.c_str());
		}

		{
		std::string load;
		ar	>> boost::serialization::make_nvp("par_load", load);
		strcpy(par_load, load.c_str());
		}

		{
		std::string seed_file;
		ar	>> boost::serialization::make_nvp("par_seed_file", seed_file);
		strcpy(par_seed_file, seed_file.c_str());
		}

		{
		std::string bubbles;
		ar	>> boost::serialization::make_nvp("par_bubbles", bubbles);
		strcpy(par_bubbles, bubbles.c_str());
		}

		ar  	>> BOOST_SERIALIZATION_NVP(par_init_grid)
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
		unsigned int no_reps;
		ar	>> boost::serialization::make_nvp("size", sim.size)
			>> boost::serialization::make_nvp("time", sim.time);

		ar	>> boost::serialization::make_nvp("no_replicators", no_reps);

		//reserve space for replicators
		sim.replicators = new class cmdv::Compart::ScmRep[sim.size*(par_splitfrom-1)+1];
		ar 	>> boost::serialization::make_nvp("replicators", boost::serialization::make_array(sim.replicators, sim.size*(par_splitfrom-1)+1) );

		//reserve space for comparts
		sim.comparts = new cmdv::Compart* [sim.size];
		for(cmdv::Compart **comp = sim.comparts, **endcomp = sim.comparts + sim.size; comp != endcomp; ++comp){
			*comp = new cmdv::Compart;
		}

		
		//add comparts
		ar	>> boost::serialization::make_nvp("cells", boost::serialization::make_array(sim.comparts, sim.size));
		//ar	>> boost::serialization::make_nvp("cells", sim.comparts);
		for(cmdv::Compart **comp = sim.comparts, **endcomp = sim.comparts + sim.size; comp != endcomp; ++comp){
			(*comp)->parent = &sim;
		}

		//set counters
		ar	>> boost::serialization::make_nvp("no_last_splits", sim.no_last_splits);
		ar	>> boost::serialization::make_nvp("no_last_replicates", sim.no_last_replicates);
		ar	>> boost::serialization::make_nvp("no_last_deaths", sim.no_last_deaths);

		ar	>> boost::serialization::make_nvp("savedir", sim.savedir);

		rnarep::CellContent::no_replicators = no_reps;
	}







	//Compart 
	template<class Archive>
	void save(Archive & ar, const cmdv::Compart & cell, unsigned int){
		double metabolism = cell.get_M();
		bool is_alive = const_cast<cmdv::Compart &>(cell).alive();

		//ar	<< boost::serialization::make_nvp("parent", cell.parent);
		ar	<< boost::serialization::make_nvp("reps", cell.reps);
		ar	<< boost::serialization::make_nvp("alive", is_alive )
			<< boost::serialization::make_nvp("reciproc_noEA", cell.reciproc_noEA)
			<< boost::serialization::make_nvp("M", metabolism);

	}

	template<class Archive>
	void load(Archive & ar, cmdv::Compart & cell, unsigned int ver){
		if(ver < 5) throw std::runtime_error("cmdv::Compart is out of date, can not load it.\n");

		//ar	&  boost::serialization::make_nvp("parent", cell.parent);

		ar	>> boost::serialization::make_nvp("reps", cell.reps);
		for(auto &rep : cell.reps) {
			rep->assignCompart(&cell);
		}

		bool is_alive;
		double metabolism;
		ar	>> boost::serialization::make_nvp("alive", is_alive )
			>> boost::serialization::make_nvp("reciproc_noEA", cell.reciproc_noEA)
			>> boost::serialization::make_nvp("M", metabolism);
	}

	// replicators
	
	template<class Archive>
	void save(Archive & ar, const cmdv::Compart::ScmRep * repl, unsigned int i){
		save(ar, static_cast<const rnarep::CellContent&>(*repl), i);
	}

	template<class Archive>
	void save(Archive & ar, cmdv::Compart::ScmRep * repl, unsigned int i){
		save(ar, static_cast<rnarep::CellContent&>(*repl), i);
	}
	
/*	template<class Archive>
	void serialize(Archive & ar, std::vector<cmdv::Compart::ScmRep *> repvec, unsigned int i){
		//save(ar, static_cast<rnarep::CellContent&>(*repl), i);
		ar & boost::serialization::make_array(&repvec[0], repvec.size());
	}
*/
}} //namespace boost::serialize

//declare version
BOOST_SERIALIZATION_SPLIT_FREE(cmdv::Compart::ScmRep)
BOOST_CLASS_VERSION(cmdv::Compart, 5)
BOOST_CLASS_VERSION(cmdv::CompartPool, 6)
BOOST_CLASS_VERSION(cmdv::Compart::ScmRep, 3)

#endif

