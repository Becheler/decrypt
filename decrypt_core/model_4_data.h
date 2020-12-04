
// Copyright 2020 Arnaud Becheler    <arnaud.becheler@gmail.com>

/***********************************************************************                                                                         *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
************************************************************************/

#ifndef __M4_REALITY_CHECK_H_INCLUDED__
#define __M4_REALITY_CHECK_H_INCLUDED__

#include "include/quetzal.h"
#include "utils.h"

#include <boost/program_options.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "sqlite3pp.h"

#include <random>
#include <algorithm>
#include <cassert>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <functional>

namespace coal = quetzal::coalescence;
namespace geo = quetzal::geography;
namespace demography = quetzal::demography;
namespace genet = quetzal::genetics;
namespace sim = quetzal::simulator;
namespace expr = quetzal::expressive;
namespace bpo = boost::program_options;


// Returns a map with the program options
auto handle_options(int argc, char* argv[])
{
  bpo::variables_map vm;
  try
  {
    bpo::options_description generalOptions{"General"};
    generalOptions.add_options()
    ("help,h", "Help screen")
    ("config", bpo::value<std::string>(), "Config file")
    ("landscape", bpo::value<std::string>()->required(), "Geospatial file in tiff format giving the friction map");

    bpo::options_description fileOptions{"File"};
    fileOptions.add_options()
    ("help", "produce help message")
    ("ploidy", bpo::value<unsigned int>()->required(), "1 for haploid, 2 for diploid.")
    ("n_replicates", bpo::value<unsigned int>()->required(), "Number of replicates to simulate. For each replicate, a forward demographic history and loci genealogies are simulated.")
    ("n_loci", bpo::value<unsigned int>()->required(), "Number of loci to simulate")
    ("sample",  bpo::value<std::string>(), "File name for the lon/lat of sampled genetic material")
    ("lon_0", bpo::value<double>()->required(), "Introduction point longitude")
    ("lat_0", bpo::value<double>()->required(), "Introduction point latitude")
    ("N_0", bpo::value<unsigned int>()->required(), "Number of gene copies at introduction point")
    ("duration", bpo::value<unsigned int>()->required(), "Number of generations to simulate")
    ("K_suit", bpo::value<unsigned int>()->required(), "Carrying capacity in suitable areas")
    ("K_max", bpo::value<unsigned int>()->required(), "Highest carrying capacity in areas with null suitability")
    ("K_min", bpo::value<unsigned int>()->required(), "Lowest carrying capacity in areas with null suitability")
    ("p_K", bpo::value<double>()->required(), "Probability to have highest carrying capacity in areas with 0 suitability")
    ("r", bpo::value<double>()->required(), "Growth rate")
    ("emigrant_rate", bpo::value<double>()->required(), "Emigrant rate between the four neighboring cells")
    ("demography_out",  bpo::value<std::string>(), "File name for the simulated demography output")
    ("database", bpo::value<std::string>(), "Filename database storing the output");

    store(parse_command_line(argc, argv, generalOptions), vm);
    if (vm.count("config"))
    {
      std::ifstream ifs{vm["config"].as<std::string>().c_str()};
      if (ifs){
        store(parse_config_file(ifs, fileOptions), vm);
      }
    }
    notify(vm);
    if (vm.count("help"))
    {
      std::cout << generalOptions << '\n';
    }
  }
  catch (const bpo::error &ex)
  {
    std::cerr << ex.what() << '\n';
  }
  return vm;
} // end of handle_options

//Initialization of the first ID
unsigned int decrypt::utils::GeneCopy::m_next_available_id = 0;

//Class for wrapping sqlite3pp code
class database_type
{
public:
  static auto already_exists(std::string const& filename)
  {
    std::string message = "Unable to create database: file " + filename +
    " already exists.\nRemove existing file or change the database name in the configuration file.";
    return std::runtime_error(message);
  }
  database_type(std::string const& filename){
    if(std::filesystem::exists(filename)) throw(already_exists(filename));
    this->m_database = sqlite3pp::database(filename.c_str());
    create_results_table();
  }

  void insert_results(std::string const& newicks)
  {
    sqlite3pp::command cmd(
      this->m_database,
      "INSERT INTO results (newicks) VALUES (?)"
    );
    cmd.binder() << newicks;
    cmd.execute();
  }

private:
  sqlite3pp::database m_database;

  void create_results_table()
  {
    sqlite3pp::command cmd(
      this->m_database,
      "CREATE TABLE IF NOT EXISTS results(id INTEGER PRIMARY KEY AUTOINCREMENT UNIQUE, newicks TEXT)"
    );
    cmd.execute();
  }

};

class SimulationContext
{
public:
  SimulationContext(bpo::variables_map const& opts, std::mt19937& gen):
  vm(opts),
  m_database(build_database()),
  m_landscape(build_landscape()),
  m_sample(build_sample()),
  m_core(build_simulation_core()),
  m_reproduction_expr(build_reproduction_function(gen)),
  m_dispersal_kernel(build_dispersal_kernel())
  {
    show_reprojected_sample();
  };

  void run(std::mt19937& gen)
  {
    expand_demography(gen);
    maybe_save_demography();
    simulate_coalescence(gen);
    save_genealogies();
  }

private:
  using time_type = int;
  using landscape_type = geo::DiscreteLandscape<std::string,time_type>;
  using coord_type = landscape_type::coord_type;
  using loader_type = genet::Loader<coord_type, quetzal::genetics::microsatellite>;
  using sample_type = loader_type::return_type;
  using demographic_policy = demography::strategy::mass_based;
  using coal_policy = coal::policies::distance_to_parent_leaf_name<coord_type, time_type>;
  using core_type = sim::SpatiallyExplicit<coord_type, time_type, demographic_policy, coal_policy>;
  using options_type = bpo::variables_map;
  // eww, but it works
  using dispersal_type = demography::strategy::mass_based::light_neighboring_migration
  <
    coord_type,
    std::function<double(coord_type)>,
    std::function<std::vector<coord_type>(coord_type)>
  >;
  using reproduction_type = std::function<unsigned int(std::mt19937&, coord_type, time_type)>;
  database_type m_database;
  landscape_type m_landscape;
  sample_type m_sample;
  core_type m_core;
  reproduction_type m_reproduction_expr;
  time_type m_t_0;
  time_type m_sample_time;
  dispersal_type m_dispersal_kernel;
  std::string m_newicks;
  bpo::variables_map vm;

  database_type build_database()
  {
    if(vm.count("database"))
    {
      try
      {
        std::string filename = vm["database"].as<std::string>();
        return database_type(filename);
      }
      catch(const std::exception& e)
      {
        std::cout << "In SimulationContext, building database: " << e.what();
      }
    }
  }

  // TODO: pas sûr de la syntaxe du try/catch
  landscape_type build_landscape()
  {
    const std::string filename = vm["landscape"].as<std::string>();
    try
    {
      return landscape_type({{"suitability", filename}}, {time_type(0)});
    }
    catch(const std::exception& e)
    {
      std::cout << "In SimulationContext, building landscape: " << e.what();
    }
  }

  // TODO haploid versus diploid. Plus, format changed: ID/coordinates/no genetics
  sample_type build_sample()
  {
    std::string datafile = vm["sample"].as<std::string>();
    quetzal::genetics::Loader<coord_type, quetzal::genetics::microsatellite> reader;
    auto sample = reader.read(datafile);
    sample.reproject(m_landscape);
    return sample;
  }

  void show_reprojected_sample()
  {
    std::cout << "Reprojected sample:\n\n" << m_sample << std::endl;
  }

  core_type build_simulation_core()
  {
    coord_type x_0(vm["lat_0"].as<double>(), vm["lon_0"].as<double>());
    x_0 = m_landscape.reproject_to_centroid(x_0);
    m_t_0 = 0;
    m_sample_time = vm["duration"].as<unsigned int>();
    unsigned int N_0 = vm["N_0"].as<unsigned int>();
    return core_type(x_0, m_t_0, N_0);
  }

  reproduction_type build_reproduction_function(std::mt19937 & gen)
  {
    using expr::literal_factory;
    using expr::use;

    // growth rate
    literal_factory<coord_type, time_type> lit;
    auto r = lit( vm["r"].as<double>() );

    // carrying capacity
    auto suitability = m_landscape["suitability"];
    unsigned int K_suit = vm["K_suit"].as<unsigned int>();
    unsigned int K_min = vm["K_min"].as<unsigned int>();
    unsigned int K_max = vm["K_max"].as<unsigned int>();
    double p_K = vm["p_K"].as<double>();
    auto K = [K_suit, K_min, K_max, p_K, &gen, suitability](coord_type const& x, time_type)
    {
      if( suitability(x,0) == 0)
      { //ocean cell
        return std::bernoulli_distribution(p_K)(gen) ? K_max : K_min;
      }else{
        return K_suit;
      }
    };

    // Retrieve population size reference to define a logistic growth process
    auto pop_sizes = m_core.pop_size_history();
    auto N = use( [pop_sizes](coord_type x, time_type t){ return pop_sizes(x,t);} );
    auto g = N * ( lit(1) + r ) / ( lit(1) + ( (r * N)/K ));
    auto reproduction = [g](auto& gen, coord_type const&x, time_type t){
      std::poisson_distribution<unsigned int> poisson(g(x,t));
      return poisson(gen);
    };
    return reproduction;
  }

  dispersal_type build_dispersal_kernel()
  {
    auto suitability = m_landscape["suitability"];
    std::function<double(coord_type)> friction = [&suitability](coord_type const& x){
      if(suitability(x,0) <= 0.1) {return 0.9;} //ocean cell
      else return 1 - suitability(x, 0);
    };
    double emigrant_rate = vm["emigrant_rate"].as<double>();
    auto env_ref = std::cref(m_landscape);
    std::function<std::vector<coord_type>(coord_type)> get_neighbors = decrypt::utils::make_neighboring_cells_functor(env_ref);
    return demographic_policy::make_light_neighboring_migration(coord_type(), emigrant_rate, friction, get_neighbors);
  }

  void expand_demography(std::mt19937 & gen)
  {
    m_core.expand_demography(m_sample_time, m_reproduction_expr, m_dispersal_kernel, gen);
  }

  void maybe_save_demography()
  {
    if(vm.count("demography_out"))
    {
      try
      {
        std::string filename = vm["demography_out"].as<std::string>();
        if(std::filesystem::exists(filename))
        {
          std::string message("Unable to save demography: file " +filename+ " already exists.");
          throw(std::runtime_error(message));
        }
        using expr::use;
        auto pop_sizes = m_core.pop_size_history();
        auto N = use( [pop_sizes](coord_type x, time_type t){ return pop_sizes(x,t);} );
        m_landscape.export_to_geotiff(N, m_t_0, m_sample_time, [&pop_sizes](time_type const& t){return pop_sizes.get().definition_space(t);}, filename);
      }
      catch(const std::exception& e)
      {
        std::cout << e.what();
      }
    }
  }

  void simulate_coalescence(std::mt19937 & gen)
  {
    using decrypt::utils::GeneCopy;
    std::vector<GeneCopy> v;
    for(auto const& it1 : m_sample.get_sampling_points())
    {
      for(unsigned int i = 0; i < m_sample.individuals_at(it1).size(); ++ i)
      {// TODO here for diploidi and names
        v.emplace_back(it1);
      }
    }
    auto get_name = [](auto const& ind, time_type){return std::to_string(ind.id());};
    auto get_position = [](auto const& ind, time_type){return ind.x();};
    std::string genealogies;
    unsigned int n_loci = vm["n_loci"].as<unsigned int>();
    for(unsigned int locus = 0; locus < n_loci ; ++locus)
    {
      genealogies.append(m_core.coalesce_to_mrca<>(v, m_sample_time, get_position, get_name, gen));
      genealogies.append("\n\n");
    }
    genealogies.pop_back();
    genealogies.pop_back();
    m_newicks = genealogies;
  }

  void save_genealogies()
  {
    m_database.insert_results(m_newicks);
  }
};

#endif
