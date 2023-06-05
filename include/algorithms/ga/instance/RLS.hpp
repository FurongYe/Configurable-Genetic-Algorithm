#pragma once
#include "algorithms/ga/geneticAlgorithm.hpp"
namespace ioha
{
  namespace alg
  {
    class RLS : public GeneticAlgorithm
    {

    public:
      RLS()
      {
        this->set_mu(1);
        this->set_lambda(1);
        this->set_crossover_mutation_r("IND");
        this->set_crossover_probability(0);
        this->set_mutation_operator("STATICSAMPLE");
        this->set_l(1);
        this->set_selection_operator("BESTPLUS");
      }

      void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, int eval_budget, int gene_budget, int independent_runs, unsigned rand_seed)
      {
        string algorithm_name = "RLS";
        ioh::logger::Analyzer logger(
            {ioh::trigger::on_improvement}, // trigger when the objective value improves
            {},                             // no additional properties
            folder_path,                    // path to store data
            folder_name,                    // name of the folder in path, which will be newly created
            algorithm_name,                 // name of the algoritm
            algorithm_name,                 // additional info about the algorithm
            false                           // where to store x positions in the data files
            );
        suite->attach_logger(logger);
        this->set_evaluation_budget(eval_budget);
        this->set_generation_budget(gene_budget);
        this->set_independent_runs(independent_runs);
        this->set_seed(rand_seed);

        this->run_N(suite);
      }
    };
  }
}