/// \file randomsearch.h
/// \brief Header file for class ghc.
///
/// Random search generates candidates solution randomly.
///
///
/// \author Furong Ye
/// \date 2020-11-23

#ifndef _RANDOMSEARCH_H_
#define _RANDOMSEARCH_H_

#include "operators/common.hpp"
#include "algorithms/base.hpp"

namespace ioha
{
  namespace alg
  {
    class RandomSearch : public ioha::IOHA
    {
    public:
      RandomSearch(){};

      ~RandomSearch() {}
      RandomSearch(const RandomSearch &) = delete;
      RandomSearch &operator=(const RandomSearch &) = delete;

      /// \fn DoRandomSearch()
      /// \brief Do function of the random search algorithm.
      ///
      /// The order of processing functions is:
      /// Loop{ random_sampling } -> ~Termination().
      /// The function is virtual, which allows users to implement their own algorithm.
      virtual void opt()
      {
        this->preparation();

        while (!this->termination())
        {
          ++this->generation_;

          for (int i = 0; i != this->get_dimension(); ++i)
          {
            if (operators::uniform_random() < 0.5)
            {
              this->solution_[i] = 1;
            }
            else
            {
              this->solution_[i] = 0;
            }
          }

          this->solution_fitness_ = this->evaluate(this->solution_);
        }
      }

      void preparation()
      {
        this->problem_->reset();

        this->solution_ = vector<int>(this->get_dimension());
        this->evaluation_ = 0;
        this->generation_ = 0;
        this->best_individual_ = vector<int>(this->get_dimension());
        if (Opt == optimizationType::MAXIMIZATION)
        {
          this->best_found_fitness_ = numeric_limits<double>::lowest();
        }
        else
        {
          this->best_found_fitness_ = numeric_limits<double>::max();
        }
      }

      void run(shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite)
      {
        for (auto &p : *suite)
        {
          this->problem_ = p;
          size_t i = 0;
          while (++i <= this->independent_runs_)
          {
            size_t i = 0;
            this->opt();
          }
          this->problem_->reset();
        }
      }

      void run(string folder_path, string folder_name, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite, int eval_budget, int independent_runs, unsigned rand_seed)
      {
        string algorithm_name = "random search";
        ioh::logger::Analyzer logger(
            {ioh::trigger::on_improvement}, // trigger when the objective value improves
            {},                             // no additional properties
            folder_path,                    // path to store data
            folder_name,                    // name of the folder in path, which will be newly created
            algorithm_name,                 // name of the algoritm
            algorithm_name,                 // additional info about the algorithm
            false                           // where to store x positions in the data files
            );
        this->problem_->attach_logger(logger);

        this->set_evaluation_budget(eval_budget);
        this->set_independent_runs(independent_runs);
        this->set_seed(rand_seed);

        this->run(suite);
      }

      void set_evaluation_budget(const size_t evaluation_budget)
      {
        this->evluation_budget_ = evaluation_budget;
      }

      void set_independent_runs(const size_t independent_runs)
      {
        this->independent_runs_ = independent_runs;
      }

      void set_solution(const vector<int> &solution)
      {
        this->solution_ = solution;
      }

      void set_solution_fitness(const double solution_fitness)
      {
        this->solution_fitness_ = solution_fitness;
      }

      void set_best_found_fitness(const double best_found_fitness)
      {
        this->best_found_fitness_ = best_found_fitness;
      }

      void set_best_individual(const vector<int> best_individual)
      {
        this->best_individual_ = best_individual;
      }

      void set_generation(const size_t generation)
      {
        this->generation_ = generation;
      }

      int get_mu() const
      {
        return 1;
      }

      int get_lambda() const
      {
        return 1;
      }

      int get_dimension() const
      {
        return this->problem_->meta_data().n_variables;
      }

      vector<int> get_solution()
      {
        return this->solution_;
      }

      double get_solution_fitness()
      {
        return this->solution_fitness_;
      }

      double get_best_found_fitness()
      {
        return this->best_found_fitness_;
      }

      vector<int> get_best_individual()
      {
        return this->best_individual_;
      }

      size_t get_generation()
      {
        return this->generation_;
      }

    private:
      vector<int> solution_;
      double solution_fitness_;
      double best_found_fitness_;
      vector<int> best_individual_;

      size_t evaluation_;       /// < evaluation times
      size_t generation_;       /// < number of iterations/generations
      size_t evluation_budget_; /// < budget for evaluations

      size_t independent_runs_; /// < number of independent runs.

      /// TODO: we assume the type of problem are integer only now.
      shared_ptr<ioh::problem::IntegerSingleObjective> problem_;
      shared_ptr<ioh::logger::Analyzer> csv_logger_;
    };
  }
}

#endif // _RANDOMSEARCH_H_