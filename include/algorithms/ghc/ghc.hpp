/// \file ghc.h
/// \brief Header file for class ghc.
///
/// (1+1) greedy hill climber algorithm goes throught the string from left to right,
/// flipping exactly one bit per iteration, and accepting the offspring if it
/// is as least as good as its parent.
///
///
/// \author Furong Ye
/// \date 2020-11-23

#ifndef _GHC_H_
#define _GHC_H_

#include "operators/common.hpp"
#include "algorithms/base.hpp"
namespace ioha
{
  namespace alg{
    class GreedyHillClimber : ioha::IOHA
    {
    public:
      GreedyHillClimber(){};

      ~GreedyHillClimber() {}
      GreedyHillClimber(const GreedyHillClimber &) = delete;
      GreedyHillClimber &operator=(const GreedyHillClimber &) = delete;

    
      /// \fn DoGreedyHillClimber()
      /// \brief Do function of the greedy hill climber algorithm.
      ///
      /// The order of processing functions is: Initialization() ->
      /// Loop{ Mutation() -> Selection()} -> ~Termination().
      /// The function is virtual, which allows users to implement their own algorithm.
      virtual void opt()
      {
        this->preparation();

        this->initialization();
        size_t flip_index = 0;
        while (!this->termination())
        {
          ++this->generation_;

          this->offspring_.clear();
          this->offspring_ = this->parent_;

          // Filp the bit at flip_index.
          this->offspring_[flip_index] = (this->offspring_[flip_index] + 1) % 2;
          if (++flip_index >= this->get_dimension())
          {
            flip_index = 0;
          }

          this->offspring_fitness_ = this->evaluate(this->offspring_);

          if (this->offspring_fitness_ >= this->parent_fitness_)
          {
            this->parent_ = this->offspring_;
            this->parent_fitness_ = this->offspring_fitness_;
          }
        }
      }
      void initialization()
      {
        int n = this->problem_->meta_data().n_variables;

        this->parent_ = vector<int>(n, 0);
        for (int i = 0; i != n; ++i)
        {
          if (operators::uniform_random() < 0.5)
          {
            this->parent_[i] = 1;
          }
        }
        this->parent_fitness_ = this->evaluate(this->parent_);
      }

      void preparation()
      {
        this->problem_->reset();

        this->parent_.clear();
        this->offspring_.clear();
        this->evaluation_ = 0;
        this->generation_ = 0;
        this->best_individual_ = vector<int>(this->problem_->meta_data().n_variables);
        if (Opt == optimizationType::MAXIMIZATION)
        {
          this->best_found_fitness_ = numeric_limits<double>::lowest();
        }
        else
        {
          this->best_found_fitness_ = numeric_limits<double>::max();
        }
      }

      double evaluate(vector<int> &x)
      {
        double result;
        result = (*this->problem_)(x);

        ++this->evaluation_;

        if (Opt == optimizationType::MAXIMIZATION)
        {
          if (result > this->best_found_fitness_)
          {
            this->best_found_fitness_ = result;
            this->best_individual_ = x;
          }
        }
        else
        {
          if (result < this->best_found_fitness_)
          {
            this->best_found_fitness_ = result;
            this->best_individual_ = x;
          }
        }
        return result;
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
        string algorithm_name = "gHC";
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

      void set_parent(const vector<int> &parent)
      {
        this->parent_ = parent;
      }

      void set_parent_fitness(const double parent_fitness)
      {
        this->parent_fitness_ = parent_fitness;
      }

      void set_offspring(const vector<int> &offspring)
      {
        this->offspring_ = offspring;
      }

      void set_offspring_fitness(const double offspring_fitness)
      {
        this->offspring_fitness_ = offspring_fitness;
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

      vector<int> get_parent()
      {
        return this->parent_;
      }

      double get_parent_fitness()
      {
        return this->parent_fitness_;
      }

      vector<int> get_offspring()
      {
        return this->offspring_;
      }

      double get_offspring_fitness()
      {
        return this->offspring_fitness_;
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
      vector<int> parent_;
      double parent_fitness_;
      vector<int> offspring_;
      double offspring_fitness_;
      double best_found_fitness_;
      vector<int> best_individual_;

      size_t evaluation_;       /// < evaluation times
      size_t generation_;       /// < number of iterations/generations
      size_t evluation_budget_; /// < budget for evaluations

      size_t independent_runs_; /// < number of independent runs.

      /// TODO: we assume the type of problem are integer only now.
      shared_ptr<ioh::problem::IntegerSingleObjective> problem_;
    };
  }
}

#endif // _GHC_H_