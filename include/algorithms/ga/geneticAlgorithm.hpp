#pragma once

/// \file geneticAlgorithm.h
/// \brief Header file for class GeneticAlgorithm.
///
/// A configurable genetic algorithm.
///
/// \author Furong Ye
/// \date 2020-07-02

#include "utils/common.hpp"
#include "operators/crossover.hpp"
#include "operators/mutation.hpp"
#include "operators/selection.hpp"

#include "algorithms/base.hpp"

#define DEFAULT_CROSSOVER_PROBABLITY_ 0
#define DEFAULT_CROSSOVER_MUTATION_RELATION_ 0

#define DEFAULT_ECDF_TARGET_LENGTH 100
#define DEFAULT_ECDF_BUDGET_LENGTH 100

namespace ioha
{
  namespace alg
  {

    /// Definition of relation between mutation and crossover.
    enum crossover_mutation_relation
    {
      IND = 1, /// < Always doing mutation, and doing crossover with probability p_c.
      OR = 0   /// < Doing crossover probability p_c, otherwise doing mutation.
    };

    class GeneticAlgorithm : public ioha::IOHA, public ioha::operators::Crossover, public ioha::operators::Mutation, public ioha::operators::Selection
    {
    public:
      GeneticAlgorithm() : crossover_probability_(DEFAULT_CROSSOVER_PROBABLITY_),
                           crossover_mutation_r_(DEFAULT_CROSSOVER_MUTATION_RELATION_)
      {
      }

      ~GeneticAlgorithm() {}
      GeneticAlgorithm(const GeneticAlgorithm &) = delete;
      GeneticAlgorithm &operator=(const GeneticAlgorithm &) = delete;

      /// \fn SetAllParameters()
      /// \brief Set all parameters of the genetic algorithm.
      /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
      void set_parameters(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para)
      {
        assert(integer_para.size() == 4);
        this->set_mu(integer_para[0]);
        this->set_lambda(integer_para[1]);
        this->set_l(integer_para[2]);
        this->set_tournament_k(integer_para[3]);

        assert(continuous_para.size() == 6);
        this->set_crossover_probability(continuous_para[0]);
        this->set_p_u(continuous_para[1]);
        this->set_mutation_rate(continuous_para[2]);
        this->set_r_n(continuous_para[3]);
        this->set_sigma_n(continuous_para[4]);
        this->set_beta_f(continuous_para[5]);

        assert(category_para.size() == 4);
        this->set_crossover_mutation_r(category_para[0]);
        this->set_crossover_operator(category_para[1]);
        this->set_mutation_operator(category_para[2]);
        this->set_selection_operator(category_para[3]);
      }

      /// \fn run()
      /// \brief Run genetic algorithm once with given parameters.
      /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}

      virtual void run(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para)
      {
        this->set_parameters(integer_para, continuous_para, category_para);
        this->opt();
      }

      /// \fn run()
      /// \brief Run genetic algorithm once with given parameters, on an IOHexperimenter suite.
      /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
      /// \param suite a shared pointer of IOHprofiler_suite

      virtual void run(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite)
      {
        this->set_parameters(integer_para, continuous_para, category_para);

        for (const auto &p : *suite)
        {
          this->problem_ = p;
          this->opt();
        }
      }

      /// \fn run_N()
      /// \brief Run genetic algorithm `independent_runs_` times with given parameters.
      /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}

      void run_N(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para)
      {
        this->set_parameters(integer_para, continuous_para, category_para);
        size_t r = 0;
        this->ecdf_sum_ = 0;
        while (r < this->independent_runs_)
        {

          this->opt();
          r++;
        }

        this->ecdf_sum_ = this->ecdf_sum_ / this->independent_runs_;
        this->ecdf_ratio_ = (double)this->ecdf_sum_ / (double)this->ecdf_budget_width_ / (double)this->ecdf_target_width_ / this->independent_runs_;
      }

      /// \fn run_N()
      /// \brief Run genetic algorithm `independent_runs_` times with given parameters, on an IOHexperimenter suite.
      /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
      /// \param suite a shared pointer of IOHprofiler_suite
      virtual void run_N(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite)
      {
        this->set_parameters(integer_para, continuous_para, category_para);
        for (const auto &p : *suite)
        {
          this->problem_ = p;
          this->run_N(integer_para, continuous_para, category_para);
        }
      }

      /// \fn run_N()
      /// \brief Run genetic algorithm `independent_runs_` times  on an IOHexperimenter suite.
      /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
      /// \param suite a shared pointer of IOHprofiler_suite

      virtual void run_N(shared_ptr<ioh::suite::Suite<ioh::problem::IntegerSingleObjective>> suite)
      {
        for (const auto &p : *suite)
        {
          this->assign_problem(p);
          size_t r = 0;
          this->ecdf_sum_ = 0;
          while (r < this->independent_runs_)
          {
            this->opt();
            r++;
          }
          this->problem_->reset();
          this->ecdf_sum_ = this->ecdf_sum_ / this->independent_runs_;
          this->ecdf_ratio_ = (double)this->ecdf_sum_ / (double)this->ecdf_budget_width_ / (double)this->ecdf_target_width_ / this->independent_runs_;
        }
      }

      virtual void opt()
      {
        double rand;
        this->preparation();

        this->initialization();
        while (!this->termination())
        {
          ++this->generation_;

          this->offspring_population_.clear();
          this->offspring_fitness_.clear();
          for (size_t i = 0; i < this->lambda_; ++i)
          {
            this->select_two_parents();
            this->offspring_population_.push_back(this->parents_population_[this->selected_parents_[0]]);

            rand = operators::uniform_random();
            this->c_flipped_index.clear();
            this->m_flipped_index.clear();
            if (rand < this->crossover_probability_)
            {
              this->do_crossover(this->offspring_population_[i], this->parents_population_[this->selected_parents_[0]], this->parents_population_[this->selected_parents_[1]]);
            }

            if (this->crossover_mutation_r_)
            {
              this->do_mutation(this->offspring_population_[i]);
            }
            else if (rand >= this->crossover_probability_)
            {
              this->do_mutation(this->offspring_population_[i]);
            }

            if (this->c_flipped_index == this->m_flipped_index)
            { /// If the flipping indexes of crossover and mutation are identical, the individual remains the same.
              this->offspring_fitness_.push_back(this->parents_fitness_[this->selected_parents_[0]]);
            }
            else if (rand < this->crossover_probability_ && this->offspring_population_[i] == this->parents_population_[this->selected_parents_[1]])
            { /// If the offspring is identical with the second parent.
              /// TODO: Do something to save time for this comparison.
              this->offspring_fitness_.push_back(this->parents_fitness_[this->selected_parents_[1]]);
            }
            else
            { /// otherwise evaluate.
              this->offspring_fitness_.push_back(this->evaluate(this->offspring_population_[i]));
            }

            if (this->termination())
              break;
          }

          if (this->termination())
            break;

          this->do_selection(this->parents_population_, this->parents_fitness_, this->offspring_population_, this->offspring_fitness_);
          this->adaptation();
        }
      }

      virtual void initialization()
      {
        int n = this->problem_->meta_data().n_variables;
        for (int i = 0; i != this->mu_; ++i)
        {
          vector<int> tmp(n, 0);
          for (int j = 0; j != n; ++j)
          {
            if (operators::uniform_random() < 0.5)
            {
              tmp[j] = 1;
            }
          }
          this->parents_population_.push_back(tmp);
          this->parents_fitness_.push_back(this->evaluate(tmp));
        }
      }

      void preparation()
      {
        this->problem_->reset();

        this->selected_parents_ = vector<size_t>(2);
        this->optimum_ = this->problem_->optimum().y;
        this->power_law_distribution_ = ioha::operators::power_law_distribution(this->problem_->meta_data().n_variables, 1.5);
        this->parents_fitness_.clear();
        this->parents_population_.clear();
        this->offspring_fitness_.clear();
        this->offspring_population_.clear();
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
        this->hitting_flag_ = false;
      }

      void select_two_parents()
      {
        this->selected_parents_[0] = static_cast<size_t>(floor(operators::uniform_random() * this->mu_));
        if (this->mu_ >= 2)
        {
          do
          {
            this->selected_parents_[1] = static_cast<size_t>(floor(operators::uniform_random() * this->mu_));
          } while (this->selected_parents_[0] == this->selected_parents_[1]);
        }
      }

      double evaluate(vector<int> &x)
      {
        double result;
        result = (*this->problem_)(x);
        // if (this->csv_logger_ != nullptr) {
        //   this->csv_logger_->do_log(this->problem_->loggerInfo()); /// TODO: we assume only using PBO suite now.
        // }
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

        // This is only used to calculate ERT
        if (this->ERT_flag_)
        {
          if (Opt == optimizationType::MAXIMIZATION)
          {
            if (!this->hitting_flag_ && this->best_found_fitness_ >= this->hiting_target_)
            {
              this->hitting_time_.push_back(this->evaluation_);
              this->hitting_flag_ = true;
            }
          }
          else
          {
            if (!this->hitting_flag_ && this->best_found_fitness_ <= this->hiting_target_)
            {
              this->hitting_time_.push_back(this->evaluation_);
              this->hitting_flag_ = true;
            }
          }
        }
        return result;
      }

      void set_crossover_probability(const double crossover_probability)
      {
        this->crossover_probability_ = crossover_probability;
      }

      void set_crossover_mutation_r(const int crossover_mutation_r)
      {
        this->crossover_mutation_r_ = crossover_mutation_r;
      }

      void set_crossover_mutation_r(string crossover_mutation_r)
      {
        transform(crossover_mutation_r.begin(), crossover_mutation_r.end(), crossover_mutation_r.begin(), ::toupper);
        if (crossover_mutation_r == "IND")
        {
          this->set_crossover_mutation_r(IND);
        }
        else if (crossover_mutation_r == "OR")
        {
          this->set_crossover_mutation_r(OR);
        }
        else
        {
          cerr << "invalid value for set_crossover_mutation_r";
          assert(false);
        }
      }

      void set_independent_runs(const size_t independent_runs)
      {
        this->independent_runs_ = independent_runs;
      }

      void set_parents_population(const vector<vector<int>> parents_population)
      {
        this->parents_population_ = parents_population;
      }

      void set_parents_fitness(const vector<double> &parents_fitness)
      {
        this->parents_fitness_ = parents_fitness;
      }

      void set_parents_population(const vector<int> &parent, const size_t index)
      {
        assert(this->parents_population_.size() > index);
        this->parents_population_[index] = parent;
      }

      void set_parents_fitness(const double fitness, const size_t index)
      {
        assert(this->parents_fitness_.size() > index);
        this->parents_fitness_[index] = fitness;
      }

      void add_parents_population(const vector<int> &parent)
      {
        this->parents_population_.push_back(parent);
      }
      void clear_parents_population()
      {
        this->parents_population_.clear();
      }

      void add_parents_fitness(const double fitness)
      {
        this->parents_fitness_.push_back(fitness);
      }

      void clear_parents_fitness()
      {
        this->parents_fitness_.clear();
      }

      void set_offspring_population(const vector<vector<int>> &offspring_population)
      {
        this->offspring_population_ = offspring_population;
      }

      void set_offspring_fitness(const vector<double> &offspring_fitness)
      {
        this->offspring_fitness_ = offspring_fitness;
      }

      void set_offspring_population(const vector<int> &offspring, const size_t index)
      {
        assert(this->offspring_population_.size() > index);
        this->offspring_population_[index] = offspring;
      }

      void set_offspring_fitness(const double fitness, const size_t index)
      {
        assert(this->offspring_fitness_.size() > index);
        this->offspring_fitness_[index] = fitness;
      }

      void add_offspring_population(const vector<int> &parent)
      {
        this->offspring_population_.push_back(parent);
      }
      void clear_offspring_population()
      {
        this->offspring_population_.clear();
      }

      void add_offspring_fitness(const double fitness)
      {
        this->offspring_fitness_.push_back(fitness);
      }

      void clear_offspring_fitness()
      {
        this->offspring_fitness_.clear();
      }

      void set_best_found_fitness(const double best_found_fitness)
      {
        this->best_found_fitness_ = best_found_fitness;
      }

      void set_best_individual(const vector<int> best_individual)
      {
        this->best_individual_ = best_individual;
      }

      double get_crossover_probability() const
      {
        return this->crossover_probability_;
      }

      int get_crossover_mutation_r() const
      {
        return this->crossover_mutation_r_;
      }

      vector<vector<int>> get_parents_population() const
      {
        return this->parents_population_;
      }

      vector<double> get_parents_fitness() const
      {
        return this->parents_fitness_;
      }

      vector<vector<int>> get_offspring_population() const
      {
        return this->offspring_population_;
      }

      vector<double> get_offspring_fitness() const
      {
        return this->offspring_fitness_;
      }

      double get_best_found_fitness() const
      {
        return this->best_found_fitness_;
      }

      vector<int> get_best_individual() const
      {
        return this->best_individual_;
      }

      size_t get_independent_runs() const
      {
        return this->independent_runs_;
      }

      //  /// \fn EstimateLinearECDF()
      //  /// \brief Calculating the area under ECDF (partition budgets and targets with linear scale) of the genetic algorithm with given parameters, on an IOHexperimenter suite with given targets.
      //  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      //  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      //  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
      //  /// \param suite a shared pointer of IOHprofiler_suite
      //  /// \param targets a vector of targets of problems of the suite.
      //  double EstimateLinearECDF(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite, vector<double> targets);
      //
      //  /// \fn EstimateLogECDF()
      //  /// \brief Calculating the area under ECDF (partition budgets and targets with log scale) of the genetic algorithm with given parameters, on an IOHexperimenter suite with given targets.
      //  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      //  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      //  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
      //  /// \param suite a shared pointer of IOHprofiler_suite
      //  /// \param targets a vector of targets of problems of the suite.
      //  double EstimateLogECDF(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite, vector<double> targets);
      //
      //  /// \fn EstimateLogERT()
      //  /// \brief Calculating ERT values of the genetic algorithm with given parameters, on an IOHexperimenter suite with given targets and budgets.
      //  /// \param integer_para {mu, lambda, the number of flipping bits for local search, the size of tournament selection}
      //  /// \param continuous_para {crossover probablity, probablity replacing bit by the bit of the other individual of uniform crossover,  mutation rate, mean of normalized bit mutation, standard deviation of normalized bit mutation, beta value of fast mutation}
      //  /// \param category_para {relation between crossover and mutation, crossover operator, mutation operator, selection operator}
      //  /// \param suite a shared pointer of IOHprofiler_suite
      //  /// \param targets a vector of targets of problems of the suite.
      //  /// \param budgets a vector of budgets of problems.
      //  vector<double> EstimateERT(const vector<int> &integer_para, const vector<double> &continuous_para, const vector<string> &category_para, shared_ptr<IOHprofiler_suite<int> > suite, const vector<double> &targets, const vector<size_t> &budgets);

    private:
     
      double crossover_probability_; /// < probability to do crossover
      int crossover_mutation_r_;     /// < a flag for correlation between crossover and mutation

      vector<vector<int>> parents_population_;
      vector<double> parents_fitness_;
      vector<vector<int>> offspring_population_;
      vector<double> offspring_fitness_;

      size_t independent_runs_; /// < number of independent runs.

      vector<size_t> selected_parents_;

      /// For calculating ECDF
      size_t ecdf_sum_;
      double ecdf_ratio_;
      size_t ecdf_budget_width_;
      size_t ecdf_target_width_;

      /// For calculating ERT
      vector<double> hitting_time_;
      double hiting_target_;
      bool hitting_flag_;
      bool ERT_flag_;
    };
  } // namespace algo

} // namespace modularGA
