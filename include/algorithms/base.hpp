#pragma once

/// \file base.hpp
/// \brief The base class for IOH algorithms.
///
/// \author Furong Ye
/// \date 2020-07-02

#include "utils/common.hpp"
#include "operators/crossover.hpp"
#include "operators/mutation.hpp"
#include "operators/selection.hpp"
#include "base.hpp"

#define DEFAULT_MU_ 1
#define DEFAULT_LAMBDA_ 1

#define DEFAULT_EVALUATION_BUDGET_ 10000
#define DEFAULT_GENERATION_BUDGET_ 10000

namespace ioha
{

    class IOHA
    {
    public:
        IOHA() : mu_(DEFAULT_MU_),
                 lambda_(DEFAULT_LAMBDA_),
                 evaluation_(0),
                 generation_(0),
                 evluation_budget_(DEFAULT_EVALUATION_BUDGET_),
                 generation_budget_(DEFAULT_GENERATION_BUDGET_),
                 optimum_(numeric_limits<double>::max()) /// < TODO: now we assume doing maximization.
        {
        }

        ~IOHA() {}

        /// \fn virtual adaptation()
        /// \brief A virtual function for adaptive methods.
        ///
        /// If you are about to implement an algorithm with adaptive paramters, you should implement the adaptive method in this function. This function will be revoked at the end of each iteration/generation.
        virtual void adaptation(){};

        /// \fn virtual initialization()
        /// \brief A virtual function for initialization.
        virtual void initialization(){};

        /// \fn virtual sampling()
        /// \brief A virtual function for Sampling new offspring
        virtual void sampling(){};

        /// \fn virtual selection()
        /// \brief A virtual function for form new parents
        virtual void selection(){};

        /// \fn virtual termination()
        /// \brief Terminate condition of the genetic algorithm.
        ///
        /// You can set up terminate condition in this function. By default, the algorithm terminates when the best fitness value is found or the budget is used out.
        virtual bool termination()
        {
            if (!this->find_optimum() && this->problem_->state().evaluations < this->evluation_budget_ && this->generation_ <= this->generation_budget_)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        /// \fn opt()
        /// \brief the loop of optimization
        virtual void opt()
        {
            this->initialization();
            while (this->termination())
            {
                this->sampling();
                this->selection();
                this->adaptation();
            }
        }

        void assign_problem(shared_ptr<ioh::problem::IntegerSingleObjective> problem_ptr)
        {
            this->problem_ = problem_ptr;
        }

        virtual double evaluate(vector<int> &x)
        {
            double result;
            result = (*this->problem_)(x);
            this->evaluation_++;
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

        void set_seed(unsigned seed)
        {
            random_gen.seed(seed);
        }

        void set_mu(const int mu)
        {
            this->mu_ = mu;
        }

        void set_lambda(const int lambda)
        {
            this->lambda_ = lambda;
        }

        void set_evaluation_budget(const size_t evaluation_budget)
        {
            this->evluation_budget_ = evaluation_budget;
        }

        void set_generation_budget(const size_t generation_budget)
        {
            this->generation_budget_ = generation_budget;
        }

        void set_generation(const size_t generation)
        {
            this->generation_ = generation;
        }

        int get_dimension() const
        {
            return this->problem_->meta_data().n_variables;
        }

        void update_generation()
        {
            this->generation_++;
        }

        int get_mu() const
        {
            return this->mu_;
        }

        int get_lambda() const
        {
            return this->lambda_;
        }

        size_t get_evaluation_budget()
        {
            return this->evluation_budget_;
        }

        size_t get_generation_budget()
        {
            return this->generation_budget_;
        }

        size_t get_generation() const
        {
            return this->generation_;
        }
        int get_evaluations()
        {
            return this->problem_->state().evaluations;
        }

        bool find_optimum()
        {
            return this->problem_->state().optimum_found;
        }

    protected:
        int mu_;     /// < parents population size
        int lambda_; /// < offspring population size

        double best_found_fitness_;
        vector<int> best_individual_;

        size_t evaluation_;
        size_t generation_;                   /// < number of iterations/generations
        size_t evluation_budget_;             /// < budget for evaluations
        size_t generation_budget_ = SIZE_MAX; /// < budget for generations

        /// TODO: we assume the type of problem are integer only now.
        std::shared_ptr<ioh::problem::IntegerSingleObjective> problem_;
        double optimum_;
    };

} // namespace modularGA
