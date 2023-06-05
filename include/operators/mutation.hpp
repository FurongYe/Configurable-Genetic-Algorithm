#pragma once

/// \file mutation.h
/// \brief Header file for class Mutation.
///
/// It contains functions of mutation operators.
///
/// \author Furong Ye
/// \date 2023-03-23

#include "common.hpp"

#define DEFAULT_MUTATION_STRENGTH_ 1
#define DEFAULT_MUTATION_RATE_ 0.01
#define DEFAULT_NORMALIZED_MUTATION_STRENGTH_MEAN_ 1
#define DEFAULT_NORMALIZED_MUTATION_STRENGTH_SD_ 0.01
#define DEFAULT_FAST_MUTATION_BETA_ 1.5

namespace ioha
{

  namespace operators
  {
      enum mutation_operator
      {
        STATICSAMPLE = 1,   /// < mutation strength is a fixed value.
        BINOMIALSAMPLE = 2, /// < sampling mutation strength from a binomial distribution.
        NORMALSAMPLE = 3,   /// < sampling mutation strength from a normal distribution.
        POWERLAWSAMPLE = 4  /// < sampling mutation strength from a power-law distribution.
      };

      class Mutation
      {
      public:
        Mutation() : mutation_operator_(2),
                     l_(DEFAULT_MUTATION_STRENGTH_),
                     mutation_rate_(DEFAULT_MUTATION_RATE_),
                     r_n_(DEFAULT_NORMALIZED_MUTATION_STRENGTH_MEAN_),
                     sigma_n_(DEFAULT_NORMALIZED_MUTATION_STRENGTH_SD_),
                     beta_f_(DEFAULT_FAST_MUTATION_BETA_) {}

        ~Mutation() {}
        Mutation(const Mutation &) = delete;
        Mutation &operator=(const Mutation &) = delete;

        int do_mutation(vector<int> &y)
        {
          int mutation_strength = this->sample_l(y.size());
          this->flip(y, mutation_strength);
          return mutation_strength;
        }
        /// \fn flip
        /// @brief flip l bits of y
        /// @param y
        /// @param l
        void flip(vector<int> &y, const int l)
        {
          size_t n = y.size();
          m_flipped_index = operators::sample_n_from_m(static_cast<size_t>(l), n);
          for (int i = 0; i != l; ++i)
          {
            y[m_flipped_index[i]] = (y[m_flipped_index[i]] + 1) % 2;
          }
        }

        static int sample_conditional_binomial(const double p, const int n)
        {
          int l = 0;
          while (l == 0)
          {
            for (int i = 0; i != n; ++i)
            {
              if (uniform_random() < p)
              {
                ++l;
              }
            }
          }
          return l;
        }

        static int sample_binomial(const double p, const int n)
        {
          int l = 0;
          for (int i = 0; i != n; ++i)
          {
            if (uniform_random() < p)
            {
              ++l;
            }
          }
          return l;
        }

        static int sample_normal(double mu, double sigma)
        {
          int l = static_cast<int>(normal_random() * sigma + mu + 0.5);
          return l;
        }

        static int sample_conditional_normal(double mu, double sigma, int upperbound)
        {
          int l;
          do
          {
            l = normal_random() * sigma + mu;
          } while (l <= 0 || l >= upperbound);
          return static_cast<int>(l + 0.5);
        }

        static int sample_logNormal(double mu, double sigma)
        {
          double l = normal_random() * sigma + mu;
          return static_cast<int>(exp(l) + 0.5);
        }

        static int sample_conditional_logNormal(double mu, double sigma, int upperbound)
        {
          double l;
          do
          {
            l = normal_random() * sigma + mu;
          } while (static_cast<int>(exp(l)) <= 0 || static_cast<int>(exp(l)) >= upperbound);
          return static_cast<int>(exp(l) + 0.5);
        }

        static int sample_from_distribution(const std::vector<double> &distribution)
        {
          double p = 0.0;
          double r = uniform_random();
          int j = 0;
          while (j < distribution.size())
          {
            if (r > p && r < p + distribution[j])
            {
              break;
            }
            p += distribution[j];
            ++j;
          }
          return j;
        }

        int sample_l(const int N)
        {
          int mutation_strength = 0;
          switch (this->mutation_operator_)
          {
          case 1:
            mutation_strength = this->l_;
            break;
          case 2:
            mutation_strength = this->sample_conditional_binomial(this->mutation_rate_, N);
            break;
          case 3:
            mutation_strength = this->sample_conditional_normal(this->r_n_, this->sigma_n_, floor(N / 2.0));
            break;
          case 4:
            mutation_strength = this->sample_from_distribution(this->power_law_distribution_);
            break;
          default:
            cerr << "unknown mutation operator" << endl;
            assert(false);
          }
          return mutation_strength;
        }

        void set_mutation_operator(const int m)
        {
          assert((m <= 4) && (m >= 1));
          this->mutation_operator_ = m;
        }

        void set_mutation_operator(string m)
        {
          transform(m.begin(), m.end(), m.begin(), ::toupper);
          if (m == "STATICSAMPLE")
          {
            this->set_mutation_operator(STATICSAMPLE);
          }
          else if (m == "BINOMIALSAMPLE")
          {
            this->set_mutation_operator(BINOMIALSAMPLE);
          }
          else if (m == "NORMALSAMPLE")
          {
            this->set_mutation_operator(NORMALSAMPLE);
          }
          else if (m == "POWERLAWSAMPLE")
          {
            this->set_mutation_operator(POWERLAWSAMPLE);
          }
          else
          {
            cerr << "invalid name for mutation operator" << endl;
            assert(false);
          }
        }

        void set_l(const int l)
        {
          this->l_ = l;
        }

        void set_mutation_rate(const double mutation_rate)
        {
          this->mutation_rate_ = mutation_rate;
        }

        void set_r_n(const double r_n)
        {
          this->r_n_ = r_n;
        }

        void set_sigma_n(const double sigma_n)
        {
          this->sigma_n_ = sigma_n;
        }

        void set_beta_f(const double beta_f)
        {
          this->beta_f_ = beta_f;
        }

        int get_mutation_operator() const
        {
          return this->mutation_operator_;
        }

        int get_l() const
        {
          return this->l_;
        }

        double get_mutation_rate() const
        {
          return this->mutation_rate_;
        }

        double get_r_n() const
        {
          return this->r_n_;
        }

        double get_sigma_n() const
        {
          return this->sigma_n_;
        }
        double get_beta_f() const
        {
          return this->beta_f_;
        }

        vector<size_t> m_flipped_index;
      
      protected:

        std::vector<double> power_law_distribution_;
      private:

        int mutation_operator_; /// < a flag for operator to be used for mutation
        int l_;                 /// < a static mutation strength
        double mutation_rate_;  /// < mutation rate for standard bit mutation
        double r_n_;            /// < mean value for normalized bit mutation
        double sigma_n_;        /// < sd for normalized bit mutation
        double beta_f_;         /// < beta for fast mutation
      };
  }

} // namespace ioha
