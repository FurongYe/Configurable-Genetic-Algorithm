#ifndef _OPERATOR_COMMON_H_
#define _OPERATOR_COMMON_H_

#include "utils/common.hpp"

using namespace std;

static default_random_engine random_gen(1);
static normal_distribution<double> normal_dis(0, 1);
static uniform_real_distribution<double> uniform_dis(0.0, 1.0);

enum optimizationType
{
    MINIMIZATION = 0,
    MAXIMIZATION = 1
};

static optimizationType Opt = MAXIMIZATION;

namespace ioha
{

    namespace operators
    {

        static double normal_random()
        {
            return normal_dis(random_gen);
        }

        static double uniform_random()
        {
            return uniform_dis(random_gen);
        }

        /// \fn normal
        /// @brief realiztion of the N(n,p)
        /// @param p
        /// @param n
        static int normal(double mu, double sigma)
        {
            int l = static_cast<int>(normal_random() * sigma + mu + 0.5);
            return l;
        }

        /// \fn normal_condition
        /// @brief realiztion of the 0 < N(n,p) < n
        /// @param p
        /// @param n
        static int normal_condition(double mu, double sigma, int upperbound)
        {
            int l;
            do
            {
                l = normal_random() * sigma + mu;
            } while (l <= 0 || l >= upperbound);
            return static_cast<int>(l + 0.5);
        }

        /// \fn binomial_condition
        /// @brief realiztion of the Bin(n,p)_>0
        /// @param p
        /// @param n
        static int binomial_condition(const double p, const int n)
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

        /// \fn binomial
        /// @brief realiztion of the Bin(n,p)
        /// @param p
        /// @param n
        static int binomial(const double p, const int n)
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

        static std::vector<double> power_law_distribution(int N, double beta_f)
        {
            auto power_law_distribution = std::vector<double>(N + 1);
            double C;
            size_t i;

            C = 0.0;
            for (i = 0; i < static_cast<size_t>(N / 2); ++i)
            {
                C += pow(i + 1, -beta_f);
            }
            for (i = 0; i < N; ++i)
            {
                power_law_distribution[i + 1] = 1 / C * pow(i + 1, -beta_f);
            }

            power_law_distribution[0] = 0.0;
            return power_law_distribution;
        }

        /// \fn sample_n_from_m
        /// \brief Sampling n different indexes from length m
        static vector<size_t> sample_n_from_m(size_t n, size_t m)
        {
            vector<size_t> sampled_number;

            if (n == 0)
            {
                std::clog << "sampled zero number" << std::endl;
            }

            size_t randPos;
            sampled_number.reserve(n);

            if (n > m / 2)
            { /// If n is larger than m/2, we sample random indexes by reordering a permutation.
                vector<size_t> population;
                population.reserve(m);
                for (size_t i = 0; i < m; ++i)
                {
                    population.push_back(i);
                }

                int temp;
                for (size_t i = m - 1; i > 0; --i)
                {
                    randPos = static_cast<size_t>(floor(uniform_random() * (i + 1)));
                    temp = population[i];
                    population[i] = population[randPos];
                    population[randPos] = temp;
                    sampled_number.push_back(population[i]);
                    if (m - i - 1 == n - 1)
                    {
                        break;
                    }
                }
                if (n == m)
                {
                    sampled_number.push_back(population[0]);
                }
            }
            else
            { /// If n is smaller than m/2, we sample indexes repeatly until getting different values.
                bool resample = false;
                for (size_t i = 0; i != n; ++i)
                {
                    do
                    {
                        resample = false;
                        randPos = static_cast<size_t>(floor(uniform_random() * m));
                        for (size_t j = 0; j != i; ++j)
                        {
                            if (randPos == sampled_number[j])
                            {
                                resample = true;
                                break;
                            }
                        }
                    } while (resample);
                    sampled_number.push_back(randPos);
                }
            }
            return sampled_number;
        }

        static void sample_n_from_m(vector<size_t> &sampled_number, size_t n, size_t m)
        {
            if (sampled_number.size() != 0)
            {
                sampled_number.clear();
            }

            if (n == 0)
            {
                std::clog << "sampled zero number" << std::endl;
            }

            size_t randPos;
            sampled_number.reserve(n);

            if (n > m / 2)
            { /// If n is larger than m/2, we sample random indexes by reordering a permutation.
                vector<size_t> population;
                population.reserve(m);
                for (size_t i = 0; i < m; ++i)
                {
                    population.push_back(i);
                }

                int temp;
                for (size_t i = m - 1; i > 0; --i)
                {
                    randPos = static_cast<size_t>(floor(uniform_random() * (i + 1)));
                    temp = population[i];
                    population[i] = population[randPos];
                    population[randPos] = temp;
                    sampled_number.push_back(population[i]);
                    if (m - i - 1 == n - 1)
                    {
                        break;
                    }
                }
                if (n == m)
                {
                    sampled_number.push_back(population[0]);
                }
            }
            else
            { /// If n is smaller than m/2, we sample indexes repeatly until getting different values.
                bool resample = false;
                for (size_t i = 0; i != n; ++i)
                {
                    do
                    {
                        resample = false;
                        randPos = static_cast<size_t>(floor(uniform_random() * m));
                        for (size_t j = 0; j != i; ++j)
                        {
                            if (randPos == sampled_number[j])
                            {
                                resample = true;
                                break;
                            }
                        }
                    } while (resample);
                    sampled_number.push_back(randPos);
                }
            }
        };
    }
}

#endif