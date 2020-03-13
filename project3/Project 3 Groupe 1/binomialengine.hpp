/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003, 2004 Ferdinando Ametrano
 Copyright (C) 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2007 StatPro Italia srl
 Copyright (C) 2007 Affine Group Limited

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

    /*\file binomialengine.hpp
    \brief Binomial option engine*/

#ifndef binomial_engine_hpp
#define binomial_engine_hpp


#include <ql/methods/lattices/binomialtree.hpp>
//#include <binomialtree.hpp>
#include <ql/methods/lattices/bsmlattice.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/vanilla/discretizedvanillaoption.hpp>
#include <ql/pricingengines/greeks.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>

namespace QuantLib {

    //! Pricing engine for vanilla options using binomial trees
    /*! \ingroup vanillaengines

        \test the correctness of the returned values is tested by
              checking it against analytic results.

        \todo Greeks are not overly accurate. They could be improved
              by building a tree so that it has three points at the
              current time. The value would be fetched from the middle
              one, while the two side points would be used for
              estimating partial derivatives.
    */
    template <class T>
    class BinomialVanillaEngine_2 : public VanillaOption::engine {
      public:
        BinomialVanillaEngine_2(
             const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
             Size timeSteps)
        : process_(process), timeSteps_(timeSteps) {
            QL_REQUIRE(timeSteps >= 2,
                       "at least 2 time steps required, "
                       << timeSteps << " provided");
            registerWith(process_);
        }
        void calculate() const;
      private:
        boost::shared_ptr<GeneralizedBlackScholesProcess> process_;
        Size timeSteps_;
    };


    // template definitions

    template <class T>
    void BinomialVanillaEngine_2<T>::calculate() const {

        DayCounter rfdc  = process_->riskFreeRate()->dayCounter();
        DayCounter divdc = process_->dividendYield()->dayCounter();
        DayCounter voldc = process_->blackVolatility()->dayCounter();
        Calendar volcal = process_->blackVolatility()->calendar();

        Real s0 = process_->stateVariable()->value();
        QL_REQUIRE(s0 > 0.0, "negative or null underlying given");
        Volatility v = process_->blackVolatility()->blackVol(
            arguments_.exercise->lastDate(), s0);
        Date maturityDate = arguments_.exercise->lastDate();
        Rate r = process_->riskFreeRate()->zeroRate(maturityDate,
            rfdc, Continuous, NoFrequency);
        Rate q = process_->dividendYield()->zeroRate(maturityDate,
            divdc, Continuous, NoFrequency);
        Date referenceDate = process_->riskFreeRate()->referenceDate();

        // binomial trees with constant coefficient
        Handle<YieldTermStructure> flatRiskFree(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(referenceDate, r, rfdc)));
        Handle<YieldTermStructure> flatDividends(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(referenceDate, q, divdc)));
        Handle<BlackVolTermStructure> flatVol(
            boost::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(referenceDate, volcal, v, voldc)));

        boost::shared_ptr<PlainVanillaPayoff> payoff =
            boost::dynamic_pointer_cast<PlainVanillaPayoff>(arguments_.payoff);
        QL_REQUIRE(payoff, "non-plain payoff given");

        Time maturity = rfdc.yearFraction(referenceDate, maturityDate);

        boost::shared_ptr<StochasticProcess1D> bs(
                         new GeneralizedBlackScholesProcess(
                                      process_->stateVariable(),
                                      flatDividends, flatRiskFree, flatVol));

        TimeGrid grid(maturity, timeSteps_);

        boost::shared_ptr<T> tree(new T(bs, maturity, timeSteps_,
                                        payoff->strike()));

        boost::shared_ptr<BlackScholesLattice<T> > lattice(
            new BlackScholesLattice<T>(tree, r, maturity, timeSteps_));

        DiscretizedVanillaOption option(arguments_, *process_, grid);

        option.initialize(lattice, maturity);

        // Partial derivatives calculated from various points in the
        // binomial tree 
        // (see J.C.Hull, "Options, Futures and other derivatives", 6th edition, pp 397/398




        // option values (p0) at this point
        option.rollback(grid[0]);
        Array va0(option.values());
        QL_ENSURE(va0.size() == 3, "Expect 3 nodes in grid at second step");
        Real p0u = va0[2]; // up
        Real p0m = va0[1]; // mid
        Real p0d = va0[0]; // down (low)
        Real s0u = lattice->underlying(0, 2); // up price
        Real s0m = lattice->underlying(0, 1); // middle price
        Real s0d = lattice->underlying(0, 0); // down (low) price

        // calculate gamma by taking the first derivate of the two deltas
        /* Real delta0u = (p0u - p0m)/(s0u-s0m);
        Real delta0d = (p0m-p0d)/(s0m-s0d);
        Real gamma0 = (delta0u - delta0d) / ((s0u-s0d)/2);
        Real delta0 = (delta0u + delta0d)/2; */

        Real h2 = s0u - s0m;
        Real h1 = s0m - s0d;
        Real f0 = p0d;
        Real f1 = p0m;
        Real f2 = p0u;

        Real gamma0 = 2*(h2*f0 - (h1+h2)*f1 + h1*f2)/((h1*h2)*(h1+h2));
        Real delta0 = (-h2/(h1*(h1+h2)))*f0 - ((h1-h2)/(h1*h2))*f1 + (h1/(h2*(h1+h2)))*f2;
        /* Real delta0 = ((2*h1 + h2)/(h1*(h1+h2)))*f0 - ((h1+h2)/(h1*h2))*f1 - (h1/(h2*(h1+h2)))*f2; */



        // Rollback to second-last step, and get option values (p1) at
        // this point
        /*
        option.rollback(grid[1]);
        Array va(option.values());
        QL_ENSURE(va.size() == 2, "Expect 2 nodes in grid at first step");
        Real p1u = va[1];
        Real p1d = va[0];
        Real s1u = lattice->underlying(1, 1); // up (high) price
        Real s1d = lattice->underlying(1, 0); // down (low) price

        Real delta = (p1u - p1d) / (s1u - s1d);

        // Finally, rollback to t=0
        option.rollback(0.0);
        Real p0 = option.presentValue();
*/
        // Store results
        results_.value = p0m;
        results_.delta = delta0;
        results_.gamma = gamma0;
        results_.theta = blackScholesTheta(process_,
                                           results_.value,
                                           results_.delta,
                                           results_.gamma);                                      
    }

}


#endif
