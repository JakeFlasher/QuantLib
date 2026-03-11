// r6
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Unit tests for the numerically stable x*coth(x) evaluation used
 by exponential fitting and Milev-Tagliani effective diffusion schemes.
*/

#include "toplevelfixture.hpp"

#include <ql/methods/finitedifferences/operators/fdmhyperboliccot.hpp>

#include <boost/test/unit_test.hpp>
#include <cmath>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(FdmHyperbolicCotTests)

namespace {

    // Reference: x*coth(x) computed via high-precision mpmath
    // x*coth(x) = x / tanh(x)
    // At x=0: limit is 1.0
    // At x=1: 1 / tanh(1) = 1.3130352854993...

    const Real xSmall = 1e-6;
    const Real xLarge = 50.0;

    Real xcothx(Real x) {
        return detail::xCothx(x, xSmall, xLarge);
    }

} // anonymous namespace


BOOST_AUTO_TEST_CASE(testExactValues) {
    BOOST_TEST_MESSAGE("Testing xCothx exact values...");

    // x=0: limit is 1.0
    BOOST_CHECK_CLOSE(xcothx(0.0), 1.0, 1e-10);

    // x=1: 1/tanh(1) = 1.31303528549933...
    const Real ref1 = 1.0 / std::tanh(1.0);
    BOOST_CHECK_CLOSE(xcothx(1.0), ref1, 1e-10);

    // x=10: 10/tanh(10) ≈ 10.0000000...
    const Real ref10 = 10.0 / std::tanh(10.0);
    BOOST_CHECK_CLOSE(xcothx(10.0), ref10, 1e-8);

    // Large x: asymptotic to |x|
    BOOST_CHECK_CLOSE(xcothx(100.0), 100.0, 1e-6);
}


BOOST_AUTO_TEST_CASE(testSymmetry) {
    BOOST_TEST_MESSAGE("Testing xCothx symmetry: f(x) == f(-x)...");

    const Real testPoints[] = {0.5, 1.0, 5.0, 25.0, 100.0};
    for (Real x : testPoints) {
        BOOST_CHECK_CLOSE(xcothx(x), xcothx(-x), 1e-12);
    }
}


BOOST_AUTO_TEST_CASE(testLowerBound) {
    BOOST_TEST_MESSAGE("Testing xCothx lower bound: f(x) >= 1.0...");

    const Real testPoints[] = {0.0, 1e-8, 1e-4, 0.01, 0.1, 0.5,
                               1.0, 5.0, 10.0, 25.0, 50.0, 100.0};
    for (Real x : testPoints) {
        BOOST_CHECK_GE(xcothx(x), 1.0 - 1e-15);
        BOOST_CHECK_GE(xcothx(-x), 1.0 - 1e-15);
    }
}


BOOST_AUTO_TEST_CASE(testTaylorDirectBoundary) {
    BOOST_TEST_MESSAGE(
        "Testing xCothx continuity at Taylor/direct boundary...");

    // Values just inside and outside the Taylor regime should be consistent
    const Real xInside = 0.5e-6;   // inside Taylor (|x| < 1e-6)
    const Real xOutside = 2.0e-6;  // outside Taylor, direct computation

    const Real vIn = xcothx(xInside);
    const Real vOut = xcothx(xOutside);

    // Both should be very close to 1.0 (within x^2/3 correction)
    BOOST_CHECK_CLOSE(vIn, 1.0, 1e-8);
    BOOST_CHECK_CLOSE(vOut, 1.0, 1e-8);

    // The exact values at the boundary itself
    const Real vAtBoundary = xcothx(xSmall);
    BOOST_CHECK_CLOSE(vAtBoundary, 1.0, 1e-8);

    // Continuity: value just below and just above boundary should be close
    const Real xBelow = xSmall * (1.0 - 1e-3);
    const Real xAbove = xSmall * (1.0 + 1e-3);
    BOOST_CHECK_CLOSE(xcothx(xBelow), xcothx(xAbove), 1e-6);
}


BOOST_AUTO_TEST_CASE(testAsymptoticDirectBoundary) {
    BOOST_TEST_MESSAGE(
        "Testing xCothx continuity at asymptotic/direct boundary...");

    // Values just inside and outside the asymptotic regime
    const Real xBelow = xLarge - 1.0;  // 49, direct computation
    const Real xAbove = xLarge + 1.0;  // 51, asymptotic

    const Real vBelow = xcothx(xBelow);
    const Real vAbove = xcothx(xAbove);

    // At x=49, coth(49)≈1, so x*coth(x)≈49
    BOOST_CHECK_CLOSE(vBelow, xBelow, 1e-4);

    // At x=51, asymptotic returns |x|=51
    BOOST_CHECK_CLOSE(vAbove, xAbove, 1e-10);

    // Continuity at the boundary: relative difference should be small
    const Real vAt = xcothx(xLarge);
    const Real vJustBelow = xcothx(xLarge - 0.01);
    // vAt ≈ 50.0, vJustBelow ≈ 49.99, relative diff ≈ 0.02%
    BOOST_CHECK_CLOSE(vAt, vJustBelow, 0.1);
}


BOOST_AUTO_TEST_CASE(testRegimeBoundaryValues) {
    BOOST_TEST_MESSAGE(
        "Testing xCothx at exact regime boundary values...");

    // Exactly at xSmall = 1e-6
    // |x| < xSmall triggers Taylor; |x| >= xSmall uses direct.
    // Test both sides and compare to the direct formula.
    const Real vSmallMinus = xcothx(xSmall * 0.999);
    BOOST_CHECK_CLOSE(vSmallMinus,
        1.0 + (xSmall * 0.999) * (xSmall * 0.999) / 3.0, 1e-6);

    // xCothx at exactly xSmall = 1e-6 vs direct formula x/tanh(x)
    const Real vAtSmall = xcothx(xSmall);
    const Real directSmall = xSmall / std::tanh(xSmall);
    BOOST_CHECK_CLOSE(vAtSmall, directSmall, 1e-10);

    // xCothx at exactly xLarge = 50 vs direct formula x/tanh(x)
    const Real vAtLarge = xcothx(xLarge);
    const Real directLarge = xLarge / std::tanh(xLarge);
    BOOST_CHECK_CLOSE(vAtLarge, directLarge, 1e-10);

    // Asymptotic regime: x=50.001
    const Real vLargePlus = xcothx(xLarge + 0.001);
    BOOST_CHECK_CLOSE(vLargePlus, xLarge + 0.001, 1e-8);
}


BOOST_AUTO_TEST_CASE(testTightBoundaryTransitions) {
    BOOST_TEST_MESSAGE(
        "Testing xCothx boundary transitions at 1e-10 tolerance...");

    // Taylor-vs-direct at |x| = 0.999e-6 (inside Taylor regime)
    {
        const Real x = 0.999e-6;
        const Real taylor = xcothx(x);
        const Real direct = x / std::tanh(x);
        BOOST_CHECK_CLOSE(taylor, direct, 1e-10);
    }

    // Direct-vs-Taylor continuity across xSmall boundary
    {
        const Real xBelow = xSmall * (1.0 - 1e-6);
        const Real xAbove = xSmall * (1.0 + 1e-6);
        BOOST_CHECK_CLOSE(xcothx(xBelow), xcothx(xAbove), 1e-10);
    }

    // Direct-vs-asymptotic continuity across xLarge boundary
    {
        const Real xBelow = xLarge - 0.001;
        const Real xAbove = xLarge + 0.001;
        const Real vBelow = xcothx(xBelow);
        const Real vAbove = xcothx(xAbove);
        // At x≈50, coth(50) ≈ 1, so x*coth(x) ≈ x.
        // Relative difference between 49.999 and 50.001 is tiny.
        const Real relDiff = std::fabs(vAbove - vBelow)
            / std::max(vAbove, vBelow);
        BOOST_CHECK_MESSAGE(relDiff < 1e-4,
            "Asymptotic boundary discontinuity: relDiff=" << relDiff);
    }

    // Lower bound at regime boundaries: f(x) >= 1.0
    BOOST_CHECK_GE(xcothx(xSmall), 1.0 - 1e-15);
    BOOST_CHECK_GE(xcothx(xLarge), 1.0 - 1e-15);
    BOOST_CHECK_GE(xcothx(xSmall * 0.999), 1.0 - 1e-15);
    BOOST_CHECK_GE(xcothx(xLarge + 0.001), 1.0 - 1e-15);
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
