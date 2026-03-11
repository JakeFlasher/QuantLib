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
    // Since |x| < xSmall triggers Taylor, 1e-6 does NOT trigger (it's ==, not <)
    // but 0.999e-6 does. Test both sides.
    const Real vSmallMinus = xcothx(xSmall * 0.999);
    BOOST_CHECK_CLOSE(vSmallMinus, 1.0 + (xSmall * 0.999) * (xSmall * 0.999) / 3.0, 1e-6);

    // Exactly at xLarge = 50
    // Since |x| > xLarge triggers asymptotic, x=50 does NOT trigger.
    // But x=50.001 does.
    const Real vLargePlus = xcothx(xLarge + 0.001);
    BOOST_CHECK_CLOSE(vLargePlus, xLarge + 0.001, 1e-8);
}


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
