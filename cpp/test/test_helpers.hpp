#pragma once
#include <cmath>
#include <iostream>

#define CHECK_NEAR(a, b, tol)                                                      \
    if (std::abs((a) - (b)) > (tol)) {                                             \
        std::cerr << "  FAIL: " << #a << " expected ~" << (b) << ", got " << (a)  \
                  << "\n";                                                          \
        return false;                                                               \
    }

#define CHECK_EQ(a, b)                                                             \
    if ((a) != (b)) {                                                              \
        std::cerr << "  FAIL: " << #a << " expected " << (b) << ", got " << (a)   \
                  << "\n";                                                          \
        return false;                                                               \
    }

#define CHECK_TRUE(expr)                                              \
    if (!(expr)) {                                                    \
        std::cerr << "  FAIL: expected true: " << #expr << "\n";     \
        return false;                                                 \
    }
