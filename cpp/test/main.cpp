#include <iostream>
#include <cstdlib>

bool run_phase2_tests();
bool run_phase1_tests();
bool run_anticycling_tests();

#define RUN(name)                                                  \
    do {                                                           \
        std::cout << "[" #name "] running...\n";                   \
        if (run_##name##_tests()) {                                \
            std::cout << "[" #name "] PASS\n\n";                   \
        } else {                                                   \
            std::cout << "[" #name "] FAIL\n\n";                   \
            all_passed = false;                                    \
        }                                                          \
    } while (false)

int main()
{
    bool all_passed = true;

    RUN(phase2);
    RUN(phase1);
    RUN(anticycling);

    if (all_passed) {
        std::cout << "All tests PASSED.\n";
        return EXIT_SUCCESS;
    } else {
        std::cout << "Some tests FAILED.\n";
        return EXIT_FAILURE;
    }
}
