//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif

#include <gtest/gtest.h>

#include <trace/trace.h>

#include <iostream>
#include <cmath>

static void fn()
{
    TRACE(__FUNCTION__);
}

TEST(trace_test, test_basic_numbers) {
    std::cerr << "cpu frequency: " << trace::time_traits< trace::rdtscp_counter >::get_cpu_frequency() << std::endl;

    {
        auto overhead = trace::get_counter_overhead< trace::default_time_traits >();
        std::cerr << "default_time_traits::get() cycle overhead (min, avg, max): " << std::get< 0 >(overhead) << ", " << std::get< 1 >(overhead) << ", " << std::get< 2 >(overhead) << std::endl;
    }
    {
        auto overhead = trace::frame_guard< trace::frame_data >::get_overhead();
        std::cerr << "frame_guard::get() cycle overhead (min, avg, max): " << std::get< 0 >(overhead) << ", " << std::get< 1 >(overhead) << ", " << std::get< 2 >(overhead) << std::endl;
    }

    //std::cerr << "TRACE() overhead " << trace::frame_registry< trace::frame_data >::instance().trace_overhead() * trace::default_time_traits::get_frequency() << std::endl;
}

TEST(trace_test, test_stack) {
    {
        TRACE("A");
        fn();
        {
            std::string b("bbbb");
            fn();
            {
                TRACE("F");
                fn();
            }

            TRACE("D");
            for (int i = 0; i < 100; ++i) {
                TRACE("C");
            }

            for(int i = 0; i < 10; ++i) {
                TRACE("B");
                std::string a("aaaa");
                fn();
            }
        }
        
    }
    
    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cout));
}

#if defined(__linux__)
TEST(trace_test, test_sleep_1) {
    trace::frame_registry< trace::frame_data >::instance().clear();

    {
        TRACE("sleep1-total");
        for (size_t i = 0; i < 1000; ++i) {
            TRACE("sleep1");
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }

    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cout));
}
#endif

TEST(trace_test, test_sleep_100) {
    trace::frame_registry< trace::frame_data >::instance().clear();

    {
        TRACE("sleep100-total");
        for (size_t i = 0; i < 10; ++i) {
            TRACE("sleep100");
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cout));
}

TEST(trace_test, test_empty) {
    trace::frame_registry< trace::frame_data >::instance().clear();
    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cout));
}

TEST(time_test, test_time_traits_sleep) {
    using Traits = trace::default_time_traits;
    auto begin = Traits::begin();
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    auto end = Traits::end();
    std::cerr << end - begin << std::endl;
}

void cycle(size_t n, volatile size_t* dummy, size_t loads = 1) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < loads; ++j)
            (*dummy)++;
    }
}

TEST(time_test, test_overhead_base) {
    trace::frame_registry< trace::frame_data >::instance().clear();
    using Traits = trace::default_time_traits;

    {
        TRACE(__FUNCTION__);

        const int N = 1000000;

        {
            for (size_t i = 0; i < N; ++i) {
            }
        }

        {
            volatile size_t dummy;
            for(size_t i = 0; i < N; ++i) {
                cycle(10, &dummy);
            }
        }

        {
            for (size_t i = 0; i < N; ++i) {
                unsigned dummy;
                __rdtscp(&dummy);
            }
        }

        {
            for (size_t i = 0; i < N; ++i) {
                unsigned dummy;
                __rdtscp(&dummy);
                __rdtscp(&dummy);
            }
        }

        {
            for (size_t i = 0; i < N; ++i) {
                unsigned dummy;
                __rdtscp(&dummy);
                __rdtscp(&dummy);
                __rdtscp(&dummy);
            }
        }
    }

    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cout));
}

TEST(time_test, test_overhead) {
    trace::frame_registry< trace::frame_data >::instance().clear();
    using Traits = trace::default_time_traits;

    {
        TRACE(__FUNCTION__);

        const int N = 1000000;

        {
            TRACE("empty");
            for (size_t i = 0; i < N; ++i) {
                TRACE("emptyloop");
            }
        }

        {
            TRACE("cycle10");
            volatile size_t dummy;
            for(size_t i = 0; i < N; ++i) {
                TRACE("cycleloop");
                cycle(10, &dummy);
            }
        }

        {
            TRACE("rdtsc1");
            for(size_t i = 0; i < N; ++i) {
                TRACE("rdtsc1loop");
                unsigned dummy;
                __rdtscp(&dummy);
            }
        }

        {
            TRACE("rdtsc2");
            for (size_t i = 0; i < N; ++i) {
                TRACE("rdtsc2loop");
                unsigned dummy;
                __rdtscp(&dummy);
                __rdtscp(&dummy);
            }

            volatile int a = 1;
        }

        {
            TRACE("rdtsc3");
            for (size_t i = 0; i < N; ++i) {
                TRACE("rdtsc3loop");
                unsigned dummy;
                __rdtscp(&dummy);
                __rdtscp(&dummy);
                __rdtscp(&dummy);
            }
        }
    }

    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cout));
}

TEST(time_test, test_differential_measure_base) {
    using Traits = trace::default_time_traits;

    const int N = 10000;

    volatile uint64_t dummy;
    int64_t begin, end, t1, t2, e = 0, f = 0;
    size_t n1 = 1, n2 = 10;

    for (size_t i = 0; i < N; ++i) {
        begin = Traits::begin();
        cycle(1, &dummy);
        end = Traits::end();
        t1 = end - begin;

        begin = Traits::begin();
        for (size_t i = 0; i < n2; ++i) {
            cycle(1, &dummy);
        }
        end = Traits::end();
        t2 = end - begin;

        e += double(n1*t2-n2*t1)/(n1-n2);
        f += double(t2-t1)/(n2-n1);
    }

    std::cerr << "f: " << f / N << ", e: " << e / N << std::endl;
}

TEST(time_test, test_differential_measure_lambda) {
    using Traits = trace::default_time_traits;
    volatile uint64_t dummy;
    auto [f, e] = trace::differential_measure< Traits >([&] {
        cycle(1, &dummy);
    });
    std::cerr << "f: " << f << ", e: " << e << std::endl;
}

template< typename Traits, typename Fn > void measure(Fn&& fn) {
    trace::frame_registry< trace::frame_data >::instance().clear();
    
    fn();

    auto [f, e] = trace::differential_measure< Traits >(fn);
    std::cerr << typeid(Traits).name() << ": differential_measure: " << f + e << ", error=" << e << ", result=" << f << std::endl;

    for (size_t i = 0; i < 100000; ++i)
    {
        TRACE("<dummy>");
        fn();
    }

    // TODO: to extract frame data, we need to call merge() here, than sort(), than we need to calculate them.
    // So, 1. there is missing model of the data, 2. there is a mess in the lifetime of the data.
    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cerr));
    trace::frame_registry< trace::frame_data >::instance().clear();
}

TEST(time_test, test_cycle_timings) {
    const size_t Loops = 32;
    size_t dummy;
    for (size_t i = 0; i < Loops; ++i) {
        std::cerr << "CYCLE " << i << std::endl;
        measure< trace::rdtscp_counter >([&] { cycle(i, &dummy); });
    }
}

// How to Benchmark Code Execution Times on Intel IA-32 and IA-64 Instruction Set Architectures
// https://www.intel.com/content/dam/www/public/us/en/documents/white-papers/ia-32-ia-64-benchmark-code-execution-paper.pdf
template< typename Counter, size_t Loops = 1000, typename Begin, typename End, size_t N = 100, typename Fn > void test_counter(Fn fn) {
    std::vector< std::array< uint64_t, N > > times(Loops);
    std::array< uint64_t, N > variances;

    auto overhead = trace::get_counter_overhead< Counter >();

    for (size_t n = 0; n < Loops; ++n) {
        for (size_t i = 0; i < N; ++i) {
            auto start = Counter::begin();
            fn(n);
            auto end = Counter::end();

            times[n][i] = end - start;
        }

        std::sort(times[n].begin(), times[n].end(), [](auto lhs, auto rhs){ return lhs < rhs; });
        
        uint64_t min_time = -1, max_time = 0, max_dev;
        double mean = 0, sd = 0, var = 0;
        for (size_t i = Begin::value() * N; i < End::value() * N; ++i) {
            min_time = std::min(times[n][i], min_time);
            max_time = std::max(times[n][i], max_time);
            mean += times[n][i];
        }
        mean /= N*(End::value() - Begin::value());
        max_dev = max_time - min_time;
        var = trace::variance(times[n].data(), Begin::value()*N, End::value()*N);
        sd = sqrt(var);
        variances[n] = var;

        auto [f, e] = trace::differential_measure< Counter, Begin, End >([&](){ fn(n); });

        std::cerr << typeid(Counter).name() << ": N=" << n << ": mean=" << mean << ", var=" << var << ", sd=" << sd << ", min=" << min_time << ", max=" << max_time
            << ", F=" << f << ", e=" << e << ", [F]=" << trace::compensate(overhead, mean, 1) << ", diff=" << trace::compensate(overhead, mean, 1) - f << "]" << std::endl;
    }
}

TEST(time_test, test_counter_stats) {
    const size_t Loops = 32;
    volatile size_t dummy;
    test_counter< trace::rdtscp_counter, Loops, trace::ratio<0>, trace::ratio<3,4> >([&](size_t n){ cycle(n, &dummy); });
}

TEST(time_test, test_histogram) {
    trace::histogram histogram;
    ASSERT_EQ(histogram.mean(), 0);
    histogram.add(10, 10);
    ASSERT_EQ(histogram.mean(), 10);
}
