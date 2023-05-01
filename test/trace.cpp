//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#include <gtest/gtest.h>

#include <trace/trace.h>

#include <iostream>

static void fn()
{
    TRACE(__FUNCTION__);
}

TEST(trace_test, test_basic_numbers) {
    std::cerr << "frequency: " << trace::default_time_traits::get_frequency() << std::endl;

    {
        auto overhead = trace::default_time_traits::get_overhead();
        std::cerr << "default_time_traits::get() cycle overhead (min, avg, max): " << std::get< 0 >(overhead) << ", " << std::get< 1 >(overhead) << ", " << std::get< 2 >(overhead) << std::endl;
    }
    {
        auto overhead = trace::frame_guard< trace::frame_data >::get_overhead();
        std::cerr << "frame_guard::get() cycle overhead (min, avg, max): " << std::get< 0 >(overhead) << ", " << std::get< 1 >(overhead) << ", " << std::get< 2 >(overhead) << std::endl;
    }

    std::cerr << "TRACE() overhead " << trace::frame_registry< trace::frame_data >::instance().trace_overhead() * trace::default_time_traits::get_frequency() << std::endl;
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

TEST(time_test, test_time_traits_sleep) {
    using Traits = trace::default_time_traits;
    auto begin = Traits::get();
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    auto end = Traits::get();
    std::cerr << Traits::diff(end, begin) << std::endl;
}

TEST(time_test, test_overhead_base) {
    trace::frame_registry< trace::frame_data >::instance().clear();
    using Traits = trace::default_time_traits;

    auto begin = Traits::get();
    {
        TRACE(__FUNCTION__);

        const int N = 1000000;

        {
            for (size_t i = 0; i < N; ++i) {
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

    auto begin = Traits::get();
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
    int64_t begin, end, t1, t2, n1 = 1, n2 = 10, e = 0, f = 0;

    for (size_t i = 0; i < N; ++i) {
        begin = Traits::get();
        dummy = Traits::get();
        end = Traits::get();
        t1 = end - begin;

        begin = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        dummy = Traits::get();
        end = Traits::get();
        t2 = end - begin;

        e += double(n1*t2-n2*t1)/(n1-n2);
        f += double(t2-t1)/(n2-n1);
    }

    std::cerr << "f: " << f / N << ", e: " << e / N << std::endl;
}

TEST(time_test, test_differential_measure_lambda) {
    using Traits = trace::default_time_traits;
    auto [f, e] = trace::differential_measure< Traits >([] { volatile uint64_t dummy = Traits::get(); });
    std::cerr << "f: " << f << ", e: " << e << std::endl;
}

template< typename Fn > void measure(const char* name, Fn&& fn) {
    trace::frame_registry< trace::frame_data >::instance().clear();
    using Traits = trace::default_time_traits;

    fn();

    auto [f, e] = trace::differential_measure< Traits >(fn);
    std::cerr << name << ": differential measure: " << f << std::endl;

    for (size_t i = 0; i < 100000; ++i)
    {
        TRACE(name);
        fn();
    }

    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cerr));
}

TEST(time_test, test_timings) {
    measure("_mm_pause()", [] { _mm_pause(); });
    measure("_mm_lfence()", [] { _mm_lfence(); });
    measure("_mm_mfence()", [] { _mm_mfence(); });
    measure("_mm_sfence()", [] { _mm_sfence(); });
    //measure("__rdtsc()", [] { volatile uint64_t dummy = __rdtsc(); });
    //measure("__rdtscp()", [] { unsigned id; volatile uint64_t dummy = __rdtscp(&id); });
    measure("default_time_traits::get()", [] { volatile uint64_t dummy = trace::default_time_traits::get(); });

#if defined(_WIN32)
    measure("qpc_time_traits::get()", [] { volatile auto value = trace::qpc_time_traits::get(); });
#endif

#if defined(__linux__)
    measure("clock_gettime_time_traits< CLOCK_MONOTONIC >::get()", [] { volatile auto value = trace::clock_gettime_time_traits< CLOCK_MONOTONIC >::get(); });
    measure("clock_gettime_time_traits< CLOCK_THREAD_CPUTIME_ID >::get()", [] { volatile auto value = trace::clock_gettime_time_traits< CLOCK_THREAD_CPUTIME_ID >::get(); });
#endif
    
}
