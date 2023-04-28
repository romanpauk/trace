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

TEST(trace_test, test_sleep_100)
{
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

TEST(trace_test, test_overhead)
{
    trace::frame_registry< trace::frame_data >::instance().clear();

    {
        TRACE("A");
        for (size_t i = 0; i < 100; ++i) {
            TRACE("B");
            for (size_t i = 0; i < 100; ++i) {
                TRACE("C");
                for (size_t i = 0; i < 100; ++i) {
                    TRACE("D");
                    std::string xxx("1234567812345678");
                }
            }
        }
    }

    {
        TRACE("A2");
        for (size_t i = 0; i < 100; ++i) {
            for (size_t i = 0; i < 100; ++i) {
                for (size_t i = 0; i < 100; ++i) {
                    std::string xxx("1234567812345678");
                }
            }
        }
    }
    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cout));
}

TEST(time_test, test) {
    using Traits = trace::default_time_traits;
    Traits::value_type begin, end;
    Traits::begin(begin);
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    Traits::end(end);
    std::cerr << Traits::diff(end, begin) << std::endl;
}
