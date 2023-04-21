//
// This file is part of trace project <https://github.com/romanpauk/trace>
//
// See LICENSE for license and copyright information
// SPDX-License-Identifier: AGPL-3.0-or-later
//

#include <gtest/gtest.h>

#include <trace/trace.h>

#include <iostream>

TEST(trace_test, test) {
    {
        TRACE("A");
        {
            TRACE("B");
            for(int i = 0; i < 10; ++i) {
                TRACE("C"); 
            }
            TRACE("D");
        }
        TRACE("E");
    }

    trace::frame_registry< trace::frame_data >::instance().for_each(trace::stream_dumper(std::cout));
}
